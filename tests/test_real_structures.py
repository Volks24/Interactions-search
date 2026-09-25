"""Offline validation against named contacts in two deposited crystal structures.

The oracle reads raw PDB coordinates with the standard library; it does not
call production contact, geometry or feature functions to derive expectations.
"""
import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from interactions_search.cli import _build_cfg
from interactions_search.io_pdb import extract_coords_from_pdb
from interactions_search.ligand_chemistry import apply_reference, load_reference
from interactions_search.ligand_hotpoints import search_hot_points, search_rings
from interactions_search.pipeline import analyze_pair

FIXTURES = Path(__file__).parent / 'fixtures/real'
MANIFEST = json.loads((FIXTURES / 'manifest.json').read_text())
CONTACTS = json.loads((FIXTURES / 'contacts.json').read_text())['contacts']


def raw_atoms(path):
    return {(int(line[22:26]), line[12:16].strip()):
            tuple(float(line[i:i+8]) for i in (30, 38, 46))
            for line in path.read_text().splitlines()
            if line.startswith(('ATOM  ', 'HETATM'))}


def independent_angle(a, b, c):
    u, v = [x-y for x, y in zip(a, b)], [x-y for x, y in zip(c, b)]
    cosine = sum(x*y for x, y in zip(u, v)) / (math.dist(a, b) * math.dist(c, b))
    return math.degrees(math.acos(max(-1, min(1, cosine))))


@pytest.fixture(scope='module')
def cfg():
    config = _build_cfg(None)
    for option in ('ligand_plot', 'vmd_output', 'cumulative_output', 'Volume_Plot', 'Bias'):
        config[option] = 'No'
    return config


def run_pair(work, pdb, cfg, *, reference=True, transform=None, displace=False, drop_atom=None):
    work.mkdir()
    receptor = work / f'{pdb}_receptor.pdb'
    ligand = work / f'{pdb}_ligand.pdb'
    for target in (receptor, ligand):
        shutil.copy(FIXTURES / target.name, target)
        if drop_atom and target == receptor:
            target.write_text(''.join(line for line in target.read_text().splitlines(keepends=True)
                                      if not (line.startswith('ATOM')
                                              and (int(line[22:26]), line[12:16].strip())
                                              == drop_atom)))
        if transform is not None or (displace and target == ligand):
            lines = []
            for line in target.read_text().splitlines(keepends=True):
                if line.startswith(('ATOM  ', 'HETATM')):
                    xyz = [float(line[i:i+8]) for i in (30, 38, 46)]
                    if transform:
                        xyz = transform(xyz)
                    if displace and target == ligand:
                        xyz = [v + 100 for v in xyz]
                    line = line[:30] + ''.join(f'{v:8.3f}' for v in xyz) + line[54:]
                lines.append(line)
            target.write_text(''.join(lines))
    component = 'STI' if pdb == '1iep' else 'BTN'
    ref = (load_reference(smiles=(FIXTURES / f'{component}.smi').read_text().strip())
           if reference else None)
    with pytest.MonkeyPatch.context() as patch:
        patch.chdir(work)
        analyze_pair(str(receptor), str(ligand), 'A', cfg, ligand_reference=ref)
    folder = work / f'{pdb}_receptor_{pdb}_ligand'
    return {suffix: pd.read_csv(folder / f'Interaction_{folder.name}_{suffix}.csv', index_col=0)
            for suffix in ('all', 'threshold', 'true')}


@pytest.fixture(scope='module')
def results(tmp_path_factory, cfg):
    work = tmp_path_factory.mktemp('real')
    return {pdb: run_pair(work / pdb, pdb, cfg) for pdb in ('1iep', '1stp')}


@pytest.mark.parametrize('case', MANIFEST['cases'], ids=lambda case: case['pdb_id'])
def test_fixture_checksums(case):
    for name, expected in case['files'].items():
        assert hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest() == expected


@pytest.mark.parametrize('contact', CONTACTS,
                         ids=lambda c: f"{c['pdb']}-{c['ligand']}-{c['residue']}-{c['atom']}")
def test_real_hbond_contacts(contact, results):
    c = contact
    case = next(case for case in MANIFEST['cases'] if case['pdb_id'].lower() == c['pdb'])
    lig = raw_atoms(FIXTURES / f"{c['pdb']}_ligand.pdb")
    rec = raw_atoms(FIXTURES / f"{c['pdb']}_receptor.pdb")
    a, b = lig[(case['ligand_resid'], c['ligand'])], rec[(c['residue'], c['atom'])]
    distance = math.dist(a, b)
    if c['type'] == 'acceptor':
        angle = independent_angle(b, a, lig[(case['ligand_resid'], c['antecedent'])])
    else:
        angle = independent_angle(a, b, rec[(c['residue'], c['antecedent'])])
    assert distance == pytest.approx(c['distance'], abs=1e-5)
    assert angle == pytest.approx(c['angle'], abs=0.001)
    expected = 'Yes' if distance < 3.2 and 100 < angle <= 180 else 'No'
    assert expected == c['status']

    def select(frame):
        return frame[(frame.Type == c['type']) & (frame.Lig == c['ligand'])
                     & (frame['Pos R'] == c['residue']) & (frame.Atom == c['atom'])]
    row = select(results[c['pdb']]['all'])
    assert len(row) == 1, f"Missing or duplicated reference contact: {c}"
    assert row.Res.iloc[0] == c['resname']
    assert row.Dist.iloc[0] == pytest.approx(distance, abs=1e-5)
    assert row.Angle.iloc[0] == pytest.approx(angle, abs=0.001)
    assert row.Interaction.iloc[0] == expected
    assert len(select(results[c['pdb']]['true'])) == (expected == 'Yes')
    assert len(select(results[c['pdb']]['threshold'])) == (distance < 3.2)


@pytest.mark.parametrize('pdb,component', [('1iep', 'STI'), ('1stp', 'BTN')])
def test_real_chemical_groups(pdb, component):
    path = FIXTURES / f'{pdb}_ligand.pdb'
    coords, _ = extract_coords_from_pdb(path)
    reference = load_reference(smiles=(FIXTURES / f'{component}.smi').read_text().strip())
    molecule = apply_reference(Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False),
                               reference)
    acceptors, donors = search_hot_points(str(path), molecule, coords, 'No', '.',
                                         reference_chemistry=True)
    acceptors = {coords[i][1] for i in acceptors}
    donors = {coords[i][1] for i in donors}
    rings, _ = search_rings(molecule, coords, 4, 0.15, reference_chemistry=True)
    if pdb == '1iep':
        # CCD neutral STI: tertiary piperazine nitrogens have no H.
        assert {'N48', 'N51', 'N3', 'O29'} <= acceptors
        assert donors == {'N13', 'N21'}
        assert 'N21' not in acceptors  # amide N
        assert {frozenset(coords[i][1] for i in ring['Atoms']) for ring in rings} == {
            frozenset(['C1', 'C2', 'N3', 'C4', 'C5', 'C6']),
            frozenset(['C7', 'N8', 'C9', 'N10', 'C11', 'C12']),
            frozenset(['C14', 'C15', 'C16', 'C17', 'C18', 'C19']),
            frozenset(['C23', 'C25', 'C26', 'C27', 'C28', 'C29']),
        }
    else:
        assert 'O3' in acceptors
        assert not {'N1', 'N2'} & acceptors  # ureido nitrogens
        assert {'N1', 'N2'} <= donors
        # The neutral acid has one OH. Its atom-name assignment is ambiguous
        # without bond-order/H evidence; do not claim a unique O11/O12 mapping.
        assert len(donors & {'O11', 'O12'}) == 1
        assert rings == []  # biotin's two five-member rings are not aromatic


@pytest.mark.parametrize('pdb', ['1iep', '1stp'])
def test_rigid_transform_preserves_real_contacts(pdb, results, tmp_path, cfg):
    moved = run_pair(tmp_path / 'moved', pdb, cfg,
                     transform=lambda xyz: [-xyz[1]+20, xyz[0]-10, xyz[2]+30])
    columns = ['Pos R', 'Res', 'Atom', 'Lig', 'Type', 'Dist', 'Angle', 'Interaction']
    keys = ['Pos R', 'Atom', 'Lig', 'Type']
    for suffix in ('all', 'threshold', 'true'):
        before = results[pdb][suffix][columns].sort_values(keys).reset_index(drop=True)
        after = moved[suffix][columns].sort_values(keys).reset_index(drop=True)
        # Receptor ring centroids are rounded to 0.001 Å by the existing code.
        pd.testing.assert_frame_equal(before, after, atol=0.001, rtol=1e-7)


@pytest.mark.parametrize('pdb', ['1iep', '1stp'])
def test_displaced_real_ligand_has_no_contacts(pdb, tmp_path, cfg):
    result = run_pair(tmp_path / 'displaced', pdb, cfg, displace=True)
    assert all(frame.empty for frame in result.values())


@pytest.mark.parametrize('pdb', ['1iep', '1stp'])
def test_pdb_only_real_workflow_remains_available(pdb, tmp_path, cfg):
    result = run_pair(tmp_path / 'pdb_only', pdb, cfg, reference=False)
    assert not result['true'].empty
    assert (result['true'].Interaction == 'Yes').all()


def test_incomplete_real_receptor_ring_is_skipped(tmp_path, cfg):
    result = run_pair(tmp_path / 'incomplete', '1iep', cfg, drop_atom=(317, 'CE1'))['true']
    assert result[(result.Type == 'aromatic') & (result['Pos R'] == 317)].empty
    assert not result[(result.Type == 'aromatic') & (result['Pos R'] == 253)].empty
    assert not result[(result.Type == 'acceptor') & (result.Lig == 'N3')].empty


@pytest.mark.parametrize('residue,label,ring_names,expected', [
    (317, 'aromatic 1 (#6)', ['C1', 'C2', 'N3', 'C4', 'C5', 'C6'], 'Yes'),
    (253, 'aromatic 2 (#6)', ['C7', 'N8', 'C9', 'N10', 'C11', 'C12'], 'Yes'),
    (382, 'aromatic 2 (#6)', ['C7', 'N8', 'C9', 'N10', 'C11', 'C12'], 'No'),
])
def test_real_aromatic_contacts(residue, label, ring_names, expected, results):
    ligand = raw_atoms(FIXTURES / '1iep_ligand.pdb')
    receptor = raw_atoms(FIXTURES / '1iep_receptor.pdb')
    lig_xyz = np.array([ligand[201, name] for name in ring_names])
    names = ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')
    rec_xyz = np.array([receptor[residue, name] for name in names])
    lc, rc = lig_xyz.mean(axis=0), rec_xyz.mean(axis=0)
    ln = np.linalg.svd(lig_xyz-lc)[2][-1]
    rn = np.linalg.svd(rec_xyz-rc)[2][-1]
    # Independent all-atom plane fit, compared with the detector's 3-atom planes.
    angle = math.degrees(math.acos(min(1, abs(float(np.dot(ln, rn))))))
    status = 'Yes' if math.dist(lc, rc) < 5.5 and (angle < 30 or angle > 60) else 'No'
    assert status == expected
    frame = results['1iep']['all']
    row = frame[(frame.Type == 'aromatic') & (frame.Lig == label) & (frame['Pos R'] == residue)]
    assert len(row) == 1
    assert row.Interaction.iloc[0] == expected
    # The historical receptor center uses four ring atoms, not all six, and
    # rounds at 0.001 Å. Allow 0.01 Å against the independent all-atom centroid.
    assert row.Dist.iloc[0] == pytest.approx(math.dist(lc, rc), abs=0.01)
    assert row.Angle.iloc[0] == pytest.approx(angle, abs=2.0)


@pytest.mark.parametrize('pdb,residue,ligand_atom,receptor_atom', [
    ('1iep', 253, 'C6', 'CE1'), ('1stp', 110, 'C9', 'CD2'),
])
def test_real_hydrophobic_contacts(pdb, residue, ligand_atom, receptor_atom, results):
    lig = raw_atoms(FIXTURES / f'{pdb}_ligand.pdb')
    rec = raw_atoms(FIXTURES / f'{pdb}_receptor.pdb')
    lig_resid = 201 if pdb == '1iep' else 300
    distance = math.dist(lig[lig_resid, ligand_atom], rec[residue, receptor_atom])
    assert distance < 4
    frame = results[pdb]['true']
    row = frame[(frame.Type == 'hydrophobic') & (frame.Lig == ligand_atom)
                & (frame['Pos R'] == residue) & (frame.Atom == receptor_atom)]
    assert len(row) == 1
    assert row.Dist.iloc[0] == pytest.approx(distance, abs=0.0005)
