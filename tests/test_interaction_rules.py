"""Cross-mode regressions: configurable cutoffs, boundaries and missing geometry."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from pydantic import ValidationError
from rdkit import Chem

from interactions_search.cli import _build_cfg
from interactions_search.config import Angles, Distances, HotspotPocket
from interactions_search.contacts import search_hydrophobic, search_pi_cation, search_salt_bridges
from interactions_search.geometry import angle_three_points
from interactions_search.hotspot_pocket import _annotate_grid
from interactions_search.interaction_rules import aromatic_status, hbond_status
from interactions_search.ligand_chemistry import load_reference
from interactions_search.pipeline import analyze_pair, analyze_probe
from interactions_search.probe import probe_interactions

ROOT = Path(__file__).resolve().parents[1]
SITE_COLUMNS = ['Pos', 'Residue', 'Atom', 'X', 'Y', 'Z']
POINT_COLUMNS = ['Type', *SITE_COLUMNS]


def make_cfg(tmp_path, distances=None, angles=None):
    data = yaml.safe_load((ROOT / 'Interacciones_variables.yml').read_text())
    data['distancias'].update(distances or {})
    data['angulos'].update(angles or {})
    for option in ('ligand_plot', 'vmd_output', 'cumulative_output', 'volume_plot', 'bias'):
        data['options'][option] = 'No'
    path = tmp_path / 'config.yml'
    path.write_text(yaml.safe_dump(data))
    return _build_cfg(path)


def pdb_line(serial, name, xyz, *, residue='ALA', record='ATOM'):
    x, y, z = xyz
    return (f'{record:<6}{serial:5d} {name:>4} {residue:>3} A{1:4d}    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {name[0]:>2}\n')


@pytest.mark.parametrize('distance,cutoff,search,missing,status,in_threshold', [
    (4.5, 5.0, 4.0, False, 'Yes', True),  # formerly discarded at 4 Å
    (4.5, 5.0, 1.0, False, 'Yes', True),  # candidate radius cannot truncate final cutoff
    (3.5, 3.2, 4.0, False, 'No', False),  # all retains rejected candidates
    (5.0, 5.0, 6.0, False, 'No', False),  # strict distance boundary
    (2.5, 3.2, 4.0, True, 'No', True),   # missing antecedent must not validate
])
def test_pair_probe_hbond_and_threshold_csv(tmp_path, monkeypatch, distance, cutoff,
                                          search, missing, status, in_threshold):
    cfg = make_cfg(tmp_path, {'Distances_Hidrogen_Bonds': cutoff,
                              'Hydrogen_Bond_Search_Distance': search,
                              'Distances_Aromatic': 2.0})
    receptor = tmp_path / 'rec.pdb'
    atoms = [('O', (distance, 0, 0)), ('CA', (distance+2.5, 1, 0)),
             ('N', (distance+3.5, 1, 0)), ('CB', (distance+2.5, 1, 1))]
    if not missing:
        atoms.append(('C', (distance+1, 0, 0)))
    receptor.write_text(''.join(pdb_line(i, name, xyz)
                                for i, (name, xyz) in enumerate(atoms, 1)) + 'END\n')
    ligand = tmp_path / 'lig.pdb'
    ligand.write_text(pdb_line(1, 'C1', (-1.45, 0, 0), residue='LIG', record='HETATM')
                      + pdb_line(2, 'N1', (0, 0, 0), residue='LIG', record='HETATM')
                      + 'CONECT    1    2\nCONECT    2    1\nEND\n')
    monkeypatch.chdir(tmp_path)
    analyze_pair(str(receptor), str(ligand), 'A', cfg, ligand_reference=load_reference(smiles='CN'))
    analyze_probe(str(receptor), 'A', [(0, 0, 0, ['donor'])], cfg, label='check')
    pair = pd.read_csv(tmp_path / 'rec_lig/Interaction_rec_lig_all.csv')
    probe = pd.read_csv(tmp_path / 'rec_probe_check/Probe_rec_all.csv')
    p = pair[(pair.Type == 'donor') & (pair.Atom == 'O')].iloc[0]
    q = probe[(probe.Type == 'donor') & (probe.Atom == 'O')].iloc[0]
    assert p.Interaction == q.Interaction == status
    assert p.Dist == q.Dist == pytest.approx(distance)
    if missing:
        assert np.isnan(p.Angle) and np.isnan(q.Angle)
    else:
        assert p.Angle == q.Angle == pytest.approx(180)
    threshold = pd.read_csv(tmp_path / 'rec_lig/Interaction_rec_lig_threshold.csv')
    assert bool(((threshold.Type == 'donor') & (threshold.Atom == 'O')).any()) == in_threshold


@pytest.mark.parametrize('kind,role,residue,atom,smiles,key', [
    ('hydrophobic', 'hydrophobic', 'ALA', 'CB', 'C', 'Distances_Hidrofobica'),
    ('salt_bridge', 'cation', 'ASP', 'OD1', 'C[NH3+]', 'Distances_Salt_Bridge'),
    ('salt_bridge', 'anion', 'LYS', 'NZ', '[O-]', 'Distances_Salt_Bridge'),
    ('pi_cation', 'aromatic', 'LYS', 'NZ', None, 'Distances_Pi_Cation'),
])
@pytest.mark.parametrize('distance,expected', [(5.9999, 1), (6.0, 0), (6.0001, 0)])
def test_distance_only_contacts_match_probe(tmp_path, kind, role, residue, atom,
                                           smiles, key, distance, expected):
    cfg = make_cfg(tmp_path, {key: 6.0})
    site = pd.DataFrame([[1, residue, atom, distance, 0, 0]], columns=SITE_COLUMNS)
    points = pd.DataFrame(columns=POINT_COLUMNS)
    if kind == 'pi_cation':
        ring = pd.DataFrame({'Caso': ['aromatic 1'] * 3,
                             'Coord X': [-1, 1, 0], 'Coord Y': [0, 0, 0], 'Coord Z': [0, 0, 0]})
        pair = search_pi_cation(site, ring, cfg[key])
    else:
        mol = Chem.MolFromSmiles(smiles)
        coords = [(i+1, a.GetSymbol()+str(i), 'LIG', 'A', 1, 0, 0, 0)
                  for i, a in enumerate(mol.GetAtoms())]
        if kind == 'hydrophobic':
            pair = search_hydrophobic(mol, coords, site, cfg['Distancia_Hidrofobica'])
        else:
            pair = search_salt_bridges(mol, coords, site, cfg[key])
    probe = probe_interactions([(0, 0, 0, [role])], site, points, cfg)
    probe = probe[probe.Type == kind]
    assert len(pair) == len(probe) == expected
    if expected:
        assert pair.Dist.iloc[0] == probe.Dist.iloc[0] == 6.0  # validated before rounding
        assert pair.Interaction.iloc[0] == probe.Interaction.iloc[0] == 'Yes'


@pytest.mark.parametrize('angle,expected', [(100, 'No'), (100.01, 'Yes'), (180, 'Yes'),
                                           (180.01, 'No'), (np.nan, 'No')])
def test_hbond_angle_boundaries(angle, expected):
    assert hbond_status(3, angle, {}) == expected


@pytest.mark.parametrize('angle,expected', [(19, 'Yes'), (20, 'No'), (45, 'No'),
                                           (70, 'No'), (71, 'Yes'), (np.nan, 'No')])
def test_custom_aromatic_angles(tmp_path, angle, expected):
    cfg = make_cfg(tmp_path, angles={'Aromatic_Parallel_Max': 20, 'Aromatic_TShaped_Min': 70})
    assert aromatic_status(5, angle, cfg) == expected
    assert aromatic_status(5, angle, cfg, angle_required=False) == 'Yes'


def test_missing_probe_geometry_is_explicit(tmp_path):
    cfg = make_cfg(tmp_path)
    site = pd.DataFrame([[1, 'ALA', 'N', 2.5, 0, 0]], columns=SITE_COLUMNS)
    rp = pd.DataFrame([['Dador', 1, 'ALA', 'N', 2.5, 0, 0]], columns=POINT_COLUMNS)
    result = probe_interactions([(0, 0, 0, ['acceptor'])], site, rp, cfg)
    assert result[result.Type == 'acceptor'].Interaction.tolist() == ['Yes']
    # Explicit H but no heavy parent: geometry is missing, not inapplicable.
    site['Atom'] = 'H'
    rp['Atom'] = 'H'
    result = probe_interactions([(0, 0, 0, ['acceptor'])], site, rp, cfg)
    assert result[result.Type == 'acceptor'].Interaction.tolist() == ['No']


def test_probe_acceptor_beyond_old_prefilter(tmp_path):
    cfg = make_cfg(tmp_path, {'Distances_Hidrogen_Bonds': 5.0})
    site = pd.DataFrame([[1, 'ALA', 'N', 4.5, 0, 0]], columns=SITE_COLUMNS)
    rp = pd.DataFrame([['Dador', 1, 'ALA', 'N', 4.5, 0, 0]], columns=POINT_COLUMNS)
    result = probe_interactions([(0, 0, 0, ['acceptor'])], site, rp, cfg)
    assert result.Interaction.tolist() == ['Yes']


def test_probe_clash_config(tmp_path):
    site = pd.DataFrame([[1, 'ALA', 'CA', 3, 0, 0]], columns=SITE_COLUMNS)
    points = pd.DataFrame(columns=POINT_COLUMNS)
    for cutoff, expected in ((2.5, 0), (3.0, 0), (3.1, 1)):
        cfg = make_cfg(tmp_path, {'Probe_Clash_Distance': cutoff})
        result = probe_interactions([(0, 0, 0, ['donor'])], site, points, cfg)
        assert len(result[result.Type == 'clash']) == expected


def test_hotspot_environment_uses_salt_cutoff_and_strict_boundary(tmp_path):
    cfg = make_cfg(tmp_path, {'Distances_Salt_Bridge': 6.0})
    site = pd.DataFrame([[1, 'ASP', 'OD1', 5.5, 0, 0], [2, 'ASP', 'OD2', 6.0, 0, 0]],
                        columns=SITE_COLUMNS)
    result = _annotate_grid(np.zeros((1, 3)), {}, pd.DataFrame(columns=POINT_COLUMNS),
                            site, cfg, HotspotPocket())
    assert result.N_Rec_Anions.iloc[0] == 1


@pytest.mark.parametrize('field', ['Hydrogen_Bond_Search_Distance', 'Distances_Salt_Bridge',
                                  'Distances_Pi_Cation', 'Probe_Clash_Distance'])
@pytest.mark.parametrize('value', [0, -1, float('inf'), float('nan')])
def test_invalid_new_distance_rejected(field, value):
    with pytest.raises(ValidationError):
        Distances(**{field: value})


def test_invalid_angle_order_rejected():
    with pytest.raises(ValidationError):
        Angles(Angle_Hidrogen_Bonds_Min=180, Angle_Hidrogen_Bonds_Max=100)
    with pytest.raises(ValidationError):
        Angles(Aromatic_Parallel_Max=80, Aromatic_TShaped_Min=60)


def test_legacy_yaml_defaults(tmp_path):
    path = tmp_path / 'old.yml'
    path.write_text('acceptors: {}\ndonors: {}\nacceptors_antecedent: {}\n')
    cfg = _build_cfg(path)
    assert cfg['Distances_Salt_Bridge'] == 4
    assert cfg['Distances_Pi_Cation'] == 5
    assert cfg['Probe_Clash_Distance'] == 2.5
    assert cfg['Aromatic_Parallel_Max'] == 30
    assert cfg['Aromatic_TShaped_Min'] == 60


def test_degenerate_angle_is_not_validated():
    angle = angle_three_points([0, 0, 0], [0, 0, 0], [1, 0, 0])
    assert np.isnan(angle)
    assert hbond_status(1, angle, {}) == 'No'


def write_ligand(path, smiles, coordinates):
    molecule = Chem.MolFromSmiles(smiles)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for i, xyz in enumerate(coordinates):
        conformer.SetAtomPosition(i, tuple(float(v) for v in xyz))
    molecule.AddConformer(conformer)
    Chem.MolToPDBFile(molecule, str(path))


def hexagon():
    theta = np.arange(6) * np.pi / 3
    return np.column_stack([1.4*np.cos(theta), 1.4*np.sin(theta), np.zeros(6)])


@pytest.mark.parametrize('kind', ['salt_bridge', 'pi_cation'])
def test_pair_extended_charged_cutoff_survives_threshold_csv(tmp_path, monkeypatch, kind):
    cfg = make_cfg(tmp_path, {'Distances_Salt_Bridge': 6, 'Distances_Pi_Cation': 6,
                              'Distances_Aromatic': 2})
    if kind == 'salt_bridge':
        smiles, coordinates = 'C[NH3+]', [(-1.45, 0, 0), (0, 0, 0)]
        residue, atom = 'ASP', 'OD1'
    else:
        smiles, coordinates = 'c1ccccc1', hexagon()
        residue, atom = 'LYS', 'NZ'
    ligand, receptor = tmp_path / 'lig.pdb', tmp_path / 'rec.pdb'
    write_ligand(ligand, smiles, coordinates)
    receptor.write_text(pdb_line(1, atom, (0, 0, 5.5), residue=residue) + 'END\n')
    monkeypatch.chdir(tmp_path)
    analyze_pair(str(receptor), str(ligand), 'A', cfg,
                 ligand_reference=load_reference(smiles=smiles))
    for suffix in ('all', 'threshold', 'true'):
        result = pd.read_csv(tmp_path / f'rec_lig/Interaction_rec_lig_{suffix}.csv')
        contacts = result[result.Type == kind]
        assert len(contacts) == 1
        assert contacts.Dist.iloc[0] == pytest.approx(5.5)


@pytest.mark.parametrize('parallel,expected', [(30, 'No'), (50, 'Yes')])
def test_pair_aromatic_angle_uses_yaml(tmp_path, monkeypatch, parallel, expected):
    cfg = make_cfg(tmp_path, angles={'Aromatic_Parallel_Max': parallel})
    ligand, receptor = tmp_path / 'lig.pdb', tmp_path / 'rec.pdb'
    write_ligand(ligand, 'c1ccccc1', hexagon())
    # A receptor ring tilted 45 degrees relative to the ligand plane.
    xyz = hexagon()
    xyz[:, 2] = xyz[:, 1] / np.sqrt(2)
    xyz[:, 1] /= np.sqrt(2)
    xyz[:, 2] += 4.5
    names = ('CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2')
    receptor.write_text(''.join(pdb_line(i, name, coord, residue='PHE')
                                for i, (name, coord) in enumerate(zip(names, xyz), 1)) + 'END\n')
    monkeypatch.chdir(tmp_path)
    analyze_pair(str(receptor), str(ligand), 'A', cfg)
    result = pd.read_csv(tmp_path / 'rec_lig/Interaction_rec_lig_all.csv')
    contact = result[result.Type == 'aromatic'].iloc[0]
    assert contact.Angle == pytest.approx(45, abs=0.1)
    assert contact.Interaction == expected
