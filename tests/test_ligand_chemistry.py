import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from interactions_search.ligand_chemistry import apply_reference, load_reference
from interactions_search.ligand_hotpoints import search_hot_points, search_rings


def pdb_molecule(smiles, explicit_h=False):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    if not explicit_h:
        mol = Chem.RemoveHs(mol)
    mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    return Chem.MolFromPDBBlock(Chem.MolToPDBBlock(mol), removeHs=False, sanitize=False)


@pytest.mark.parametrize('smiles', ['CC(=O)N', 'CC(=O)[O-]', 'C[NH3+]', 'c1cc[nH]c1'])
@pytest.mark.parametrize('explicit_h', [False, True])
def test_reference_preserves_pose_and_atom_identity(smiles, explicit_h):
    reference = load_reference(smiles=smiles)
    pdb = pdb_molecule(smiles, explicit_h)
    result = apply_reference(pdb, reference)
    assert Chem.MolToSmiles(Chem.RemoveHs(result)) == Chem.MolToSmiles(reference)
    np.testing.assert_array_equal(result.GetConformer().GetPositions(),
                                  pdb.GetConformer().GetPositions())
    assert [(a.GetSymbol(), a.GetPDBResidueInfo().GetName(),
             a.GetPDBResidueInfo().GetSerialNumber()) for a in result.GetAtoms()] == [
        (a.GetSymbol(), a.GetPDBResidueInfo().GetName(),
         a.GetPDBResidueInfo().GetSerialNumber()) for a in pdb.GetAtoms()]


def test_repair_single_bonds_and_charges():
    pdb = pdb_molecule('CC(=O)[O-]')
    for bond in pdb.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    for atom in pdb.GetAtoms():
        atom.SetFormalCharge(0)
    result = apply_reference(pdb, load_reference(smiles='CC(=O)[O-]'))
    assert Chem.MolToSmiles(result) == 'CC(=O)[O-]'


@pytest.mark.parametrize('smiles,acceptors,donors', [('CC(=O)N', 1, 1),
                                                   ('C[NH3+]', 0, 1),
                                                   ('c1cc[nH]c1', 0, 1)])
def test_reference_features_do_not_use_pdb_inference(smiles, acceptors, donors, monkeypatch):
    def unexpected(*args):
        pytest.fail('OpenBabel must not override explicit reference chemistry')
    monkeypatch.setattr('interactions_search.ligand_hotpoints._openbabel_h_flags', unexpected)
    mol = apply_reference(pdb_molecule(smiles), load_reference(smiles=smiles))
    acc, don = search_hot_points('unused.pdb', mol, [], 'No', '.', reference_chemistry=True)
    assert len(acc) == acceptors
    assert len(don) == donors


def ring_coords(mol):
    n = mol.GetNumAtoms()
    return [(i, f'C{i}', 'LIG', 'A', 1, np.cos(i * 2*np.pi/n),
             np.sin(i * 2*np.pi/n), 0.0) for i in range(n)]


def test_five_member_ring_and_legacy_filter():
    mol = Chem.MolFromSmiles('c1cc[nH]c1')
    coords = ring_coords(mol)
    assert len(search_rings(mol, coords, 4, 0.15)[0]) == 1
    assert search_rings(mol, coords, 5, 0.15)[0] == []
    assert len(search_rings(mol, coords, 4, 0.15, reference_chemistry=True)[0]) == 1


def test_reference_rejects_planar_nonaromatic_ring():
    mol = Chem.MolFromSmiles('C1CCCCC1')
    coords = ring_coords(mol)
    assert len(search_rings(mol, coords, 4, 0.15)[0]) == 1
    assert search_rings(mol, coords, 4, 0.15, reference_chemistry=True)[0] == []


def test_sdf_reference(tmp_path):
    path = tmp_path / 'ref.sdf'
    with Chem.SDWriter(str(path)) as writer:
        writer.write(Chem.AddHs(Chem.MolFromSmiles('CC(=O)N')))
    assert Chem.MolToSmiles(load_reference(sdf=path)) == 'CC(N)=O'
    with Chem.SDWriter(str(path)) as writer:
        writer.write(Chem.MolFromSmiles('CC'))
        writer.write(Chem.MolFromSmiles('CO'))
    with pytest.raises(ValueError, match='exactamente una'):
        load_reference(sdf=path)


def test_invalid_references():
    assert load_reference() is None
    with pytest.raises(ValueError):
        load_reference(smiles='not a smiles')
    with pytest.raises(ValueError):
        load_reference(smiles='')
    with pytest.raises(ValueError):
        apply_reference(pdb_molecule('CCC'), load_reference(smiles='CCO'))
    with pytest.raises(ValueError):
        apply_reference(pdb_molecule('CC'), load_reference(smiles='CCO'))


def test_incompatible_explicit_hydrogens():
    with pytest.raises(ValueError, match='hidrógenos explícitos'):
        apply_reference(pdb_molecule('CO', explicit_h=True), load_reference(smiles='C[O-]'))


@pytest.mark.parametrize('reference_option', ['--ligand-smiles', '--ligand-sdf'])
def test_cli_reference_end_to_end(tmp_path, monkeypatch, reference_option):
    from interactions_search.cli import main

    receptor = tmp_path / 'receptor.pdb'
    shutil.copy(Path(__file__).parent / 'fixtures/receptor_mini.pdb', receptor)
    ligand = tmp_path / 'ligand.pdb'
    Chem.MolToPDBFile(pdb_molecule('CC(=O)N', explicit_h=True), str(ligand))
    ref = 'CC(=O)N'
    if reference_option == '--ligand-sdf':
        ref = str(tmp_path / 'reference.sdf')
        with Chem.SDWriter(ref) as writer:
            writer.write(Chem.MolFromSmiles('CC(=O)N'))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['interactions-search', '-r', str(receptor),
                                     '-l', str(ligand), '-c', 'A', reference_option, ref])
    before = ligand.read_text()
    main()
    assert ligand.read_text() == before
    assert (tmp_path / 'receptor_ligand/Interaction_receptor_ligand_all.csv').exists()
    assert (tmp_path / 'receptor_ligand/ligand_acceptors.png').exists()
    record = json.loads((tmp_path / 'receptor_ligand/run_metadata.json').read_text())
    assert record['cli_arguments'][reference_option[2:].replace('-', '_')] == ref
    if reference_option == '--ligand-sdf':
        assert record['inputs_before']['source_ligand_sdf']['exists']


@pytest.mark.parametrize('args', [
    ['-r', 'r.pdb', '-l', 'a.pdb', 'b.pdb', '--ligand-smiles', 'CC'],
    ['-r', 'r.pdb', '--probe', '0', '0', '0', '--ligand-smiles', 'CC'],
    ['-r', 'r.pdb', '-l', 'a.pdb', '--ligand-smiles', 'CC', '--ligand-sdf', 'x.sdf'],
])
def test_cli_rejects_ambiguous_reference_usage(args, monkeypatch):
    from interactions_search.cli import main

    monkeypatch.setattr(sys, 'argv', ['interactions-search', '-c', 'A', *args])
    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code == 2


@pytest.mark.parametrize('smiles,n_acceptors,n_donors', [
    ('CN(C)C', 1, 0),        # neutral tertiary amine: no donor H
    ('C[NH+](C)C', 0, 1),    # protonated tertiary amine
    ('C[N+](C)(C)C', 0, 0),  # quaternary ammonium: no lone pair or donor H
    ('CC(=O)O', 1, 1),       # neutral carboxylic acid
    ('CC(=O)[O-]', 2, 0),    # carboxylate
    ('CS(=O)(=O)N', 2, 1),  # sulfonamide N is not an acceptor
    ('CC#N', 1, 0),
    ('c1ccncc1', 1, 0),
    ('c1cc[nH+]cc1', 0, 1),
    ('c1cc[nH]c1', 0, 1),
    ('CC(=O)N(C)C', 1, 0),  # tertiary amide
    ('CS(=O)(=O)C', 2, 0),
])
@pytest.mark.parametrize('explicit_h', [False, True])
def test_protonation_and_functional_group_roles(smiles, n_acceptors, n_donors, explicit_h):
    molecule = apply_reference(pdb_molecule(smiles, explicit_h), load_reference(smiles=smiles))
    acceptors, donors = search_hot_points('unused.pdb', molecule, [], 'No', '.',
                                         reference_chemistry=True)
    assert len(acceptors) == n_acceptors
    assert len(donors) == n_donors
    assert all(molecule.GetAtomWithIdx(i).GetTotalNumHs(includeNeighbors=True) > 0
               for i in donors)


@pytest.mark.parametrize('explicit_h', [False, True])
def test_equivalent_oxygens_preserve_available_pdb_evidence(explicit_h):
    pdb = pdb_molecule('CC(=O)O', explicit_h)
    carbonyl_oxygen = next(a.GetIdx() for a in pdb.GetAtoms()
                           if a.GetAtomicNum() == 8
                           and any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()))
    # Reverse reference atom order to challenge first-match assignment.
    reference = Chem.MolFromSmiles('OC(C)=O')
    result = apply_reference(pdb, reference)
    assert any(b.GetBondType() == Chem.BondType.DOUBLE
               for b in result.GetAtomWithIdx(carbonyl_oxygen).GetBonds())
    assert result.GetAtomWithIdx(carbonyl_oxygen).GetTotalNumHs(includeNeighbors=True) == 0
