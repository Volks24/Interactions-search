import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from interactions_search.cli import _build_cfg
from interactions_search.config import load_config
from interactions_search.interaction_rules import REASON_CODES, evaluate_aromatic, evaluate_hbond
from interactions_search.ligand_chemistry import load_reference
from interactions_search.pipeline import analyze_pair, analyze_probe, analyze_site_bias
from interactions_search.probe import probe_interactions

FIXTURES = Path(__file__).parent / 'fixtures'


@pytest.fixture
def cfg():
    config = _build_cfg(None)
    for option in ('ligand_plot', 'vmd_output', 'cumulative_output', 'Volume_Plot', 'Bias'):
        config[option] = 'No'
    return config


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    receptor, ligand = tmp_path / 'rec.pdb', tmp_path / 'lig.pdb'
    shutil.copy(FIXTURES / 'receptor_mini.pdb', receptor)
    shutil.copy(FIXTURES / 'ligand_mini.pdb', ligand)
    return receptor, ligand


@pytest.mark.parametrize('kind,distance,angle,required,status,reason', [
    ('hbond', 3.5, 80, True, 'No', 'distance_outside_cutoff;angle_outside_range'),
    ('hbond', 3, np.nan, True, 'No', 'missing_required_geometry'),
    ('hbond', 3, 100, True, 'No', 'angle_outside_range'),
    ('hbond', 3, 150, True, 'Yes', 'distance_and_angle_pass'),
    ('hbond', 3, np.nan, False, 'Yes', 'distance_only_no_probe_orientation'),
    ('hbond', 4, np.nan, False, 'No', 'distance_outside_cutoff'),
    ('aromatic', 5, 45, True, 'No', 'angle_outside_range'),
    ('aromatic', 5, np.nan, True, 'No', 'missing_required_geometry'),
    ('aromatic', 5, np.nan, False, 'Yes', 'distance_only_no_probe_orientation'),
    ('aromatic', 5, 20, True, 'Yes', 'distance_and_angle_pass'),
])
def test_all_failed_and_passed_criteria_are_explained(
        kind, distance, angle, required, status, reason):
    evaluate = evaluate_hbond if kind == 'hbond' else evaluate_aromatic
    assert evaluate(distance, angle, {}, angle_required=required) == (status, reason)
    assert set(reason.split(';')) <= REASON_CODES.keys()


@pytest.mark.parametrize('distance,angle', [(3.1999, 180), (3, 100.004)])
def test_probe_reason_uses_unrounded_geometry(cfg, distance, angle):
    site = pd.DataFrame([[1, 'ALA', 'O', 0, 0, 0], [1, 'ALA', 'C', 1, 0, 0]],
                        columns=['Pos', 'Residue', 'Atom', 'X', 'Y', 'Z'])
    points = site.iloc[[0]].copy()
    points.insert(0, 'Type', 'Aceptor')
    x, y = distance * math.cos(math.radians(angle)), distance * math.sin(math.radians(angle))
    result = probe_interactions([(x, y, 0, ['donor'])], site, points, cfg)
    contact = result[result.Type == 'donor'].iloc[0]
    assert contact.Interaction == 'Yes'
    assert contact.Reason == 'distance_and_angle_pass'
    assert contact.Dist == round(distance, 3)
    assert contact.Angle == round(angle, 2)


def test_pair_records_config_hashes_history_and_explanations(inputs, cfg, tmp_path):
    receptor, ligand = inputs
    cfg['Distances_Hidrogen_Bonds'] = 4.5
    analyze_pair(str(receptor), str(ligand), 'A', cfg)
    folder = tmp_path / 'rec_lig'
    metadata = json.loads((folder / 'run_metadata.json').read_text())
    assert metadata['status'] == 'completed'
    assert metadata['mode'] == 'ligand'
    assert metadata['interaction_csv_schema_version'] == 2
    assert metadata['effective_hbond_search_distance_A'] == 4.5
    assert metadata['parameters']['chain_receptor'] == 'A'
    assert metadata['inputs_before']['receptor']['sha256'] == hashlib.sha256(
        receptor.read_bytes()).hexdigest()
    assert metadata['inputs_before'] == metadata['inputs_after']
    assert len(metadata['software']['source_sha256']) == 64
    assert metadata['software']['python']
    assert metadata['software']['dependencies']['rdkit']
    assert metadata['duration_seconds'] >= 0
    assert metadata['finished_at_utc'] >= metadata['started_at_utc']
    assert _build_cfg(folder / 'config_used.yml') == cfg
    history = folder / 'run_history' / metadata['run_id']
    assert (history / 'run_metadata.json').read_bytes() == (
        folder / 'run_metadata.json').read_bytes()
    assert (history / 'config_used.yml').read_bytes() == (folder / 'config_used.yml').read_bytes()
    assert (folder / metadata['config_snapshot']['file']).is_file()
    for suffix in ('all', 'threshold', 'true'):
        name = f'Interaction_rec_lig_{suffix}.csv'
        frame = pd.read_csv(folder / name)
        assert frame.Reason.notna().all()
        assert set(';'.join(frame.Reason).split(';')) - {''} <= REASON_CODES.keys()
        assert metadata['outputs_written'][name]['sha256'] == hashlib.sha256(
            (folder / name).read_bytes()).hexdigest()


def test_failed_rerun_does_not_claim_stale_outputs(inputs, cfg, tmp_path, monkeypatch):
    receptor, ligand = inputs
    analyze_pair(str(receptor), str(ligand), 'A', cfg)
    folder = tmp_path / 'rec_lig'
    first = json.loads((folder / 'run_metadata.json').read_text())
    monkeypatch.setattr('interactions_search.pipeline.Chem.MolFromPDBFile', lambda *a, **k: None)
    with pytest.raises(ValueError, match='RDKit'):
        analyze_pair(str(receptor), str(ligand), 'A', cfg)
    second = json.loads((folder / 'run_metadata.json').read_text())
    assert second['status'] == 'failed'
    assert second['error']['type'] == 'ValueError'
    assert second['run_id'] != first['run_id']
    assert 'Interaction_rec_lig_all.csv' not in second['outputs_written']
    archived_path = folder / 'run_history' / first['run_id'] / 'run_metadata.json'
    archived = json.loads(archived_path.read_text())
    assert archived == first
    assert len(list((folder / 'run_history').iterdir())) == 2


def test_original_and_cleaned_ligand_are_distinguishable(inputs, cfg, tmp_path):
    receptor, ligand = inputs
    ligand.write_text(ligand.read_text() + 'REMARK dummy bias marker CM removed\n')
    before = hashlib.sha256(ligand.read_bytes()).hexdigest()
    analyze_pair(str(receptor), str(ligand), 'A', cfg)
    metadata = json.loads((tmp_path / 'rec_lig/run_metadata.json').read_text())
    assert metadata['inputs_before']['ligand']['sha256'] == before
    assert metadata['inputs_after']['ligand']['sha256'] != before
    assert metadata['outputs_written']['lig_old.pdb']['sha256'] == before


def test_reference_template_is_saved_and_reloadable(inputs, cfg, tmp_path):
    receptor, ligand = inputs
    reference = load_reference(smiles='OCCC')
    analyze_pair(str(receptor), str(ligand), 'A', cfg,
                 ligand_reference=reference, legacy_rings=True)
    folder = tmp_path / 'rec_lig'
    record = json.loads((folder / 'run_metadata.json').read_text())
    saved = folder / record['ligand_reference']['file']
    restored = load_reference(sdf=saved)
    assert [a.GetAtomicNum() for a in restored.GetAtoms()] == [
        a.GetAtomicNum() for a in reference.GetAtoms()]
    assert [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType()))
            for b in restored.GetBonds()] == [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType()))
        for b in reference.GetBonds()]
    assert record['parameters']['legacy_rings'] is True


def test_probe_records_effective_roles_and_clash_reasons(inputs, cfg, tmp_path, capsys):
    receptor, _ = inputs
    cfg['Probe_Clash_Distance'] = 3.0
    analyze_probe(str(receptor), 'A', [(0, 0, 0, ['acceptor'])], cfg, label='record')
    folder = tmp_path / 'rec_probe_record'
    record = json.loads((folder / 'run_metadata.json').read_text())
    assert record['parameters']['points'] == [[0, 0, 0, ['acceptor']]]
    assert record['runtime_config']['Probe_Clash_Distance'] == 3.0
    frame = pd.read_csv(folder / 'Probe_rec_all.csv')
    assert set(frame[frame.Type == 'clash'].Reason) == {'steric_clash'}
    assert set(frame[frame.Type == 'acceptor'].Reason) == {'distance_only_no_probe_orientation'}
    assert 'steric clash (< 3 Å)' in capsys.readouterr().out


def test_site_bias_records_method_and_radius(inputs, cfg, tmp_path):
    receptor, _ = inputs
    analyze_site_bias(str(receptor), 'A', (0, 0, 0), 8.0, cfg, method='atom')
    folder = tmp_path / 'rec_site_0.0_0.0_0.0'
    record = json.loads((folder / 'run_metadata.json').read_text())
    assert record['mode'] == 'site_bias'
    assert record['parameters']['method'] == 'atom'
    assert record['parameters']['radius'] == 8.0
    assert load_config(folder / 'config_used.yml').donors == cfg['Dadores_Prot']


def test_empty_results_still_record_success_and_schema(inputs, cfg, tmp_path):
    receptor, _ = inputs
    analyze_probe(str(receptor), 'A', [(100, 100, 100, ['donor'])], cfg, label='empty')
    folder = tmp_path / 'rec_probe_empty'
    record = json.loads((folder / 'run_metadata.json').read_text())
    frame = pd.read_csv(folder / 'Probe_rec_all.csv')
    assert record['status'] == 'completed'
    assert frame.empty and 'Reason' in frame.columns
