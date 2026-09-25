"""Tests del modo sondeo (probe.py / analyze_probe) contra receptor_mini.pdb.

receptor_mini: ALA 10 con N en (-2.5, 0, 0) y CB en (-2, 2.3, 0); LEU 20 con
CB/CG/CD1/CD2 alrededor de x = 6-7.5. Sin H explícitos, sin aromáticos ni cargas.
"""
from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from interactions_search.cli import _build_cfg
from interactions_search.pipeline import analyze_probe
from interactions_search.probe import PROBE_TYPES, probe_interactions, read_probe_file
from interactions_search.receptor_site import Coordenadas_interes_receptor, active_site_residues

RECEPTOR = Path(__file__).parent / "fixtures" / "receptor_mini.pdb"


@pytest.fixture(scope="module")
def cfg():
    return _build_cfg(None)


@pytest.fixture(scope="module")
def site(cfg):
    structure = PDBParser(QUIET=True).get_structure("pdb", str(RECEPTOR))
    df_site = active_site_residues(structure, [0.0, 0.0, 0.0], "A", 30.0, "")
    rec_points = Coordenadas_interes_receptor(cfg["Aceptores_Prot"], cfg["Dadores_Prot"], df_site)
    return df_site, rec_points


def test_acceptor_probe_finds_backbone_n(cfg, site):
    df = probe_interactions([(0.0, 0.0, 0.0, ["acceptor"])], *site, cfg)
    hb = df[(df["Type"] == "acceptor") & (df["Interaction"] == "Yes")]
    assert list(zip(hb["Pos R"], hb["Atom"])) == [(10, "N")]
    assert hb["Dist"].iloc[0] == pytest.approx(2.5)


def test_clash_reported_and_not_validated(cfg, site):
    # (0,0,0) queda a 1.49 Å del CA de ALA10
    df = probe_interactions([(0.0, 0.0, 0.0, ["acceptor"])], *site, cfg)
    clash = df[df["Type"] == "clash"]
    assert "CA" in set(clash["Atom"])
    assert (clash["Interaction"] == "Clash").all()


def test_hydrophobic_probe_collapses_by_residue(cfg, site):
    df = probe_interactions([(6.5, 1.0, 3.0, ["hydrophobic"])], *site, cfg)
    hpho = df[df["Type"] == "hydrophobic"]
    assert list(hpho["Pos R"]) == [20]
    assert "," in hpho["Atom"].iloc[0]


def test_far_probe_is_empty(cfg, site):
    df = probe_interactions([(100.0, 100.0, 100.0, list(PROBE_TYPES))], *site, cfg)
    assert df.empty


def test_read_probe_file_bpf_and_csv(tmp_path):
    bpf = tmp_path / "pts.bpf"
    bpf.write_text("x\ty\tz\tVset\tr\ttype\n1.0\t2.0\t3.0\t-2.72\t1.2\tdon\n4\t5\t6\t-2.0\t2.0\taro\n")
    assert read_probe_file(bpf) == [(1.0, 2.0, 3.0, "donor"), (4.0, 5.0, 6.0, "aromatic")]

    csv = tmp_path / "pts.csv"
    csv.write_text("x,y,z,type\n1,2,3,cation\n7,8,9\n")
    assert read_probe_file(csv) == [(1.0, 2.0, 3.0, "cation"), (7.0, 8.0, 9.0, None)]


def test_analyze_probe_writes_outputs(cfg, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    analyze_probe(str(RECEPTOR), "A", [(0.0, 0.0, 0.0, list(PROBE_TYPES))], cfg, label="t")
    folder = tmp_path / "receptor_mini_probe_t"
    for name in ("Probe_receptor_mini_all.csv", "Probe_receptor_mini_true.csv",
                 "receptor_mini_probe_points.pdb", "receptor_mini.pdb"):
        assert (folder / name).exists(), name
