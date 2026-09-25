"""Tests del modo --hotspots (hotspot_pocket.py) con un directorio de hotspots
sintético sobre receptor_mini.pdb (ALA 10 alrededor de x = -2, LEU 20 alrededor
de x = 4-7.5)."""
import json
from pathlib import Path

import numpy as np
import pytest

from interactions_search.cli import _build_cfg
from interactions_search.config import HotspotPocket
from interactions_search.hotspot_pocket import (
    analyze_hotspot_pockets,
    group_sites,
    load_hotspots,
    read_dx,
    write_dx,
)

RECEPTOR = Path(__file__).parent / "fixtures" / "receptor_mini.pdb"
_HEADER = ("ws_id,count,x,y,z,center_x_A,center_y_A,center_z_A,R90,R90_A,hits,occ_prob,"
           "p_bulk,WFP,DG\n")


def _cluster_row(ws_id, xyz, dg, r90=1.0):
    x, y, z = xyz
    return f"{ws_id},10,{x},{y},{z},{x},{y},{z},{r90},{r90},5,0.1,0.001,50,{dg}\n"


def _points_pdb(points_by_id):
    lines = []
    for ws_id, pts in points_by_id.items():
        for x, y, z in pts:
            lines.append(f"HETATM{len(lines) + 1:>5}  O   WAT A{ws_id:>4}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  1.00           O\n")
    return "".join(lines)


@pytest.fixture
def hotspot_dir(tmp_path):
    root = tmp_path / "results_global"
    # acc1 junto al N de ALA10, don1 cerca de acc1 (mismo sitio), hyd1 junto a LEU20 y
    # hyd2 lejos (sitio propio).
    acc = root / "acceptors"
    acc.mkdir(parents=True)
    (acc / "clusters.csv").write_text(_HEADER + _cluster_row(1, (-1.0, -1.5, 0.0), -3.0))
    # el último punto queda fuera de R90
    (acc / "cluster_points.pdb").write_text(_points_pdb(
        {1: [(-1.0, -1.5, 0.0), (-1.2, -1.4, 0.2), (5.0, 5.0, 5.0)]}))
    don = root / "donors"
    don.mkdir()
    (don / "clusters.csv").write_text(_HEADER + _cluster_row(1, (1.5, -1.0, 0.0), -2.5))
    hyd = root / "hydrophobics"
    hyd.mkdir()
    (hyd / "clusters.csv").write_text(_HEADER + _cluster_row(1, (6.5, 1.0, 2.5), -2.0)
                                      + _cluster_row(2, (40.0, 40.0, 40.0), -2.2))
    return root


def test_load_hotspots_filters_points_by_r90(hotspot_dir):
    hs, points, dx_maps = load_hotspots(hotspot_dir)
    assert sorted(hs["HS_ID"]) == ["acc1", "don1", "hyd1", "hyd2"]
    assert len(points["acc1"]) == 2  # el punto a (5, 5, 5) queda fuera de R90
    assert points["don1"].shape == (1, 3)  # sin cluster_points.pdb: solo el centro
    assert dx_maps == {}


def test_group_sites_single_linkage(hotspot_dir):
    hs, _, _ = load_hotspots(hotspot_dir)
    hs = group_sites(hs, link_distance=8.0)
    site_of = dict(zip(hs["HS_ID"], hs["Site"]))
    assert site_of["acc1"] == site_of["don1"] == site_of["hyd1"]
    assert site_of["hyd2"] != site_of["acc1"]
    assert site_of["acc1"] == 1  # el de ΔG sumado más favorable


def test_dx_roundtrip(tmp_path):
    values = np.arange(24, dtype=float).reshape(2, 3, 4) * -0.1
    write_dx(tmp_path / "g.dx", np.array([1.0, 2.0, 3.0]), 0.5, values)
    dx = read_dx(tmp_path / "g.dx")
    assert np.allclose(dx["origin"], [1, 2, 3])
    assert np.allclose(dx["delta"], [0.5, 0.5, 0.5])
    assert np.allclose(dx["values"], values)


def test_analyze_hotspot_pockets(hotspot_dir, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    summary = analyze_hotspot_pockets(str(RECEPTOR), "A", hotspot_dir, _build_cfg(None),
                                      HotspotPocket(min_hotspots=2))
    built = summary[summary["Built"] == "Yes"]
    assert list(built["Site"]) == [1]
    assert set(built["Residues"].iloc[0].split(",")) == {"ALA10", "LEU20"}
    site_dir = tmp_path / "receptor_mini_hotspot_pockets" / "site_1"
    for name in ("hotspots.csv", "residues.csv", "grid.csv", "grid.pdb", "hotspots.pdb",
                 "pocket_mask.dx", "vmd_site_1.tcl", "receptor_mini.pdb"):
        assert (site_dir / name).exists(), name
    record = json.loads((site_dir.parent / 'run_metadata.json').read_text())
    assert record['mode'] == 'hotspots' and record['status'] == 'completed'
    assert any(name.startswith('hotspots/') for name in record['inputs_before'])
    assert 'site_1/grid.csv' in record['outputs_written']
    from interactions_search.config import load_config
    assert load_config(site_dir.parent / 'config_used.yml').hotspot_pocket.min_hotspots == 2
