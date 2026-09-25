"""Pockets a partir de hotspots de dinámica (clusters de aceptor/donor/
hidrofóbico de una simulación con cosolvente): agrupa los hotspots cercanos en
sitios, arma el pocket de cada sitio (residuos a residue_cutoff Å de los
puntos del hotspot) y lo rellena con una grilla estilo AutoDock anotada con
propiedades.

Entrada: un directorio con subcarpetas acceptors/, donors/ y hydrophobics/
(alcanza con las que existan), cada una con:
  clusters.csv        ws_id, x, y, z, R90_A, DG, occ_prob, ... (un cluster por fila)
  cluster_points.pdb  posiciones crudas de la sonda, resid = ws_id (opcional)
  grid_*.dx           ΔG por vóxel de ese tipo, kcal/mol, <= 0 (opcional)

Por sitio:
  1. Hotspots: los centros se agrupan por single-linkage (link_distance). Cada
     hotspot aporta sus puntos de cluster_points.pdb dentro de R90 del centro
     (el 90% central; si no hay puntos, solo el centro).
  2. Residuos del pocket: los que tienen un átomo pesado a <= residue_cutoff Å
     de algún punto de un hotspot del sitio.
  3. Grilla (grid_spacing Å) en la caja que envuelve residuos + puntos. Un punto
     queda si: está dentro de la envolvente convexa de los átomos pesados de
     los residuos del pocket; no choca con la proteína (> grid_clash Å de todo
     átomo pesado); está enterrado (>= grid_buriedness de 30 rayos chocan con
     la proteína dentro de 10 Å, estilo LIGSITE); y está conectado (vecindad
     26) a algún punto de hotspot del sitio.
  4. Anotación de cada punto: ΔG de cada .dx en ese punto, Best_Type (tipo de
     menor ΔG si llega a grid_dg_threshold, si no 'none'), y el entorno
     químico del receptor: cuántos donores/aceptores/átomos hidrofóbicos/
     anillos/cationes/aniones del receptor quedan a distancia de interacción
     (mismos umbrales que el modo sondeo, solo distancia)."""
from __future__ import annotations

import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy import ndimage
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import Delaunay, QhullError, cKDTree

from interactions_search.contacts import (
    _CATION_REC_ATOMS,
    _HYDROPHOBIC_ATOMS,
    _SALT_NEG_ATOMS,
    _SALT_POS_ATOMS,
)
from interactions_search.interaction_rules import distance_cutoff
from interactions_search.provenance import record_analysis
from interactions_search.receptor_site import Coordenadas_interes_receptor

__all__ = ["read_dx", "write_dx", "load_hotspots", "group_sites", "analyze_hotspot_pockets"]

# carpeta -> (tipo, prefijo del id de hotspot, resname en los PDB de salida)
_HS_FOLDERS = {
    'acceptors':    ('acceptor',    'acc', 'ACC'),
    'donors':       ('donor',       'don', 'DON'),
    'hydrophobics': ('hydrophobic', 'hyd', 'HPH'),
}
_TYPE_RESNAME = {t: rn for t, _, rn in _HS_FOLDERS.values()}
_TYPE_RESNAME['none'] = 'NON'
_TYPE_COLORID = {'ACC': 1, 'DON': 0, 'HPH': 3, 'NON': 8}  # red, blue, orange, white

_WATER_NAMES = {'HOH', 'WAT', 'TIP', 'TIP3', 'SOL', 'DOD'}

# Enterramiento: 30 direcciones ~uniformes (espiral de Fibonacci), muestreadas
# cada 1.5 Å hasta 10 Å; un rayo "choca" si pasa a < 2 Å de un átomo pesado.
_RAY_STEPS = np.arange(1.5, 10.01, 1.5)
_RAY_HIT = 2.0


def _fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]


_RAY_DIRS = _fibonacci_sphere(30)


# ──────────────────────────────────────────────────────────────────────────────
# OpenDX
# ──────────────────────────────────────────────────────────────────────────────

def read_dx(path):
    """Lee una grilla OpenDX regular (ortogonal). Devuelve dict con origin (3,),
    delta (3,), values (nx, ny, nz) -- el índice z es el que varía más rápido."""
    text = Path(path).read_text()
    header, data = text.split('data follows', 1)
    counts = [int(c) for c in re.search(r'gridpositions counts\s+(\d+)\s+(\d+)\s+(\d+)',
                                        header).groups()]
    origin = np.array(re.search(r'origin\s+(\S+)\s+(\S+)\s+(\S+)', header).groups(), dtype=float)
    deltas = [np.array(m, dtype=float)
              for m in re.findall(r'delta\s+(\S+)\s+(\S+)\s+(\S+)', header)]
    delta = np.array([deltas[0][0], deltas[1][1], deltas[2][2]])
    values = np.array(re.split(r'attribute|object', data, maxsplit=1)[0].split(), dtype=float)
    return {'origin': origin, 'delta': delta, 'values': values.reshape(counts)}


def write_dx(path, origin, delta, values):
    nx, ny, nz = values.shape
    flat = values.ravel()
    with open(path, 'w') as f:
        f.write(f'object 1 class gridpositions counts {nx} {ny} {nz}\n')
        f.write(f'origin {origin[0]:.6f} {origin[1]:.6f} {origin[2]:.6f}\n')
        f.write(f'delta {delta:.6f} 0 0\ndelta 0 {delta:.6f} 0\ndelta 0 0 {delta:.6f}\n')
        f.write(f'object 2 class gridconnections counts {nx} {ny} {nz}\n')
        f.write(f'object 3 class array type double rank 0 items {flat.size} data follows\n')
        for i in range(0, flat.size, 3):
            f.write(' '.join(f'{v:.6e}' for v in flat[i:i + 3]) + '\n')
        f.write('attribute "dep" string "positions"\n')
        f.write('object "regular positions regular connections" class field\n')
        f.write('component "positions" value 1\ncomponent "connections" value 2\n')
        f.write('component "data" value 3\n')


def _dx_lookup(dx, points):
    """Valor del vóxel más cercano para cada punto (NaN fuera de la grilla)."""
    idx = np.rint((points - dx['origin']) / dx['delta']).astype(int)
    shape = np.array(dx['values'].shape)
    inside = np.all((idx >= 0) & (idx < shape), axis=1)
    out = np.full(len(points), np.nan)
    out[inside] = dx['values'][tuple(idx[inside].T)]
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Hotspots
# ──────────────────────────────────────────────────────────────────────────────

def _read_cluster_points(path):
    pts = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(('ATOM', 'HETATM')):
                pts.setdefault(int(line[22:26]), []).append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return {k: np.array(v) for k, v in pts.items()}


def load_hotspots(results_dir):
    """Devuelve (hs, points, dx_maps):
      hs      DataFrame: HS_ID, Type, ws_id, X, Y, Z, R90, DG, occ_prob
      points  {HS_ID: array (n, 3)} -- puntos del cluster dentro de R90 del centro
      dx_maps {Type: dict de read_dx}"""
    results_dir = Path(results_dir)
    rows, points, dx_maps = [], {}, {}
    for folder, (tipo, prefix, _) in _HS_FOLDERS.items():
        d = results_dir / folder
        if not (d / 'clusters.csv').exists():
            continue
        df = pd.read_csv(d / 'clusters.csv')
        raw = _read_cluster_points(d / 'cluster_points.pdb') \
            if (d / 'cluster_points.pdb').exists() else {}
        for _, r in df.iterrows():
            hs_id = f'{prefix}{int(r["ws_id"])}'
            center = np.array([r['x'], r['y'], r['z']], dtype=float)
            pts = raw.get(int(r['ws_id']))
            if pts is not None:
                pts = pts[np.linalg.norm(pts - center, axis=1) <= r['R90_A']]
            points[hs_id] = pts if pts is not None and len(pts) else center[None, :]
            rows.append([hs_id, tipo, int(r['ws_id']), *center, float(r['R90_A']),
                         float(r['DG']), float(r.get('occ_prob', np.nan))])
        dx_files = sorted(d.glob('*.dx'))
        if dx_files:
            dx_maps[tipo] = read_dx(dx_files[0])
    if not rows:
        raise FileNotFoundError(f"No se encontró ningún clusters.csv en {results_dir}/"
                                f"{{{','.join(_HS_FOLDERS)}}}")
    hs = pd.DataFrame(rows, columns=['HS_ID', 'Type', 'ws_id', 'X', 'Y', 'Z', 'R90', 'DG',
                                     'occ_prob'])
    return hs, points, dx_maps


def group_sites(hs, link_distance):
    """Agrega la columna Site: single-linkage de los centros con corte
    link_distance. Sitios numerados por ΔG sumado (el más favorable = 1)."""
    hs = hs.copy()
    if len(hs) == 1:
        hs['Site'] = 1
        return hs
    labels = fcluster(linkage(hs[['X', 'Y', 'Z']].values, 'single'), link_distance, 'distance')
    order = pd.Series(hs['DG'].values).groupby(labels).sum().sort_values().index
    rank = {lab: n for n, lab in enumerate(order, start=1)}
    hs['Site'] = [rank[lab] for lab in labels]
    return hs.sort_values(['Site', 'DG']).reset_index(drop=True)


# ──────────────────────────────────────────────────────────────────────────────
# Receptor
# ──────────────────────────────────────────────────────────────────────────────

def _receptor_atoms(receptor_pdb, chain, exclude_resnames):
    """Átomos de la cadena como DataFrame con las columnas que espera
    Coordenadas_interes_receptor() (Serial, Pos, Residue, Atom, X, Y, Z) +
    Element. Excluye agua y los resnames de exclude_resnames (ej. ligandos
    guardados como ATOM, que si no ocuparían el pocket)."""
    structure = PDBParser(QUIET=True).get_structure('pdb', receptor_pdb)
    exclude = {n.upper() for n in exclude_resnames} | _WATER_NAMES
    rows = []
    for residue in structure[0][chain]:
        resname = residue.get_resname()
        if resname in exclude:
            continue
        for atom in residue:
            x, y, z = (round(float(c), 3) for c in atom.get_coord())
            rows.append([atom.get_serial_number(), residue.get_id()[1], resname,
                         atom.get_name(), x, y, z, (atom.element or '').upper()])
    return pd.DataFrame(rows, columns=['Serial', 'Pos', 'Residue', 'Atom', 'X', 'Y', 'Z',
                                       'Element'])


def _pocket_residues(site_hs, points, heavy, heavy_tree, cutoff):
    """Residuos con algún átomo pesado a <= cutoff de un punto de los hotspots."""
    per_res = {}
    for _, h in site_hs.iterrows():
        pts = points[h['HS_ID']]
        idx = sorted({i for lst in heavy_tree.query_ball_point(pts, cutoff) for i in lst})
        if not idx:
            continue
        dists, _ = cKDTree(pts).query(heavy.iloc[idx][['X', 'Y', 'Z']].values)
        for i, d in zip(idx, dists):
            atom = heavy.iloc[i]
            key = (int(atom['Pos']), atom['Residue'])
            info = per_res.setdefault(key, {'Min_Dist': np.inf, 'atoms': set(), 'hs': set(),
                                            'types': set()})
            info['Min_Dist'] = min(info['Min_Dist'], d)
            info['atoms'].add(atom['Atom'])
            info['hs'].add(h['HS_ID'])
            info['types'].add(h['Type'])
    rows = [[pos, res, round(float(v['Min_Dist']), 3), len(v['atoms']),
             ','.join(sorted(v['hs'])), ','.join(sorted(v['types']))]
            for (pos, res), v in sorted(per_res.items())]
    return pd.DataFrame(rows, columns=['Pos', 'Residue', 'Min_Dist', 'N_Atoms', 'Hotspots',
                                       'Types'])


# ──────────────────────────────────────────────────────────────────────────────
# Grilla del pocket
# ──────────────────────────────────────────────────────────────────────────────

def _buriedness(points, heavy_tree):
    samples = (points[:, None, None, :]
               + _RAY_DIRS[None, :, None, :] * _RAY_STEPS[None, None, :, None])
    d, _ = heavy_tree.query(samples.reshape(-1, 3), distance_upper_bound=_RAY_HIT, workers=-1)
    hits = np.isfinite(d).reshape(len(points), len(_RAY_DIRS), len(_RAY_STEPS)).any(axis=2)
    return hits.mean(axis=1)


def _pocket_grid(site_pts, res_xyz, heavy_tree, hp):
    """Devuelve (origin, shape, mask 3D bool, buriedness 3D) de la grilla del pocket."""
    spacing = hp.grid_spacing
    lo = np.minimum(res_xyz.min(axis=0), site_pts.min(axis=0))
    hi = np.maximum(res_xyz.max(axis=0), site_pts.max(axis=0))
    shape = tuple(np.ceil((hi - lo) / spacing).astype(int) + 1)
    idx = np.indices(shape).reshape(3, -1).T
    xyz = lo + idx * spacing

    mask = np.zeros(len(xyz), dtype=bool)
    try:
        mask = Delaunay(res_xyz).find_simplex(xyz) >= 0
    except QhullError:
        pass  # < 4 átomos o coplanares: no hay volumen que encerrar
    cand = np.where(mask)[0]
    if len(cand):
        d, _ = heavy_tree.query(xyz[cand], distance_upper_bound=hp.grid_clash, workers=-1)
        cand = cand[~np.isfinite(d)]
    bur = np.zeros(len(xyz))
    if len(cand):
        bur[cand] = _buriedness(xyz[cand], heavy_tree)
        cand = cand[bur[cand] >= hp.grid_buriedness]
    mask = np.zeros(len(xyz), dtype=bool)
    mask[cand] = True
    mask = mask.reshape(shape)

    # Conectividad: solo componentes (vecindad 26) que tocan algún punto de hotspot
    # (a <= 2 Å: el centro del hotspot puede caer en un vóxel descartado por choque).
    labels, n = ndimage.label(mask, structure=np.ones((3, 3, 3)))
    if n:
        kept = np.argwhere(mask)
        d, _ = cKDTree(site_pts).query(lo + kept * spacing, distance_upper_bound=2.0)
        seeds = set(labels[tuple(kept[np.isfinite(d)].T)])
        mask = np.isin(labels, list(seeds))
    return lo, shape, mask, bur.reshape(shape)


def _annotate_grid(xyz, dx_maps, receptor_points, rec_atoms, cfg, hp):
    out = pd.DataFrame({'X': xyz[:, 0].round(3), 'Y': xyz[:, 1].round(3),
                        'Z': xyz[:, 2].round(3)})
    dg_cols = []
    for tipo in ('acceptor', 'donor', 'hydrophobic'):
        if tipo in dx_maps:
            out[f'DG_{tipo}'] = _dx_lookup(dx_maps[tipo], xyz).round(3)
            dg_cols.append(f'DG_{tipo}')
    if dg_cols:
        dg = out[dg_cols].to_numpy()
        valid = ~np.all(np.isnan(dg), axis=1)
        best = np.full(len(out), np.nan)
        best[valid] = np.nanmin(dg[valid], axis=1)
        best_type = np.array(['none'] * len(out), dtype=object)
        ok = valid & (best <= hp.grid_dg_threshold)
        best_type[ok] = [dg_cols[i][3:] for i in np.nanargmin(dg[ok], axis=1)]
        out['Best_Type'], out['Best_DG'] = best_type, best.round(3)
    else:
        out['Best_Type'], out['Best_DG'] = 'none', np.nan

    def count_near(sub_xyz, cutoff):
        if len(sub_xyz) == 0:
            return np.zeros(len(xyz), dtype=int)
        # query_ball_point includes the radius; interaction cutoffs are strict.
        return cKDTree(sub_xyz).query_ball_point(
            xyz, np.nextafter(float(cutoff), -np.inf), return_length=True)

    def atoms_in(table):
        mask = [a in table.get(r, set()) for r, a in zip(rec_atoms['Residue'], rec_atoms['Atom'])]
        return rec_atoms[mask][['X', 'Y', 'Z']].to_numpy(dtype=float)

    rp = receptor_points
    rp_xyz = lambda t: rp[rp['Type'] == t][['X', 'Y', 'Z']].to_numpy(dtype=float)  # noqa: E731
    out['N_Rec_Donors']      = count_near(rp_xyz('Dador'), distance_cutoff('donor', cfg))
    out['N_Rec_Acceptors']   = count_near(rp_xyz('Aceptor'), distance_cutoff('acceptor', cfg))
    out['N_Rec_Hydrophobic'] = count_near(atoms_in(_HYDROPHOBIC_ATOMS),
                                          distance_cutoff('hydrophobic', cfg))
    out['N_Rec_Aromatic']    = count_near(rp_xyz('aromatic'), distance_cutoff('aromatic', cfg))
    out['N_Rec_Cations']     = count_near(atoms_in({**_SALT_POS_ATOMS, **_CATION_REC_ATOMS}),
                                          distance_cutoff('salt_bridge', cfg))
    out['N_Rec_Anions']      = count_near(atoms_in(_SALT_NEG_ATOMS),
                                        distance_cutoff('salt_bridge', cfg))
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Salidas
# ──────────────────────────────────────────────────────────────────────────────

def _write_points_pdb(path, xyz, resnames, occupancy, beta):
    with open(path, 'w') as f:
        for i, (p, rn, occ, b) in enumerate(zip(xyz, resnames, occupancy, beta), start=1):
            f.write(f"HETATM{i % 100000:>5} {'C':<4} {rn:>3} X{i % 10000:>4}    "
                    f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}{occ:6.2f}{b:6.2f}          {'C':>2}\n")
        f.write('END\n')


def _write_site_tcl(path, receptor_file, chain, residues):
    resids = ' '.join(str(p) for p in residues['Pos'])
    with open(path, 'w') as f:
        f.write('display projection orthographic\n')
        f.write(f'set molReceptor [mol new "{receptor_file}"]\n')
        f.write('mol modselect 0 $molReceptor all\n')
        f.write('mol modstyle 0 $molReceptor Lines 3\n')
        f.write('mol modcolor 0 $molReceptor ColorID 6\n')
        if resids:
            f.write('mol addrep $molReceptor\n')
            f.write(f'mol modselect 1 $molReceptor "chain {chain} and resid {resids}"\n')
            f.write('mol modstyle 1 $molReceptor Licorice 0.2\n')
        f.write('\n# Grilla del pocket: un punto por vóxel, color por Best_Type\n')
        f.write('set molGrid [mol new "grid.pdb"]\n')
        for rep, (rn, cid) in enumerate(_TYPE_COLORID.items()):
            if rep:
                f.write('mol addrep $molGrid\n')
            f.write(f'mol modselect {rep} $molGrid "resname {rn}"\n')
            f.write(f'mol modstyle {rep} $molGrid Points 3\n')
            f.write(f'mol modcolor {rep} $molGrid ColorID {cid}\n')
        f.write('\n# Hotspots (centros), radio = tipo, color por tipo\n')
        f.write('set molHS [mol new "hotspots.pdb"]\n')
        for rep, rn in enumerate(('ACC', 'DON', 'HPH')):
            if rep:
                f.write('mol addrep $molHS\n')
            f.write(f'mol modselect {rep} $molHS "resname {rn}"\n')
            f.write(f'mol modstyle {rep} $molHS VDW 0.6 12\n')
            f.write(f'mol modcolor {rep} $molHS ColorID {_TYPE_COLORID[rn]}\n')
        f.write('\n# Forma del pocket: isosuperficie (wireframe) de la máscara\n')
        f.write('set molMask [mol new "pocket_mask.dx" type dx]\n')
        f.write('mol modstyle 0 $molMask Isosurface 0.5 0 0 1 1 1\n')
        f.write('mol modcolor 0 $molMask ColorID 2\n')
        f.write('display resetview\n')


@record_analysis('hotspots')
def analyze_hotspot_pockets(receptor_pdb, chain, hotspot_dir, cfg, hp, exclude_resnames=(), *,
                            run_context=None):
    """Arma un pocket por sitio de hotspots (ver docstring del módulo). cfg es el
    dict de _build_cfg() (distancias de interacción, tablas del YAML), hp la
    sección hotspot_pocket del YAML. Escribe en <receptor>_hotspot_pockets/ y
    devuelve el DataFrame de resumen por sitio."""
    receptor = Path(receptor_pdb).stem
    folder = Path(f'{receptor}_hotspot_pockets')
    folder.mkdir(exist_ok=True)

    hs, points, dx_maps = load_hotspots(hotspot_dir)
    hs = group_sites(hs, hp.link_distance)

    rec_atoms = _receptor_atoms(receptor_pdb, chain, exclude_resnames)
    heavy = rec_atoms[rec_atoms['Element'] != 'H'].reset_index(drop=True)
    heavy_tree = cKDTree(heavy[['X', 'Y', 'Z']].to_numpy(dtype=float))

    by_type = ', '.join(f'{t}: {n}' for t, n in hs['Type'].value_counts().items())
    print(f"\n  Hotspots: {len(hs)} ({by_type}) -> {hs['Site'].nunique()} site(s) "
          f"(link {hp.link_distance} Å)")
    if exclude_resnames:
        print(f"  Excluded from receptor: {', '.join(exclude_resnames)}")

    summary = []
    for site, site_hs in hs.groupby('Site'):
        n_by_type = site_hs['Type'].value_counts()
        center = site_hs[['X', 'Y', 'Z']].mean().to_numpy()
        row = {'Site': site, 'N_Hotspots': len(site_hs),
               'N_Acceptor': int(n_by_type.get('acceptor', 0)),
               'N_Donor': int(n_by_type.get('donor', 0)),
               'N_Hydrophobic': int(n_by_type.get('hydrophobic', 0)),
               'DG_Sum': round(site_hs['DG'].sum(), 3), 'DG_Best': round(site_hs['DG'].min(), 3),
               'Center_X': round(center[0], 3), 'Center_Y': round(center[1], 3),
               'Center_Z': round(center[2], 3), 'Hotspots': ','.join(site_hs['HS_ID'])}
        if len(site_hs) < hp.min_hotspots:
            summary.append({**row, 'Built': 'No'})
            continue

        site_dir = folder / f'site_{site}'
        site_dir.mkdir(exist_ok=True)
        site_pts = np.vstack([points[h] for h in site_hs['HS_ID']])

        residues = _pocket_residues(site_hs, points, heavy, heavy_tree, hp.residue_cutoff)
        res_heavy = heavy[heavy['Pos'].isin(residues['Pos'])]
        res_xyz = res_heavy[['X', 'Y', 'Z']].to_numpy(dtype=float)

        if len(res_xyz) >= 4:
            origin, shape, mask, bur = _pocket_grid(site_pts, res_xyz, heavy_tree, hp)
        else:
            origin, shape = site_pts.min(axis=0), (1, 1, 1)
            mask, bur = np.zeros(shape, dtype=bool), np.zeros(shape)
        kept = np.argwhere(mask)
        xyz = origin + kept * hp.grid_spacing

        # Puntos de interés del receptor (donores/aceptores/anillos) solo de los
        # residuos cercanos al sitio: todo lo que pueda quedar a distancia de
        # interacción de la grilla.
        near = rec_atoms[rec_atoms['Pos'].isin(
            heavy.iloc[sorted({i for lst in heavy_tree.query_ball_point(site_pts, 15.0)
                               for i in lst})]['Pos'])].reset_index(drop=True)
        receptor_points = Coordenadas_interes_receptor(cfg['Aceptores_Prot'], cfg['Dadores_Prot'],
                                                       near)
        grid = _annotate_grid(xyz, dx_maps, receptor_points, near, cfg, hp)
        grid.insert(3, 'Buriedness', bur[tuple(kept.T)].round(3) if len(kept) else [])

        site_hs.to_csv(site_dir / 'hotspots.csv', index=False)
        residues.to_csv(site_dir / 'residues.csv', index=False)
        grid.to_csv(site_dir / 'grid.csv', index=False)
        _write_points_pdb(site_dir / 'grid.pdb', xyz,
                          [_TYPE_RESNAME[t] for t in grid['Best_Type']],
                          grid['Buriedness'], grid['Best_DG'].fillna(0.0))
        _write_points_pdb(site_dir / 'hotspots.pdb', site_hs[['X', 'Y', 'Z']].to_numpy(),
                          [_TYPE_RESNAME[t] for t in site_hs['Type']],
                          site_hs['R90'], site_hs['DG'])
        write_dx(site_dir / 'pocket_mask.dx', origin, hp.grid_spacing, mask.astype(float))
        shutil.copy(receptor_pdb, site_dir / Path(receptor_pdb).name)
        _write_site_tcl(site_dir / f'vmd_site_{site}.tcl', Path(receptor_pdb).name, chain, residues)

        type_frac = grid['Best_Type'].value_counts(normalize=True) if len(grid) else pd.Series()
        row.update({'Built': 'Yes', 'N_Residues': len(residues),
                    'Residues': ','.join(f"{r}{p}" for p, r in zip(residues['Pos'],
                                                                   residues['Residue'])),
                    'Grid_Points': len(grid),
                    'Volume_A3': round(len(grid) * hp.grid_spacing ** 3, 1),
                    **{f'Frac_{t}': round(float(type_frac.get(t, 0.0)), 3)
                       for t in ('acceptor', 'donor', 'hydrophobic', 'none')}})
        summary.append(row)

        print(f"  {'─'*74}")
        print(f"  Site {site}: {len(site_hs)} hotspots (acc {row['N_Acceptor']}, "
              f"don {row['N_Donor']}, hyd {row['N_Hydrophobic']}), ΣΔG {row['DG_Sum']:.2f}, "
              f"best {row['DG_Best']:.2f} kcal/mol")
        print(f"    Residues ({len(residues)}): {row['Residues']}")
        print(f"    Pocket grid: {len(grid)} points, {row['Volume_A3']:.1f} Å³  "
              f"(acc {row['Frac_acceptor']:.0%}, don {row['Frac_donor']:.0%}, "
              f"hyd {row['Frac_hydrophobic']:.0%}, none {row['Frac_none']:.0%})")

    df_summary = pd.DataFrame(summary)
    df_summary.to_csv(folder / 'sites_summary.csv', index=False)
    skipped = df_summary[df_summary['Built'] == 'No']
    if not skipped.empty:
        print(f"  {'─'*74}")
        print(f"  {len(skipped)} site(s) with < {hp.min_hotspots} hotspots skipped: "
              f"{', '.join(skipped['Hotspots'])}")
    print(f"  -> {folder}/sites_summary.csv\n")
    return df_summary
