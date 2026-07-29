"""Detección de bolsillos hidrofóbicos: agrupa los contactos hidrofóbicos por
fragmento del ligando (conectividad de enlaces) y evalúa si un fragmento
está realmente envuelto por 3+ residuos (Coverage_R) o solo tocado
superficialmente, además de su volumen (convex hull) y densidad de contacto
local (estilo fpocket)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem

from interactions_search.contacts import _HPHO_LIG_SMARTS, _HYDROPHOBIC_ATOMS
from interactions_search.geometry import convex_hull_volume

__all__ = ["search_hydrophobic_pockets"]

_POCKET_SUMMARY_COLS = ['Pocket', 'Fragment_Atoms', 'N_Ligand_Atoms', 'Residues',
                        'N_Residues', 'Coverage_R', 'Volume_A3', 'Density_Score',
                        'Is_Pocket', 'X', 'Y', 'Z']
_POCKET_DETAIL_COLS  = ['Pocket', 'Pos R', 'Res', 'Atom', 'Lig_Atom', 'Lig_Serial', 'Dist']
_POCKET_COLORIDS     = [3, 9, 11, 4, 7, 10, 14, 17]  # orange, pink, purple2, yellow, green, cyan, ...


def _ligand_hydrophobic_fragments(mol, hpho_idx):
    """Componentes conexos por enlace (grafo de RDKit) dentro del set de átomos
    hidrofóbicos del ligando: un anillo o cadena contigua contactada = un
    fragmento. Se usa conectividad real, no cercanía espacial, para que 'mismo
    fragmento del ligando' tenga sentido químico."""
    idx_list = sorted(hpho_idx)
    parent = {i: i for i in idx_list}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in hpho_idx and b in hpho_idx:
            union(a, b)

    fragments = {}
    for i in idx_list:
        fragments.setdefault(find(i), []).append(i)
    return list(fragments.values())


def _hydrophobic_density_score(frag_contacts, pdb_coords, rec_rows, radius):
    """Densidad hidrofóbica local (estilo fpocket, 'mean local hydrophobic
    density'): cada contacto átomo(ligando)-átomo(receptor) del fragmento se
    representa por su punto medio; para cada punto medio se cuenta cuántos
    otros puntos medios del mismo fragmento caen dentro de 'radius' Å, y se
    devuelve el promedio de esos conteos.

    Es una métrica distinta de Coverage_R (dirección de los residuos) y de
    Volume_A3 (tamaño de la cavidad que los envuelve): mide qué tan apretados
    están los contactos átomo-átomo entre sí. Dos fragmentos con el mismo
    Coverage_R y Volume_A3 pueden tener un encaje hidrofóbico apretado (muchos
    contactos cercanos entre sí → score alto) o disperso dentro de la misma
    cavidad (score bajo). 0.0 si hay menos de 2 contactos (no hay par de
    puntos para medir distancia)."""
    if len(frag_contacts) < 2:
        return 0.0
    midpoints = []
    for lig_idx, lig_serial, lig_name, dist, rec_idx in frag_contacts:
        lig_xyz = np.array([pdb_coords[lig_idx][5], pdb_coords[lig_idx][6], pdb_coords[lig_idx][7]],
                           dtype=float)
        rec_xyz = rec_rows.loc[rec_idx][['X', 'Y', 'Z']].values.astype(float)
        midpoints.append((lig_xyz + rec_xyz) / 2)
    midpoints = np.array(midpoints)
    dist_matrix = np.linalg.norm(midpoints[:, None, :] - midpoints[None, :, :], axis=-1)
    np.fill_diagonal(dist_matrix, np.inf)
    return round(float((dist_matrix < radius).sum(axis=1).mean()), 3)


def search_hydrophobic_pockets(mol, pdb_coords, DF_Active_Site, Distancia_Hidrofobica,
                                min_residues=3, coverage_threshold=0.5, density_radius=5.0):
    """Detecta pockets hidrofóbicos reales: fragmentos contiguos del ligando
    (por conectividad) contactados por 3+ residuos distintos con buena
    cobertura espacial alrededor del fragmento — consistente con sitios de
    unión bien definidos en estructuras cristalográficas.

    Cobertura (Coverage_R): módulo del vector resultante normalizado de las
    direcciones residuo→fragmento (0 a 1). R bajo = los residuos rodean el
    fragmento desde varias direcciones (pocket real). R alto (cercano a 1) =
    todos los residuos del mismo lado (contacto superficial, no un pocket
    envolvente), aunque haya 3+ residuos.

    Retorna (df_summary, df_detail, pocket_hulls): df_summary tiene una fila por
    fragmento candidato (pase o no el filtro), con X/Y/Z = centroide de los
    átomos del ligando efectivamente en contacto (frag_center, usado también
    para Coverage_R), Volume_A3 = volumen (Å³) de la envolvente convexa
    (ConvexHull) de todos los átomos de los residuos que contactan el pocket
    (NaN si hay menos de 4 átomos o si son coplanares/degenerados), y
    Density_Score = densidad hidrofóbica local estilo fpocket, ver
    _hydrophobic_density_score (0.0 si el fragmento tiene < 2 contactos);
    df_detail solo tiene los contactos átomo-residuo de los fragmentos que sí
    califican como pocket (Is_Pocket == 'Yes'), para alimentar la visualización VMD;
    pocket_hulls es {pocket_n: (points, ConvexHull)} para los pockets con
    volumen calculable, reutilizado por plot_hull_volume() para no volver a
    correr Qhull."""
    pattern  = Chem.MolFromSmarts(_HPHO_LIG_SMARTS)
    hpho_idx = {i for match in mol.GetSubstructMatches(pattern) for i in match}
    if not hpho_idx:
        return pd.DataFrame(columns=_POCKET_SUMMARY_COLS), pd.DataFrame(columns=_POCKET_DETAIL_COLS), {}

    rec_mask = DF_Active_Site.apply(
        lambda r: r['Atom'] in _HYDROPHOBIC_ATOMS.get(r['Residue'], set()), axis=1)
    rec_rows = DF_Active_Site[rec_mask]
    if rec_rows.empty:
        return pd.DataFrame(columns=_POCKET_SUMMARY_COLS), pd.DataFrame(columns=_POCKET_DETAIL_COLS), {}
    rec_coords = np.array(rec_rows[['X', 'Y', 'Z']], dtype=float)

    # Contactos átomo(ligando)-átomo(receptor) crudos, sin colapsar por residuo
    raw_contacts = []  # (lig_atom_idx, lig_serial, lig_name, dist, rec_row_index)
    for i in hpho_idx:
        if i >= len(pdb_coords):
            continue
        lig_serial, lig_name = pdb_coords[i][0], pdb_coords[i][1]
        lig_xyz = np.array([pdb_coords[i][5], pdb_coords[i][6], pdb_coords[i][7]], dtype=float)
        dists = np.linalg.norm(rec_coords - lig_xyz, axis=1)
        for k in np.where(dists < Distancia_Hidrofobica)[0]:
            raw_contacts.append((i, lig_serial, lig_name, float(dists[k]), rec_rows.index[k]))

    if not raw_contacts:
        return pd.DataFrame(columns=_POCKET_SUMMARY_COLS), pd.DataFrame(columns=_POCKET_DETAIL_COLS), {}

    fragments = _ligand_hydrophobic_fragments(mol, hpho_idx)

    summary_rows, detail_rows = [], []
    pocket_hulls = {}
    for pocket_n, frag_atoms in enumerate(fragments, start=1):
        frag_set = set(frag_atoms)
        frag_contacts = [c for c in raw_contacts if c[0] in frag_set]
        if not frag_contacts:
            continue

        residues = {}  # Pos -> {'Residue': str, 'xyz': [[x,y,z], ...]}
        contacted_lig_atoms = set()
        for lig_idx, lig_serial, lig_name, dist, rec_idx in frag_contacts:
            rec_row = rec_rows.loc[rec_idx]
            pos = int(rec_row['Pos'])
            residues.setdefault(pos, {'Residue': rec_row['Residue'], 'xyz': []})
            residues[pos]['xyz'].append([rec_row['X'], rec_row['Y'], rec_row['Z']])
            contacted_lig_atoms.add(lig_idx)

        n_residues = len(residues)
        if n_residues < min_residues:
            continue

        frag_xyz = np.array([[pdb_coords[i][5], pdb_coords[i][6], pdb_coords[i][7]]
                             for i in contacted_lig_atoms], dtype=float)
        frag_center = frag_xyz.mean(axis=0)

        vectors = []
        for pos, info in residues.items():
            res_center = np.mean(info['xyz'], axis=0)
            v = res_center - frag_center
            norm = np.linalg.norm(v)
            if norm > 1e-6:
                vectors.append(v / norm)
        resultant = float(np.linalg.norm(np.sum(vectors, axis=0)) / len(vectors)) if vectors else 1.0
        is_pocket = (n_residues >= min_residues) and (resultant < coverage_threshold)

        residues_str    = ','.join(f"{info['Residue']}{pos}" for pos, info in sorted(residues.items()))
        frag_atoms_str  = ','.join(sorted({pdb_coords[i][1] for i in frag_atoms if i < len(pdb_coords)}))

        # Volumen: envolvente convexa de TODOS los átomos de los residuos que
        # contactan (no solo los apolares) — el mismo conjunto de átomos que
        # scripting_vmd_pockets selecciona para la representación Surf.
        pocket_points = DF_Active_Site[DF_Active_Site['Pos'].isin(residues.keys())][['X', 'Y', 'Z']] \
            .values.astype(float)
        volume, hull = convex_hull_volume(pocket_points)
        if hull is not None:
            pocket_hulls[pocket_n] = (pocket_points, hull)

        density = _hydrophobic_density_score(frag_contacts, pdb_coords, rec_rows, density_radius)

        summary_rows.append([pocket_n, frag_atoms_str, len(contacted_lig_atoms),
                             residues_str, n_residues, round(resultant, 3), volume, density,
                             'Yes' if is_pocket else 'No', *np.round(frag_center, 3)])

        if is_pocket:
            for lig_idx, lig_serial, lig_name, dist, rec_idx in frag_contacts:
                rec_row = rec_rows.loc[rec_idx]
                detail_rows.append([pocket_n, int(rec_row['Pos']), rec_row['Residue'], rec_row['Atom'],
                                    lig_name, lig_serial, round(dist, 3)])

    df_summary = pd.DataFrame(summary_rows, columns=_POCKET_SUMMARY_COLS) if summary_rows \
        else pd.DataFrame(columns=_POCKET_SUMMARY_COLS)
    df_detail  = pd.DataFrame(detail_rows, columns=_POCKET_DETAIL_COLS) if detail_rows \
        else pd.DataFrame(columns=_POCKET_DETAIL_COLS)
    return df_summary, df_detail, pocket_hulls
