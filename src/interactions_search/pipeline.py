"""Orquesta el pipeline completo de análisis para un par receptor-ligando:
hot-points del ligando, sitio activo del receptor, búsqueda de contactos,
validación por ángulo, y escritura de todas las salidas (CSV, PNG, .tcl)."""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from rdkit import Chem

from interactions_search.bias import (
    export_bpf,
    export_bpf_pdb,
    export_bpf_pdb_points,
    export_bpf_pdb_receptor,
    export_bpf_points,
    export_bpf_receptor,
)
from interactions_search.chi_angles import compute_active_site_chi_angles, compute_pocket_chi_angles
from interactions_search.config import load_config
from interactions_search.contacts import (
    Busqueda_Antecesor_Lig,
    Interaccion_Aromatica,
    residuos_contacto,
    search_hydrophobic,
    search_pi_cation,
    search_salt_bridges,
)
from interactions_search.geometry import angle_three_points, convex_hull_volume
from interactions_search.ideal_sites import ideal_site_points
from interactions_search.io_pdb import extract_coords_from_pdb, remove_bias
from interactions_search.interaction_rules import (
    distance_cutoff,
    evaluate_aromatic,
    evaluate_hbond,
    hbond_search_cutoff,
    neighbors,
    receptor_acceptor_angle,
    within_distance,
)
from interactions_search.ligand_chemistry import apply_reference
from interactions_search.ligand_hotpoints import (
    generate_df_ligand,
    search_hot_points,
    search_rings,
    visualize_rings,
)
from interactions_search.plotting import (
    plot_chi_profile,
    plot_hull_surface,
    plot_hull_volume,
    plot_ramachandran,
)
from interactions_search.pockets import search_hydrophobic_pockets
from interactions_search.probe import PROBE_TYPES, probe_interactions, write_probe_pdb
from interactions_search.provenance import record_analysis
from interactions_search.ramachandran import compute_active_site_phi_psi
from interactions_search.receptor_site import (
    Coordenadas_interes_receptor,
    active_site_residues,
)
from interactions_search.vmd import (
    scripting_vmd,
    scripting_vmd_combined,
    scripting_vmd_hydrophobic,
    scripting_vmd_pockets,
    scripting_vmd_probe,
)

__all__ = ["carga_variables", "add_interaction_coords", "print_summary", "analyze_pair",
           "analyze_site_bias", "analyze_probe"]


def carga_variables(config_path=None):
    cfg = load_config(config_path)
    return (
        cfg.options.ligand_plot,
        cfg.options.vmd_output,
        cfg.options.cumulative_output,
        cfg.options.interaction_coord,
        cfg.options.volume_plot,
        cfg.options.bias,
        cfg.options.bias_validated_only,
        cfg.distancias.Distances_Hidrogen_Bonds,
        cfg.distancias.Distances_Aromatic,
        cfg.distancias.Distances_Hidrofobica,
        cfg.distancias.centroid_distance,
        cfg.angulos.Angle_Hidrogen_Bonds_Min,
        cfg.angulos.Angle_Hidrogen_Bonds_Max,
        cfg.aromaticidad.Ring_Planarity_RMSD_Max,
        cfg.pockets.min_residues,
        cfg.pockets.coverage_threshold,
        cfg.pockets.density_radius,
        cfg.acceptors,
        cfg.donors,
        cfg.acceptors_antecedent,
        cfg.special,
    )


# ──────────────────────────────────────────────────────────────────────────────
# Coordenadas de las interacciones
# ──────────────────────────────────────────────────────────────────────────────

_COORD_SOURCES = {'receptor', 'ligand', 'center'}


def add_interaction_coords(DF, receptor_points, DF_Lig, DF_Lig_All, DF_Active_Site, coord_source):
    """Agrega columnas X, Y, Z a DF con la coordenada 3D de cada interacción,
    según coord_source (options.interaction_coord en el YAML):
      'receptor' -> átomo/centroide del receptor
      'ligand'   -> átomo/centroide del ligando
      'center'   -> punto medio entre ambos
    'aromatic'/'pi_cation' usan el centroide del anillo (columna 'Lig' = 'Caso'
    del anillo); el resto ubica el átomo puntual vía LigID (serial), no por
    nombre, porque el nombre puede repetirse en el ligando. 'hydrophobic' puede
    traer varios átomos del receptor colapsados en 'Atom' (ej. 'CD1,CD2,CG'); se
    promedian. Filas cuyo átomo/anillo no se puede resolver quedan con X/Y/Z NaN."""
    if coord_source not in _COORD_SOURCES:
        raise ValueError(f"coord_source inválido: {coord_source!r} (usar {_COORD_SOURCES})")
    need_rec = coord_source in ('receptor', 'center')
    need_lig = coord_source in ('ligand', 'center')

    xs, ys, zs = [], [], []
    for _, r in DF.iterrows():
        tipo, pos_r = r['Type'], r['Pos R']
        rec_xyz = lig_xyz = None

        if need_rec:
            if tipo == 'aromatic':
                sub = receptor_points[(receptor_points['Pos'] == pos_r) &
                                       (receptor_points['Atom'] == 'center')]
            else:
                atoms = str(r['Atom']).split(',')
                sub = DF_Active_Site[(DF_Active_Site['Pos'] == pos_r) &
                                      (DF_Active_Site['Atom'].isin(atoms))]
            if not sub.empty:
                rec_xyz = sub[['X', 'Y', 'Z']].astype(float).mean(axis=0).values

        if need_lig:
            if tipo in ('aromatic', 'pi_cation'):
                ring = DF_Lig[DF_Lig['Caso'] == r['Lig']]
                if not ring.empty:
                    lig_xyz = ring[['Coord X', 'Coord Y', 'Coord Z']].astype(float).mean(axis=0).values
            else:
                lig_row = DF_Lig_All[DF_Lig_All['Atom ID'] == r['LigID']]
                if not lig_row.empty:
                    lig_xyz = lig_row[['X', 'Y', 'Z']].astype(float).values[0]

        if coord_source == 'receptor':
            xyz = rec_xyz
        elif coord_source == 'ligand':
            xyz = lig_xyz
        else:
            xyz = (rec_xyz + lig_xyz) / 2 if rec_xyz is not None and lig_xyz is not None else None

        if xyz is None:
            xs.append(np.nan); ys.append(np.nan); zs.append(np.nan)
        else:
            xs.append(round(float(xyz[0]), 3))
            ys.append(round(float(xyz[1]), 3))
            zs.append(round(float(xyz[2]), 3))

    DF = DF.copy()
    DF['X'], DF['Y'], DF['Z'] = xs, ys, zs
    # Reordenar: X,Y,Z entre 'Angle' e 'Interaction' (en vez de al final).
    cols = [c for c in DF.columns if c not in ('X', 'Y', 'Z')]
    insert_at = cols.index('Interaction')
    cols[insert_at:insert_at] = ['X', 'Y', 'Z']
    return DF[cols]


# ──────────────────────────────────────────────────────────────────────────────
# Resumen en consola
# ──────────────────────────────────────────────────────────────────────────────

_TYPE_LABELS = {
    'acceptor':    'H-bond lig→acceptor',
    'donor':       'H-bond lig→donor',
    'aromatic':    'Aromatic',
    'hydrophobic': 'Hydrophobic',
    'salt_bridge': 'Salt bridge',
    'pi_cation':   'π-cation',
}


def _append_cumulative_csv(row, filename):
    """Agrega una fila (Receptor+Ligand identifican el par) a un CSV acumulado
    en el directorio actual. Escribe el header solo si el archivo no existe."""
    path = Path(filename)
    pd.DataFrame([row]).to_csv(path, mode='a', header=not path.exists(), index=False)


def print_summary(receptor, ligand, DF_validated, df_pockets_summary=None):
    bar = '═' * 76
    print(f'\n{bar}')
    print(f'  Receptor : {Path(receptor).stem}')
    print(f'  Ligand   : {Path(ligand).stem}')
    if df_pockets_summary is not None:
        n_pockets = int((df_pockets_summary['Is_Pocket'] == 'Yes').sum()) if not df_pockets_summary.empty else 0
        print(f'  Hydrophobic pockets     : {n_pockets}')
    print(f'  Validated interactions  : {len(DF_validated)}')
    if DF_validated.empty:
        print('  (none)')
    else:
        for t, n in DF_validated['Type'].value_counts().items():
            print(f'    {_TYPE_LABELS.get(t, t):<22}: {n}')
        print(f'  {"─"*74}')
        print(f'  {"Type":<22} {"Residue":>9}  {"Atom":<12} {"Dist":>6}  {"Angle":>7}  {"Lig"}')
        print(f'  {"─"*22} {"─"*9}  {"─"*12} {"─"*6}  {"─"*7}  {"─"*18}')
        for _, row in DF_validated.iterrows():
            label   = _TYPE_LABELS.get(row['Type'], row['Type'])
            res_str = f"{row['Res']}{int(row['Pos R'])}"
            ang_str = f"{float(row['Angle']):.1f}°" if float(row['Angle']) != 0 else '  —'
            lig_str = str(row['Lig'])
            print(f"  {label:<22} {res_str:>9}  {str(row['Atom']):<12} {row['Dist']:>5.2f}Å  {ang_str:>7}  {lig_str}")
    print(f'{bar}\n')


def _validate_interaction_angles(interactions, ligand_points, ligand_atoms,
                                 active_site, acceptor_antecedents):
    """Fill the angle column for directional and aromatic interaction rows."""
    for row in range(interactions.shape[0]):
        kind = interactions.iloc[row, 5]
        if kind == 'acceptor':
            antecedent = Busqueda_Antecesor_Lig(interactions.iloc[row, 8], ligand_atoms)
            acceptor = np.array(ligand_points[
                ligand_points['Atom ID'] == interactions.iloc[row, 8]].iloc[0, [1, 2, 3]])
            receptor_atom = active_site[
                (active_site['Pos'] == interactions.iloc[row, 0])
                & (active_site['Atom'] == interactions.iloc[row, 2])]
            donor = np.array(receptor_atom[['X', 'Y', 'Z']]).reshape(-1)
            interactions.iloc[row, 6] = float(angle_three_points(donor, acceptor, antecedent))
        elif kind == 'donor':
            donor = np.array(ligand_points[
                ligand_points['Atom ID'] == interactions.iloc[row, 8]].iloc[0, [1, 2, 3]])
            receptor_atom = active_site[
                (active_site['Pos'] == interactions.iloc[row, 0])
                & (active_site['Atom'] == interactions.iloc[row, 2])]
            interactions.iloc[row, 6] = (
                receptor_acceptor_angle(donor, receptor_atom.iloc[0], active_site,
                                        acceptor_antecedents)
                if len(receptor_atom) == 1 else np.nan)
        elif kind == 'aromatic':
            receptor_ring = active_site[active_site['Pos'] == interactions.iloc[row, 0]]
            ligand_ring = ligand_points[ligand_points['Caso'] == interactions.iloc[row, 4]]
            interactions.iloc[row, 6] = Interaccion_Aromatica(receptor_ring, ligand_ring)
    return interactions


def _classify_interactions(interactions, cfg):
    """Set interaction status and explain the criteria for each candidate."""
    reasons = []
    for row in range(interactions.shape[0]):
        kind = interactions.iloc[row, 5]
        distance = float(interactions.iloc[row, 3])
        angle = float(interactions.iloc[row, 6])
        if kind in ('hydrophobic', 'salt_bridge', 'pi_cation'):
            reason = 'distance_pass'  # validated on raw distances by their detectors
        elif kind in ('acceptor', 'donor'):
            status, reason = evaluate_hbond(distance, angle, cfg)
            interactions.iloc[row, 7] = status
        elif kind == 'aromatic':
            status, reason = evaluate_aromatic(distance, angle, cfg)
            interactions.iloc[row, 7] = status
        else:
            raise ValueError(f'Tipo de interacción sin criterio: {kind}')
        reasons.append(reason)
    interactions['Reason'] = pd.Series(reasons, index=interactions.index, dtype='object')
    return interactions


# ──────────────────────────────────────────────────────────────────────────────
# Análisis de un par receptor-ligando
# ──────────────────────────────────────────────────────────────────────────────

_ALL_TYPES = ['acceptor', 'donor', 'aromatic', 'hydrophobic', 'salt_bridge', 'pi_cation']


@record_analysis('ligand')
def analyze_pair(receptor_pdb, Ligand_imput, chain_receptor, cfg, *, ligand_reference=None,
                 legacy_rings=False, run_context=None):
    """Ejecuta el pipeline completo para un par receptor-ligando."""
    ligand_plot           = cfg['ligand_plot']
    vmd_output            = cfg['vmd_output']
    cumulative_output     = cfg['cumulative_output']
    Interaction_Coord_Source = cfg['Interaction_Coord_Source']
    Volume_Plot           = cfg['Volume_Plot']
    Bias                  = cfg['Bias']
    Bias_Validated_Only   = cfg['Bias_Validated_Only']
    Distances_Aromatic    = cfg['Distances_Aromatic']
    Distancia_Hidrofobica = cfg['Distancia_Hidrofobica']
    Distancia_Centro_Activo = cfg['Distancia_Centro_Activo']
    Ring_Planarity_RMSD_Max = cfg['Ring_Planarity_RMSD_Max']
    Pocket_Min_Residues   = cfg['Pocket_Min_Residues']
    Pocket_Coverage_Threshold = cfg['Pocket_Coverage_Threshold']
    Pocket_Density_Radius = cfg['Pocket_Density_Radius']
    Aceptores_Prot        = cfg['Aceptores_Prot']
    Dadores_Prot          = cfg['Dadores_Prot']
    Aceptot_antecedent    = cfg['Aceptot_antecedent']

    threshold_PH = hbond_search_cutoff(cfg)
    numero_anillo_aromatico = 5 if legacy_rings else 4

    receptor = Path(receptor_pdb).stem
    ligand   = Path(Ligand_imput).stem
    folder   = f'{receptor}_{ligand}'
    Path(folder).mkdir(exist_ok=True)

    # ── Limpieza del ligando ──────────────────────────────────────
    remove_bias(Ligand_imput, folder)

    # ── Ligando: hot-points ───────────────────────────────────────
    mol = Chem.MolFromPDBFile(Ligand_imput, removeHs=False,
                            sanitize=ligand_reference is None)
    if mol is None:
        raise ValueError(f"RDKit no pudo leer el ligando: {Ligand_imput}")

    if ligand_reference is not None:
        mol = apply_reference(mol, ligand_reference)

    pdb_coords, CM = extract_coords_from_pdb(Ligand_imput)
    if ligand_reference is not None and (
            len(pdb_coords) != mol.GetNumAtoms()
            or any(atom.GetPDBResidueInfo() is None
                   or atom.GetPDBResidueInfo().GetSerialNumber() != pdb_coords[i][0]
                   for i, atom in enumerate(mol.GetAtoms()))):
        raise ValueError("Los átomos leídos por RDKit no corresponden a las filas del PDB.")
    acceptor_atoms, donor_atoms = search_hot_points(
        Ligand_imput, mol, pdb_coords, ligand_plot, folder,
        reference_chemistry=ligand_reference is not None)

    # ── Anillos aromáticos ────────────────────────────────────────
    aromatic_rings_data, rings_data = search_rings(mol, pdb_coords, numero_anillo_aromatico,
                                                    Ring_Planarity_RMSD_Max,
                                                    reference_chemistry=ligand_reference is not None)
    if ligand_plot == 'Yes':
        visualize_rings(mol, aromatic_rings_data, Ligand_imput, folder)
    DF_Aro = pd.DataFrame(rings_data, columns=['Átomo', 'Coord X', 'Coord Y', 'Coord Z', 'Caso'])

    # ── DataFrame de puntos del ligando ──────────────────────────
    # acceptor_atoms/donor_atoms ya son índices únicos dentro de pdb_coords: se usan
    # directamente (en vez de rebuscar por nombre de átomo, que puede repetirse).
    coordenadas = []
    for idx in acceptor_atoms:
        p = pdb_coords[idx]
        coordenadas.append([p[1], p[5], p[6], p[7], 'acceptor', p[0]])
    for idx in donor_atoms:
        p = pdb_coords[idx]
        coordenadas.append([p[1], p[5], p[6], p[7], 'donor', p[0]])
    DF_Lig = pd.DataFrame(coordenadas,
                          columns=['Átomo', 'Coord X', 'Coord Y', 'Coord Z', 'Caso', 'Atom ID'])
    # Frames vacíos (sin aceptores/donores o sin anillos) no tienen dtypes declarados
    # (mismo caso que la concatenación de tipos de interacción más abajo): concatenarlos
    # igual dispara el FutureWarning de pandas sobre inferencia de dtype en columnas
    # vacías/all-NA. Se excluyen antes de concatenar.
    lig_frames = [f for f in (DF_Lig, DF_Aro) if not f.empty]
    DF_Lig = pd.concat(lig_frames, ignore_index=True) if lig_frames else DF_Lig

    # ── Receptor: sitio activo ────────────────────────────────────
    pdb_parser = PDBParser(QUIET=True)
    structure  = pdb_parser.get_structure('pdb', receptor_pdb)
    DF_Active_Site = active_site_residues(structure, CM, chain_receptor,
                                           Distancia_Centro_Activo, ligand)
    receptor_points = Coordenadas_interes_receptor(Aceptores_Prot, Dadores_Prot, DF_Active_Site)

    # Volumen del sitio activo completo (envolvente convexa de todos sus átomos),
    # independiente del volumen por pocket hidrofóbico calculado más abajo.
    Site_Volume, site_hull = convex_hull_volume(DF_Active_Site[['X', 'Y', 'Z']].values.astype(float))

    # ── DataFrame de interacciones ────────────────────────────────
    DF_Interacciones = pd.DataFrame({c: pd.Series(dtype=t) for c, t in [
        ('Pos R', 'int'), ('Res', 'object'), ('Atom', 'object'), ('Dist', 'float64'),
        ('Lig', 'object'), ('Type', 'object'), ('Angle', 'float64'), ('Interaction', 'object'),
        ('LigID', 'float64')]})

    DF_Interacciones = residuos_contacto('Dador',   'acceptor', receptor_points,
                                          DF_Lig, DF_Interacciones, threshold_PH)
    DF_Interacciones = residuos_contacto('Aceptor', 'donor',    receptor_points,
                                          DF_Lig, DF_Interacciones, threshold_PH)

    # Aromáticas
    aromatic_lig_df  = DF_Lig.query('Caso.str.contains("aromatic")', engine='python')
    Sub_Set_Receptor = receptor_points.query('Type == "aromatic"')
    for cas in aromatic_lig_df['Caso'].unique():
        Sub_Set_Ligando = aromatic_lig_df.query('Caso == @cas')
        ring_center     = np.mean(np.array(Sub_Set_Ligando.iloc[:, [1,2,3]]), axis=0)
        Matriz_receptor = np.array(Sub_Set_Receptor.iloc[:, [4,5,6]])
        for idx, distance in neighbors(ring_center, Matriz_receptor, Distances_Aromatic):
            closest = Sub_Set_Receptor.iloc[idx]
            DF_Interacciones.loc[len(DF_Interacciones)] = [
                closest.iloc[1], closest.iloc[2], closest.iloc[3],
                distance, Sub_Set_Ligando.iloc[0, 4], 'aromatic', 0.0, 0, np.nan]

    DF_Lig_All = generate_df_ligand(pdb_coords)
    DF_Interacciones = DF_Interacciones.drop_duplicates()

    # ── Nuevos tipos de interacción ───────────────────────────────
    df_hpho = search_hydrophobic(mol, pdb_coords, DF_Active_Site, Distancia_Hidrofobica)
    df_salt = search_salt_bridges(mol, pdb_coords, DF_Active_Site,
                                 distance_cutoff('salt_bridge', cfg))
    df_pica = search_pi_cation(DF_Active_Site, aromatic_lig_df,
                              distance_cutoff('pi_cation', cfg))
    df_pocket_summary, df_pocket_detail, pocket_hulls = search_hydrophobic_pockets(
        mol, pdb_coords, DF_Active_Site, Distancia_Hidrofobica,
        Pocket_Min_Residues, Pocket_Coverage_Threshold, Pocket_Density_Radius)
    # Frames vacíos (sin matches) no tienen dtypes declarados (columns=_DF_COLS sin data);
    # concatenarlos junto con DF_Interacciones (tipado) dispara el FutureWarning de pandas
    # sobre inferencia de dtype sobre columnas vacías/all-NA. Se excluyen antes de concatenar.
    frames = [f for f in (DF_Interacciones, df_hpho, df_salt, df_pica) if not f.empty]
    if frames:
        DF_Interacciones = pd.concat(frames, ignore_index=True).drop_duplicates()

    # ── Validación geométrica y clasificación ─────────────────────
    DF_Interacciones = _validate_interaction_angles(
        DF_Interacciones, DF_Lig, DF_Lig_All, DF_Active_Site, Aceptot_antecedent)
    DF_Interacciones = _classify_interactions(DF_Interacciones, cfg)

    DF_Interacciones = DF_Interacciones.drop_duplicates()

    # ── Coordenadas X,Y,Z (receptor/ligand/center según config) ────
    # Se agrega antes de partir en all/threshold/true para que las tres salidas
    # compartan exactamente las mismas columnas. Necesita LigID (se descarta recién
    # al escribir cada CSV, más abajo).
    DF_Interacciones = add_interaction_coords(DF_Interacciones, receptor_points, DF_Lig,
                                              DF_Lig_All, DF_Active_Site, Interaction_Coord_Source)

    # ── Salidas CSV ───────────────────────────────────────────────
    # LigID (serial de átomo, uso interno) se excluye de los CSV; X/Y/Z sí quedan.
    DF_Interacciones.drop(columns=['LigID']).to_csv(f'{folder}/Interaction_{receptor}_{ligand}_all.csv')
    # Distance-only types were already filtered on raw distances in their
    # detectors: do not re-filter rounded/averaged output values at the boundary.
    distance_mask = DF_Interacciones['Type'].isin(['hydrophobic', 'salt_bridge', 'pi_cation'])
    for kind in ('acceptor', 'donor', 'aromatic'):
        distance_mask |= ((DF_Interacciones['Type'] == kind)
                          & within_distance(DF_Interacciones['Dist'], distance_cutoff(kind, cfg)))
    DF_dist = DF_Interacciones[distance_mask]
    DF_dist.drop(columns=['LigID']).to_csv(f'{folder}/Interaction_{receptor}_{ligand}_threshold.csv')
    DF_true = DF_Interacciones[DF_Interacciones['Interaction'] == 'Yes']
    DF_true.drop(columns=['LigID']).to_csv(f'{folder}/Interaction_{receptor}_{ligand}_true.csv')

    if Bias == 'Yes':
        bias_df_true = DF_true if Bias_Validated_Only == 'Yes' else None
        export_bpf(DF_Lig, f'{folder}/{receptor}_{ligand}.bpf', bias_df_true)
        export_bpf_pdb(DF_Lig, f'{folder}/{receptor}_{ligand}_bias.pdb', bias_df_true)

    df_pocket_summary.to_csv(f'{folder}/Pockets_{receptor}_{ligand}.csv', index=False)

    df_pocket_chi = compute_pocket_chi_angles(df_pocket_summary, DF_Active_Site)
    df_pocket_chi.to_csv(f'{folder}/Pockets_{receptor}_{ligand}_chi.csv', index=False)

    df_site_chi = compute_active_site_chi_angles(DF_Active_Site)
    df_site_chi.to_csv(f'{folder}/ActiveSite_{receptor}_{ligand}_chi.csv', index=False)

    df_phi_psi = compute_active_site_phi_psi(structure, chain_receptor, DF_Active_Site)
    df_phi_psi.to_csv(f'{folder}/ActiveSite_{receptor}_{ligand}_ramachandran.csv', index=False)

    if Volume_Plot == 'Yes':
        plot_ramachandran(df_phi_psi, f'Ramachandran — active site ({receptor}/{ligand})',
                          f'{folder}/ActiveSite_{receptor}_{ligand}_ramachandran.png')
        for chi_name in ('chi1', 'chi2', 'chi3', 'chi4', 'chi5'):
            if df_site_chi[chi_name].notna().any():
                plot_chi_profile(df_site_chi, chi_name,
                                 f'{chi_name} — active site ({receptor}/{ligand})',
                                 f'{folder}/ActiveSite_{receptor}_{ligand}_{chi_name}.png')
        if site_hull is not None:
            site_points = DF_Active_Site[['X', 'Y', 'Z']].values.astype(float)
            site_title  = f'Active site — {Site_Volume:.1f} Å³'
            plot_hull_volume(site_points, site_hull, site_title,
                             f'{folder}/ActiveSite_{receptor}_{ligand}_volume.png')
            plot_hull_surface(site_points, site_hull, site_title,
                              f'{folder}/ActiveSite_{receptor}_{ligand}_volume_solid.png')
        if not df_pocket_summary.empty:
            qualifying = df_pocket_summary[df_pocket_summary['Is_Pocket'] == 'Yes']
            for _, prow in qualifying.iterrows():
                hull_data = pocket_hulls.get(prow['Pocket'])
                if hull_data is None:
                    continue  # volumen no calculable (< 4 puntos o geometría degenerada)
                points, hull = hull_data
                pocket_title = f"Pocket {prow['Pocket']} — {prow['Volume_A3']:.1f} Å³"
                plot_hull_volume(points, hull, pocket_title,
                                 f"{folder}/Pocket_{prow['Pocket']}_{receptor}_{ligand}_volume.png")
                plot_hull_surface(points, hull, pocket_title,
                                  f"{folder}/Pocket_{prow['Pocket']}_{receptor}_{ligand}_volume_solid.png")

    shutil.copy(Ligand_imput, f'{folder}/{Path(Ligand_imput).name}')
    shutil.copy(receptor_pdb, f'{folder}/{Path(receptor_pdb).name}')

    if vmd_output == 'Yes':
        scripting_vmd(DF_true, receptor_points, aromatic_lig_df, DF_Lig,
                      receptor_pdb, chain_receptor, Ligand_imput, folder)
        scripting_vmd_hydrophobic(DF_true, DF_Active_Site, DF_Lig_All,
                                  receptor_pdb, chain_receptor, Ligand_imput, folder)
        scripting_vmd_pockets(df_pocket_detail, receptor_pdb, chain_receptor, Ligand_imput, folder)
        scripting_vmd_combined(DF_true, receptor_points, DF_Lig, df_pocket_detail,
                               receptor_pdb, chain_receptor, Ligand_imput, folder)

    print_summary(receptor_pdb, Ligand_imput, DF_true, df_pocket_summary)

    # ── Resumen del par (dentro de la carpeta) ────────────────────
    counts_dist = dict(DF_dist['Type'].value_counts())
    counts_true = dict(DF_true['Type'].value_counts())
    dat = {'Receptor': receptor, 'Ligand': ligand,
           'Total_all': len(DF_Interacciones), 'Total_dist': len(DF_dist),
           'Total_true': len(DF_true)}
    for t in _ALL_TYPES:
        dat[f'dist_{t}'] = counts_dist.get(t, 0)
        dat[f'true_{t}'] = counts_true.get(t, 0)
    dat['pocket_hydrophobic'] = int((df_pocket_summary['Is_Pocket'] == 'Yes').sum()) \
        if not df_pocket_summary.empty else 0
    dat['ActiveSite_Volume_A3'] = Site_Volume
    pd.DataFrame([dat]).to_csv(f'{folder}/summary.csv', index=False)

    pd.DataFrame([{'Receptor': receptor, 'Ligand': ligand,
                   'CM X': CM[0], 'CM Y': CM[1], 'CM Z': CM[2]}]).to_csv(
        f'{folder}/CM.csv', index=False)

    # ── Acumulados globales (fuera de la carpeta del par) ─────────
    # Cada fila queda identificada por Receptor+Ligand, así que corridas
    # sucesivas o en batch no se pisan entre sí.
    if cumulative_output == 'Yes':
        _append_cumulative_csv(dat, 'Interactions_close.csv')
        _append_cumulative_csv({'Receptor': receptor, 'Ligand': ligand,
                                 'CM X': CM[0], 'CM Y': CM[1], 'CM Z': CM[2]}, 'CM_all.csv')


# ──────────────────────────────────────────────────────────────────────────────
# Modo especial: bias points del receptor alrededor de una coordenada
# ──────────────────────────────────────────────────────────────────────────────

@record_analysis('site_bias')
def analyze_site_bias(receptor_pdb, chain_receptor, point, radius, cfg, method='atom', *,
                      run_context=None):
    """Modo especial sin ligando: dada una coordenada arbitraria `point`
    (x, y, z) -- por ejemplo, donde se espera que caiga un ligando -- busca
    todos los residuos del receptor dentro de `radius` Å (reutilizando
    active_site_residues(), igual que analyze_pair() pero centrado en `point`
    en vez del centro de masa de un ligando real), y exporta sus puntos de
    bias como .bpf (GOLD) + PDB dummy para VMD. No usa nada del ligando: es
    exclusivamente para cuando no hay una pose de ligando disponible.

    method:
      'atom'  (default) -- un punto de bias por átomo/anillo del receptor,
              en su propia coordenada (Coordenadas_interes_receptor()).
      'ideal' -- abanico de puntos "ideales": geometría de H-bond (distancia +
              ángulo + diedro) hacia donde debería caer el átomo COMPLEMENTARIO
              del ligando, más el abanico de posiciones de stacking para
              anillos aromáticos (ideal_sites.py, puerto de
              ideal_interaction_sites.py). Requiere un PDB del receptor
              protonado (nomenclatura Maestro/Amber) para los grupos donores;
              los aceptores basados en átomos pesados y los aromáticos no lo
              necesitan.

    Útil para generar bias points del receptor antes de tener una pose de
    ligando (ej. para guiar un docking ciego con GOLD en una región conocida
    del sitio activo)."""
    if method not in ('atom', 'ideal'):
        raise ValueError(f"method inválido: {method!r} (usar 'atom' o 'ideal')")

    receptor = Path(receptor_pdb).stem
    px, py, pz = point
    folder = f'{receptor}_site_{px:.1f}_{py:.1f}_{pz:.1f}'
    Path(folder).mkdir(exist_ok=True)

    pdb_parser = PDBParser(QUIET=True)
    structure  = pdb_parser.get_structure('pdb', receptor_pdb)
    # lig='' -> no se excluye ningún residuo por nombre (no hay ligando real acá).
    DF_Active_Site = active_site_residues(structure, list(point), chain_receptor, radius, '')

    suffix   = '' if method == 'atom' else '_ideal'
    bpf_path = f'{folder}/{receptor}_site{suffix}.bpf'
    pdb_path = f'{folder}/{receptor}_site{suffix}_bias.pdb'

    if method == 'atom':
        receptor_points = Coordenadas_interes_receptor(
            cfg['Aceptores_Prot'], cfg['Dadores_Prot'], DF_Active_Site)
        export_bpf_receptor(receptor_points, bpf_path)
        export_bpf_pdb_receptor(receptor_points, pdb_path)
        n_points = len(receptor_points)
        counts   = receptor_points['Type'].value_counts() if n_points else {}
        result   = receptor_points
    else:
        resids = DF_Active_Site['Pos'].unique().tolist()
        rows   = ideal_site_points(structure, chain_receptor, resids)
        export_bpf_points(rows, bpf_path)
        export_bpf_pdb_points(rows, pdb_path)
        n_points = len(rows)
        counts   = pd.Series([r[3] for r in rows]).value_counts() if rows else {}
        result   = rows

    shutil.copy(receptor_pdb, f'{folder}/{Path(receptor_pdb).name}')

    print(f"\n  Site bias [{method}] @ ({px:.2f}, {py:.2f}, {pz:.2f}), radius {radius} Å")
    if n_points == 0:
        print(f"  [WARN] No bias points found within {radius} Å.")
    else:
        for t, n in counts.items():
            print(f"    {t:<10}: {n}")
    print(f"  -> {bpf_path}")
    print(f"  -> {pdb_path}\n")

    return result


# ──────────────────────────────────────────────────────────────────────────────
# Modo sondeo: interacciones simuladas de puntos arbitrarios
# ──────────────────────────────────────────────────────────────────────────────

_PROBE_TYPE_LABELS = {**_TYPE_LABELS,
                      'acceptor': 'H-bond (probe acceptor)',
                      'donor':    'H-bond (probe donor)'}


@record_analysis('probe')
def analyze_probe(receptor_pdb, chain_receptor, points, cfg, label=None, *, run_context=None):
    """Modo sondeo sin ligando: para cada punto (x, y, z, tipos) simula qué
    interacciones haría un átomo/grupo del ligando de esos tipos ubicado ahí
    (ver probe.py). El sitio activo se arma una sola vez alrededor del
    centroide de todos los puntos, con radio centroid_distance + la distancia
    del punto más alejado, así cubre cada punto con el mismo margen que
    analyze_pair() usa alrededor de un ligando real.

    Escribe en <receptor>_probe[_label]/: Probe_<rec>_all.csv (todos los
    contactos candidatos, incluidos choques), Probe_<rec>_true.csv (solo
    Interaction == 'Yes'), <rec>_probe_points.pdb (los puntos como átomos
    dummy) y, si vmd_output == 'Yes', vmd_probe_<rec>.tcl."""
    receptor = Path(receptor_pdb).stem
    folder = f'{receptor}_probe' + (f'_{label}' if label else '')
    Path(folder).mkdir(exist_ok=True)

    coords = np.array([p[:3] for p in points], dtype=float)
    center = coords.mean(axis=0)
    radius = cfg['Distancia_Centro_Activo'] + float(np.linalg.norm(coords - center, axis=1).max())

    pdb_parser = PDBParser(QUIET=True)
    structure  = pdb_parser.get_structure('pdb', receptor_pdb)
    # lig='' -> no se excluye ningún residuo por nombre (no hay ligando real).
    DF_Active_Site  = active_site_residues(structure, list(center), chain_receptor, radius, '')
    receptor_points = Coordenadas_interes_receptor(cfg['Aceptores_Prot'], cfg['Dadores_Prot'],
                                                   DF_Active_Site)

    df = probe_interactions(points, DF_Active_Site, receptor_points, cfg)
    df_true = df[df['Interaction'] == 'Yes']
    df.to_csv(f'{folder}/Probe_{receptor}_all.csv', index=False)
    df_true.to_csv(f'{folder}/Probe_{receptor}_true.csv', index=False)

    probe_pdb = f'{folder}/{receptor}_probe_points.pdb'
    write_probe_pdb(points, probe_pdb)
    shutil.copy(receptor_pdb, f'{folder}/{Path(receptor_pdb).name}')
    if cfg['vmd_output'] == 'Yes':
        scripting_vmd_probe(df_true, receptor_pdb, probe_pdb, chain_receptor, folder)

    bar = '═' * 76
    print(f'\n{bar}')
    print(f'  Receptor : {receptor}   ({len(points)} probe point(s))')
    for n, (x, y, z, types) in enumerate(points, start=1):
        sub   = df[df['Probe'] == n]
        yes   = sub[sub['Interaction'] == 'Yes']
        clash = sub[sub['Interaction'] == 'Clash']
        tipos = 'all' if set(types) == set(PROBE_TYPES) else ','.join(types)
        print(f'  {"─"*74}')
        print(f'  Probe {n} @ ({x:.2f}, {y:.2f}, {z:.2f})  as: {tipos}')
        if not clash.empty:
            res = ', '.join(sorted({f"{r['Res']}{r['Pos R']}:{r['Atom']}"
                                    for _, r in clash.iterrows()}))
            print(f"    [WARN] steric clash (< {distance_cutoff('clash', cfg):g} Å): {res}")
        if yes.empty:
            print('    (no interactions)')
        for _, row in yes.iterrows():
            lbl     = _PROBE_TYPE_LABELS.get(row['Type'], row['Type'])
            res_str = f"{row['Res']}{int(row['Pos R'])}"
            ang_str = f"{row['Angle']:.1f}°" if pd.notna(row['Angle']) else '  —'
            print(f"    {lbl:<24} {res_str:>9}  {str(row['Atom']):<12} "
                  f"{row['Dist']:>5.2f}Å  {ang_str:>7}")
    print(f'{bar}')
    print(f'  -> {folder}/Probe_{receptor}_all.csv')
    print(f'  -> {folder}/Probe_{receptor}_true.csv\n')

    return df
