"""Modo sondeo (probe): simula qué interacciones haría un átomo/grupo del
ligando si estuviera en una coordenada dada, sin necesidad de un ligando real.

Cada punto de sondeo se evalúa con uno o varios "roles" (PROBE_TYPES): como
aceptor, donor, anillo aromático, carbono apolar, catión o anión. Para cada rol
se buscan los átomos/grupos complementarios del receptor con los mismos
criterios de distancia (y de ángulo cuando el punto alcanza para definirlo) que
usa analyze_pair() para un ligando real:

  acceptor    -> donores del receptor (YAML `donors`). Si el donor es un H
                 explícito, se valida también el ángulo D-H···sonda en el H.
  donor       -> aceptores del receptor (YAML `acceptors`), ángulo
                 sonda-Aceptor-Antecedente, igual que en analyze_pair().
  aromatic    -> centroides de TYR/PHE/TRP (Distances_Aromatic) y cationes
                 ARG/LYS/HIS (π-catión). Un punto no define un plano, así que
                 no hay ángulo entre anillos: 'Angle' es el ángulo entre la
                 normal del anillo del receptor y el vector centroide->sonda
                 (0° = sonda sobre la cara del anillo, 90° = en su plano).
                 Validado solo por distancia.
  hydrophobic -> átomos apolares del receptor (_HYDROPHOBIC_ATOMS), colapsados
                 por residuo como en search_hydrophobic().
  cation      -> carboxilatos ASP/GLU (puente salino) y anillos aromáticos
                 del receptor (π-catión).
  anion       -> ARG/LYS/HIP (puente salino).

Además, para todo punto se reportan los átomos pesados del receptor a menos de
Probe_Clash_Distance Å (Type 'clash', Interaction 'Clash'): una sonda ahí caería dentro
de la proteína y sus "interacciones" no son realistas."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from interactions_search.contacts import (
    _CATION_REC_ATOMS,
    _HYDROPHOBIC_ATOMS,
    _SALT_NEG_ATOMS,
    _SALT_POS_ATOMS,
)
from interactions_search.geometry import angle_three_points
from interactions_search.interaction_rules import (
    atoms_near,
    distance_cutoff,
    evaluate_aromatic,
    evaluate_hbond,
    hbond_search_cutoff,
    neighbors,
    receptor_acceptor_angle,
)

__all__ = ["PROBE_TYPES", "read_probe_file", "probe_interactions", "write_probe_pdb"]

PROBE_TYPES = ('acceptor', 'donor', 'aromatic', 'hydrophobic', 'cation', 'anion')

_PROBE_COLS = ['Probe', 'Probe_Type', 'Probe_X', 'Probe_Y', 'Probe_Z',
               'Pos R', 'Res', 'Atom', 'Dist', 'Type', 'Angle', 'Interaction',
               'Rec_X', 'Rec_Y', 'Rec_Z', 'Reason']

# Mismos 6 átomos cuyo centroide calcula get_aromatic_coord() (el anillo de 6 en TRP).
_AROMATIC_RING_ATOMS = {
    'PHE': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'),
    'TYR': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'),
    'TRP': ('CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'),
}

# Tipos de los archivos que genera el propio paquete (.bpf y PDB dummy de bias.py):
# 'don' = ahí iría un DONOR del ligando -> se sondea como donor, etc.
_BPF_TO_PROBE = {'don': 'donor', 'acc': 'acceptor', 'aro': 'aromatic'}
_PDB_RESNAME_TO_PROBE = {'DON': 'donor', 'ACC': 'acceptor', 'ARO': 'aromatic',
                         'HPH': 'hydrophobic', 'CAT': 'cation', 'ANI': 'anion'}
_PROBE_TO_PDB_RESNAME = {v: k for k, v in _PDB_RESNAME_TO_PROBE.items()}


def _parse_probe_type(token):
    token = token.strip().lower()
    if token in _BPF_TO_PROBE:
        return _BPF_TO_PROBE[token]
    return token if token in PROBE_TYPES else None


def read_probe_file(path):
    """Lee puntos de sondeo de un archivo. Devuelve [(x, y, z, tipo_o_None)].

    - .pdb: coordenadas de cada ATOM/HETATM; el resname DON/ACC/ARO/HPH/CAT/ANI
      (ver bias.py / write_probe_pdb) fija el tipo, cualquier otro -> None.
    - cualquier otro (.bpf, .csv, .txt): una fila por punto 'x y z [...] [tipo]'
      separada por espacios, tabs, comas o ';'. Filas que no empiezan con 3
      números (headers, comentarios) se ignoran. La última columna, si es un
      tipo conocido (don/acc/aro de un .bpf, o un nombre de PROBE_TYPES), fija
      el tipo del punto.
    Los puntos con tipo None se sondean con los tipos pasados por --probe-type."""
    path = Path(path)
    points = []
    with open(path) as fh:
        if path.suffix.lower() == '.pdb':
            for line in fh:
                if line.startswith(('ATOM', 'HETATM')):
                    x, y, z = (float(line[i:i + 8]) for i in (30, 38, 46))
                    points.append((x, y, z, _PDB_RESNAME_TO_PROBE.get(line[17:20].strip())))
            return points
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            tokens = [t for t in re.split(r'[\s,;]+', line) if t]
            try:
                x, y, z = (float(t) for t in tokens[:3])
            except ValueError:
                continue
            ptype = _parse_probe_type(tokens[-1]) if len(tokens) > 3 else None
            points.append((x, y, z, ptype))
    return points


def _is_hydrogen(atom_name):
    # '1HH1'/'HH11'/'H' -> H. Ningún átomo pesado estándar de proteína empieza con H.
    return atom_name.lstrip('0123456789').startswith('H')


def _xyz(df):
    return df[['X', 'Y', 'Z']].to_numpy(dtype=float)


def _ring_normal(DF_Active_Site, pos, resname):
    """Normal unitaria del plano de mejor ajuste del anillo aromático del residuo
    (None si faltan átomos)."""
    ring = DF_Active_Site[(DF_Active_Site['Pos'] == pos) &
                          (DF_Active_Site['Atom'].isin(_AROMATIC_RING_ATOMS.get(resname, ())))]
    if len(ring) < 3:
        return None
    coords = _xyz(ring)
    _, _, vh = np.linalg.svd(coords - coords.mean(axis=0))
    return vh[-1]


def _face_angle(normal, center, p):
    """Ángulo (0-90°) entre la normal del anillo y el vector centroide->sonda."""
    if normal is None:
        return np.nan
    v = p - center
    if np.linalg.norm(v) == 0:
        return np.nan
    cos = abs(np.dot(normal, v)) / np.linalg.norm(v)
    return float(np.degrees(np.arccos(np.clip(cos, 0.0, 1.0))))


def _row(tipo, rec, dist, angle, interaction, rec_xyz, atom=None, *, reason='distance_pass'):
    return {'Pos R': int(rec['Pos']), 'Res': rec['Residue'],
            'Atom': atom if atom is not None else rec['Atom'],
            'Dist': round(float(dist), 3), 'Type': tipo,
            'Angle': round(float(angle), 2) if not np.isnan(angle) else np.nan,
            'Interaction': interaction,
            'Rec_X': round(float(rec_xyz[0]), 3), 'Rec_Y': round(float(rec_xyz[1]), 3),
            'Rec_Z': round(float(rec_xyz[2]), 3), 'Reason': reason}


def _probe_as_acceptor(p, receptor_points, DF_Active_Site, cfg):
    donors = receptor_points[receptor_points['Type'] == 'Dador']
    rows = []
    if donors.empty:
        return rows
    for k, distance in neighbors(p, _xyz(donors), hbond_search_cutoff(cfg)):
        rec = donors.iloc[k]
        d_xyz = _xyz(donors.iloc[[k]])[0]
        angle = np.nan
        if _is_hydrogen(rec['Atom']):
            # Ángulo D-H···A en el H: átomo pesado padre = el más cercano del residuo.
            res_atoms = DF_Active_Site[DF_Active_Site['Pos'] == rec['Pos']]
            heavy = res_atoms[~res_atoms['Atom'].map(_is_hydrogen)]
            if not heavy.empty:
                heavy_xyz = _xyz(heavy)
                parent = heavy_xyz[np.argmin(np.linalg.norm(heavy_xyz - d_xyz, axis=1))]
                angle = float(angle_three_points(parent, d_xyz, p))
        status, reason = evaluate_hbond(distance, angle, cfg,
                                        angle_required=_is_hydrogen(rec['Atom']))
        rows.append(_row('acceptor', rec, distance, angle, status, d_xyz, reason=reason))
    return rows


def _probe_as_donor(p, receptor_points, DF_Active_Site, cfg):
    acceptors = receptor_points[receptor_points['Type'] == 'Aceptor']
    rows = []
    if acceptors.empty:
        return rows
    for k, distance in neighbors(p, _xyz(acceptors), hbond_search_cutoff(cfg)):
        rec = acceptors.iloc[k]
        a_xyz = _xyz(acceptors.iloc[[k]])[0]
        angle = receptor_acceptor_angle(p, rec, DF_Active_Site, cfg['Aceptot_antecedent'])
        status, reason = evaluate_hbond(distance, angle, cfg)
        rows.append(_row('donor', rec, distance, angle, status, a_xyz, reason=reason))
    return rows


def _rings_near(p, receptor_points, DF_Active_Site, cutoff, tipo, cfg):
    rings = receptor_points[receptor_points['Type'] == 'aromatic']
    rows = []
    if rings.empty:
        return rows
    centers = _xyz(rings)
    for k, distance in neighbors(p, centers, cutoff):
        rec = rings.iloc[k]
        normal = _ring_normal(DF_Active_Site, rec['Pos'], rec['Residue'])
        angle = _face_angle(normal, centers[k], p)
        status, reason = (evaluate_aromatic(distance, angle, cfg, angle_required=False)
                          if tipo == 'aromatic' else ('Yes', 'distance_pass'))
        rows.append(_row(tipo, rec, distance, angle, status, centers[k], atom='center',
                         reason=reason))
    return rows


def _atoms_near(p, DF_Active_Site, atom_table, cutoff, tipo):
    return [_row(tipo, rec, distance, np.nan, 'Yes', xyz)
            for rec, distance, xyz in atoms_near(p, DF_Active_Site, atom_table, cutoff)]


def _collapse_by_residue(rows):
    """Como _collapse_same_residue_contacts: 1 fila por residuo, átomos unidos
    por coma, distancia y coordenada del receptor promediadas."""
    if not rows:
        return rows
    df = pd.DataFrame(rows)
    out = df.groupby('Pos R', as_index=False).agg(
        Res=('Res', 'first'), Atom=('Atom', lambda s: ','.join(sorted(s))),
        Dist=('Dist', 'mean'), Type=('Type', 'first'), Angle=('Angle', 'first'),
        Interaction=('Interaction', 'first'), Reason=('Reason', 'first'),
        Rec_X=('Rec_X', 'mean'), Rec_Y=('Rec_Y', 'mean'), Rec_Z=('Rec_Z', 'mean'))
    return out.round({'Dist': 3, 'Rec_X': 3, 'Rec_Y': 3, 'Rec_Z': 3}).to_dict('records')


def _clashes(p, DF_Active_Site, cfg):
    heavy = DF_Active_Site[~DF_Active_Site['Atom'].map(_is_hydrogen)]
    if heavy.empty:
        return []
    coords = _xyz(heavy)
    return [_row('clash', heavy.iloc[k], distance, np.nan, 'Clash', coords[k], reason='steric_clash')
            for k, distance in neighbors(p, coords, distance_cutoff('clash', cfg))]


def probe_interactions(points, DF_Active_Site, receptor_points, cfg):
    """points: [(x, y, z, [tipos])] con tipos ⊆ PROBE_TYPES. Devuelve un
    DataFrame (_PROBE_COLS) con una fila por contacto candidato; Interaction
    'Yes' si cumple los criterios, 'No' si no (H-bond fuera de distancia/ángulo
    final), 'Clash' para choques estéricos."""
    rows = []
    for n, (x, y, z, types) in enumerate(points, start=1):
        p = np.array([x, y, z], dtype=float)
        per_type = []
        for t in types:
            if t == 'acceptor':
                found = _probe_as_acceptor(p, receptor_points, DF_Active_Site, cfg)
            elif t == 'donor':
                found = _probe_as_donor(p, receptor_points, DF_Active_Site, cfg)
            elif t == 'aromatic':
                found = (_rings_near(p, receptor_points, DF_Active_Site,
                                     distance_cutoff('aromatic', cfg), 'aromatic', cfg)
                         + _atoms_near(p, DF_Active_Site, _CATION_REC_ATOMS,
                                       distance_cutoff('pi_cation', cfg), 'pi_cation'))
            elif t == 'hydrophobic':
                found = _collapse_by_residue(_atoms_near(
                    p, DF_Active_Site, _HYDROPHOBIC_ATOMS, distance_cutoff('hydrophobic', cfg),
                    'hydrophobic'))
            elif t == 'cation':
                found = (_atoms_near(p, DF_Active_Site, _SALT_NEG_ATOMS,
                                     distance_cutoff('salt_bridge', cfg), 'salt_bridge')
                         + _rings_near(p, receptor_points, DF_Active_Site,
                                       distance_cutoff('pi_cation', cfg), 'pi_cation', cfg))
            elif t == 'anion':
                found = _atoms_near(p, DF_Active_Site, _SALT_POS_ATOMS,
                                    distance_cutoff('salt_bridge', cfg), 'salt_bridge')
            else:
                raise ValueError(f"tipo de sonda inválido: {t!r} (usar {PROBE_TYPES})")
            per_type += [{**r, 'Probe_Type': t} for r in found]
        per_type += [{**r, 'Probe_Type': '-'} for r in _clashes(p, DF_Active_Site, cfg)]
        for r in per_type:
            r.update({'Probe': n, 'Probe_X': x, 'Probe_Y': y, 'Probe_Z': z})
        rows += per_type

    if not rows:
        return pd.DataFrame(columns=_PROBE_COLS)
    return pd.DataFrame(rows)[_PROBE_COLS]


def write_probe_pdb(points, filepath):
    """PDB dummy con un átomo por punto de sondeo (resname según tipo si el
    punto tiene uno solo -- ACC/DON/ARO/HPH/CAT/ANI --, PRB si se sondeó con
    varios), resid = número de sonda. Mismo formato que bias._write_bpf_pdb."""
    with open(filepath, 'w') as f:
        for i, (x, y, z, types) in enumerate(points, start=1):
            resname = _PROBE_TO_PDB_RESNAME.get(types[0], 'PRB') if len(types) == 1 else 'PRB'
            f.write(
                f"HETATM{i:>5} {'C':<4} {resname:>3} X{i:>4}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {'C':>2}\n"
            )
        f.write('END\n')
