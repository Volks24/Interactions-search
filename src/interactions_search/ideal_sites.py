"""Puntos de bias "ideales" del receptor: en vez de un punto por átomo en su
propia coordenada (ver receptor_site.py / bias.py:_collect_bias_points_receptor),
calcula dónde debería caer el átomo COMPLEMENTARIO del ligando -- geometría de
H-bond ideal (distancia + ángulo + diedro) -- y, para grupos con varias
posiciones de lone pair o de stacking aromático, el abanico completo de
posiciones en vez de una sola.

Puerto a Python 3 de `ideal_interaction_sites.py` (script standalone en la raíz
del repo, no integrado al paquete). Solo se usa del lado del receptor -- ver
`analyze_site_bias(..., method='ideal')` en pipeline.py -- nunca del lado del
ligando: `bias.py` sigue generando el bias del ligando a partir de sus
hot-points reales (search_hot_points/search_rings), donde "el punto ideal" no
aplica porque ya se conoce la posición real del átomo.

Requiere hidrógenos explícitos con nomenclatura Maestro/Amber (HNE, HH11,
HH12, HH21, HH22 para ARG; HD21/HD22 ASN; HE21/HE22 GLN; HG/HG1/HH SER/THR/TYR;
HE1 TRP; HZ1/HZ2/HZ3 LYS; HD1/HE2 HIS con protonación HIE/HID/HIP explícita)
para los grupos donores del receptor -- si faltan, esos puntos puntuales se
omiten en silencio. Los aceptores basados en átomos pesados (carbonilo,
carboxilato, amida, imidazol) y los anillos aromáticos no necesitan H."""
from __future__ import annotations

import math as mt

import numpy as np

__all__ = ["ideal_site_points"]

HBOND_DIST = 1.9    # Å, distancia ideal O/N···H
AROM_DIST = 3.8      # Å, distancia normal al plano del anillo (stacking apilado)
_PARALLEL_DISP = 1.5  # Å, desplazamiento paralelo (stacking parallel-displaced)

# type: 'acc' -> este átomo del receptor es aceptor, el punto generado es
#                donde iría un DONOR del ligando ('don' en el .bpf)
#       'don' -> este átomo del receptor es donor, el punto generado es
#                donde iría un ACEPTOR del ligando ('acc' en el .bpf)
_ROLE_TO_SITE_TYPE = {'acc': 'don', 'don': 'acc'}


def _sidechain_interactions():
    """[(resname, (atom1, atom2, atom3), (dist, angulo, diedro_o_alpha_beta), role)]
    atom1/atom2/atom3 definen el plano/eje de referencia (atom3 es el átomo de
    interés real, ej. OE1, ND2, HH11); diedro_o_alpha_beta es un float o
    'alpha'/'beta' para el caso especial de alcoholes (SER/THR/TYR), cuyos
    lone pairs se ubican relativos al diedro C-C-O-H real medido, no fijo."""
    return [
        # Aceptores del receptor (carboxilato/amida/imidazol) -> abanico de
        # 5 posiciones donor de ligando por lone pair, ideal_interactions()
        ('GLU', ('CG', 'CD', 'OE1'), (HBOND_DIST, 180, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE1'), (HBOND_DIST, 210, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE1'), (HBOND_DIST, 240, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE1'), (HBOND_DIST, 150, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE1'), (HBOND_DIST, 120, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE2'), (HBOND_DIST, 180, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE2'), (HBOND_DIST, 210, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE2'), (HBOND_DIST, 240, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE2'), (HBOND_DIST, 150, 0), 'acc'),
        ('GLU', ('CG', 'CD', 'OE2'), (HBOND_DIST, 120, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD1'), (HBOND_DIST, 180, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD1'), (HBOND_DIST, 210, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD1'), (HBOND_DIST, 240, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD1'), (HBOND_DIST, 150, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD1'), (HBOND_DIST, 120, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD2'), (HBOND_DIST, 180, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD2'), (HBOND_DIST, 210, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD2'), (HBOND_DIST, 240, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD2'), (HBOND_DIST, 150, 0), 'acc'),
        ('ASP', ('CB', 'CG', 'OD2'), (HBOND_DIST, 120, 0), 'acc'),
        ('ASN', ('CB', 'CG', 'OD1'), (HBOND_DIST, 180, 0), 'acc'),
        ('ASN', ('CB', 'CG', 'OD1'), (HBOND_DIST, 210, 0), 'acc'),
        ('ASN', ('CB', 'CG', 'OD1'), (HBOND_DIST, 240, 0), 'acc'),
        ('ASN', ('CB', 'CG', 'OD1'), (HBOND_DIST, 150, 0), 'acc'),
        ('ASN', ('CB', 'CG', 'OD1'), (HBOND_DIST, 120, 0), 'acc'),
        ('GLN', ('CG', 'CD', 'OE1'), (HBOND_DIST, 180, 0), 'acc'),
        ('GLN', ('CG', 'CD', 'OE1'), (HBOND_DIST, 210, 0), 'acc'),
        ('GLN', ('CG', 'CD', 'OE1'), (HBOND_DIST, 240, 0), 'acc'),
        ('GLN', ('CG', 'CD', 'OE1'), (HBOND_DIST, 150, 0), 'acc'),
        ('GLN', ('CG', 'CD', 'OE1'), (HBOND_DIST, 120, 0), 'acc'),
        # Alcoholes como aceptores: 2 lone pairs (alpha/beta), posición
        # relativa al diedro C-C-O-H real (requiere el H para orientarlos)
        ('SER', ('CA', 'CB', 'OG'), (HBOND_DIST, 109.5, 'alpha'), 'acc'),
        ('SER', ('CA', 'CB', 'OG'), (HBOND_DIST, 109.5, 'beta'), 'acc'),
        ('THR', ('CA', 'CB', 'OG1'), (HBOND_DIST, 109.5, 'alpha'), 'acc'),
        ('THR', ('CA', 'CB', 'OG1'), (HBOND_DIST, 109.5, 'beta'), 'acc'),
        ('TYR', ('CE1', 'CZ', 'OH'), (HBOND_DIST, 109.5, 'alpha'), 'acc'),
        ('TYR', ('CE1', 'CZ', 'OH'), (HBOND_DIST, 109.5, 'beta'), 'acc'),
        ('HID', ('CG', 'CD2', 'NE2'), (HBOND_DIST, 125.4, 180), 'acc'),
        ('HIE', ('NE2', 'CE1', 'ND1'), (HBOND_DIST, 125.4, 180), 'acc'),
        # Donores del receptor (requieren H explícito) -> extensión colineal
        # de la posición ideal del aceptor de ligando
        ('ASN', ('CG', 'ND2', 'HD21'), (HBOND_DIST, 180, 0), 'don'),
        ('ASN', ('CG', 'ND2', 'HD22'), (HBOND_DIST, 180, 0), 'don'),
        ('GLN', ('CD', 'NE2', 'HE21'), (HBOND_DIST, 180, 0), 'don'),
        ('GLN', ('CD', 'NE2', 'HE22'), (HBOND_DIST, 180, 0), 'don'),
        ('ARG', ('CZ', 'NE', 'HNE'), (HBOND_DIST, 180, 0), 'don'),
        ('ARG', ('CZ', 'NH1', 'HH11'), (HBOND_DIST, 180, 0), 'don'),
        ('ARG', ('CZ', 'NH1', 'HH12'), (HBOND_DIST, 180, 0), 'don'),
        ('ARG', ('CZ', 'NH2', 'HH21'), (HBOND_DIST, 180, 0), 'don'),
        ('ARG', ('CZ', 'NH2', 'HH22'), (HBOND_DIST, 180, 0), 'don'),
        ('ARG', ('CZ', 'NE', 'HE'), (HBOND_DIST, 180, 0), 'don'),
        ('SER', ('CB', 'OG', 'HG'), (HBOND_DIST, 180, 0), 'don'),
        ('THR', ('CB', 'OG1', 'HG1'), (HBOND_DIST, 180, 0), 'don'),
        ('TYR', ('CZ', 'OH', 'HH'), (HBOND_DIST, 180, 0), 'don'),
        ('TRP', ('CD1', 'NE1', 'HE1'), (HBOND_DIST, 180, 0), 'don'),
        ('LYS', ('CE', 'NZ', 'HZ1'), (HBOND_DIST, 180, 0), 'don'),
        ('LYS', ('CE', 'NZ', 'HZ2'), (HBOND_DIST, 180, 0), 'don'),
        ('LYS', ('CE', 'NZ', 'HZ3'), (HBOND_DIST, 180, 0), 'don'),
        ('HIE', ('CE1', 'NE2', 'HE2'), (HBOND_DIST, 180, 0), 'don'),
        ('HID', ('CE1', 'ND1', 'HD1'), (HBOND_DIST, 180, 0), 'don'),
        ('HIP', ('CE1', 'NE2', 'HE2'), (HBOND_DIST, 180, 0), 'don'),
        ('HIP', ('CE1', 'ND1', 'HD1'), (HBOND_DIST, 180, 0), 'don'),
    ]


def _backbone_interactions():
    """Igual que _sidechain_interactions pero para el esqueleto (aplica a
    cualquier residuo, ideal_interactions_bb()): carbonilo -> abanico de 5
    posiciones donor; NH -> 1 posición aceptor (colineal, no hay abanico
    porque no hay lone pairs múltiples que orientar)."""
    return [
        ('CA', 'C', 'O', HBOND_DIST, 180, 0, 'acc'),
        ('CA', 'C', 'O', HBOND_DIST, 210, 0, 'acc'),
        ('CA', 'C', 'O', HBOND_DIST, 240, 0, 'acc'),
        ('CA', 'C', 'O', HBOND_DIST, 150, 0, 'acc'),
        ('CA', 'C', 'O', HBOND_DIST, 120, 0, 'acc'),
        ('CA', 'N', 'H', HBOND_DIST, 180, 0, 'don'),
    ]


def _aromatic_defs():
    """[(resname, atomos_del_anillo, distancia_normal, desplazamiento_paralelo)]"""
    ring6 = ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')
    ring_his = ('CG', 'ND1', 'CD2', 'CE1', 'NE2')
    return [
        ('PHE', ring6, AROM_DIST, _PARALLEL_DISP),
        ('TYR', ring6, AROM_DIST, _PARALLEL_DISP),
        ('TRP', ('CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'), AROM_DIST, _PARALLEL_DISP),
        ('HIS', ring_his, AROM_DIST, _PARALLEL_DISP),
        ('HIE', ring_his, AROM_DIST, _PARALLEL_DISP),
        ('HID', ring_his, AROM_DIST, _PARALLEL_DISP),
        ('HIP', ring_his, AROM_DIST, _PARALLEL_DISP),
    ]


def _normal_vector(a, b, c, magnitude):
    n = np.cross(a - c, b - c)
    return n / np.linalg.norm(n) * magnitude


def _rotation_matrix(axis, theta):
    axis = np.asarray(axis, dtype=float)
    axis = axis / mt.sqrt(np.dot(axis, axis))
    a = mt.cos(theta / 2.0)
    b, c, d = -axis * mt.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
    ])


def _rotate_around_axis(axis, theta, v):
    return np.dot(_rotation_matrix(axis, theta), v)


def _dihedral(p1, p2, p3, p4):
    v1 = -1.0 * (p2 - p1)
    v2 = p3 - p2
    v3 = p4 - p3
    v2 = v2 / np.linalg.norm(v2)
    v = v1 - np.dot(v1, v2) * v2
    w = v3 - np.dot(v3, v2) * v2
    x = np.dot(v, w)
    y = np.dot(np.cross(v2, v), w)
    return np.degrees(np.arctan2(y, x))


def _pos4(p1, p2, p3, r4, a4, d4):
    """Coordenadas del punto 4 tal que dist(3,4)=r4, ángulo(2,3,4)=a4,
    diedro(1,2,3,4)=d4. Puerto directo de pos4() en ideal_interaction_sites.py."""
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    x3, y3, z3 = p3

    xejx = (y3 - y2) * (z1 - z2) - (z3 - z2) * (y1 - y2)
    yejx = -1 * ((x3 - x2) * (z1 - z2) - (z3 - z2) * (x1 - x2))
    zejx = (x3 - x2) * (y1 - y2) - (y3 - y2) * (x1 - x2)
    rejx = mt.sqrt(xejx ** 2 + yejx ** 2 + zejx ** 2)
    l1, m1, n1 = xejx / rejx, yejx / rejx, zejx / rejx

    r23 = mt.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2 + (z3 - z2) ** 2)

    xejz = yejx * (z3 - z2) - zejx * (y3 - y2)
    yejz = -1 * (xejx * (z3 - z2) - zejx * (x3 - x2))
    zejz = xejx * (y3 - y2) - yejx * (x3 - x2)
    rejz = mt.sqrt(xejz ** 2 + yejz ** 2 + zejz ** 2)
    l3, m3, n3 = xejz / rejz, yejz / rejz, zejz / rejz

    l2, m2, n2 = (x3 - x2) / r23, (y3 - y2) / r23, (z3 - z2) / r23

    d4r = d4 * mt.pi / 180
    a4r = (180 - a4) * mt.pi / 180

    z = r4 * mt.sin(a4r) * mt.cos(d4r)
    x = r4 * mt.sin(a4r) * mt.sin(d4r)
    y = r4 * mt.cos(a4r) + r23

    x4 = l1 * x + l2 * y + l3 * z + x2
    y4 = m1 * x + m2 * y + m3 * z + y2
    z4 = n1 * x + n2 * y + n3 * z + z2
    return np.array([x4, y4, z4])


def ideal_site_points(structure, chain_id, resids):
    """Para cada residuo en resids (números de posición) de structure[0][chain_id],
    calcula los puntos de bias "ideales" (ver docstring del módulo). Devuelve
    [(x, y, z, tipo)] con tipo en {'acc', 'don', 'aro'} -- mismo significado
    que en bias.py (tipo del átomo de LIGANDO que idealmente iría ahí), listo
    para export_bpf_points()/export_bpf_pdb_points()."""
    chain = structure[0][chain_id]
    points = []

    for resid in resids:
        try:
            residue = chain[(' ', int(resid), ' ')]
        except KeyError:
            continue
        resname = residue.get_resname()

        for atom in residue:
            name = atom.get_name()

            for r_name, (a1, a2, a3), (dist, ang, dih), role in _sidechain_interactions():
                if resname != r_name or name != a3:
                    continue
                try:
                    p1 = residue[a1].get_coord()
                    p2 = residue[a2].get_coord()
                    p3 = atom.get_coord()
                except KeyError:
                    continue

                site_type = _ROLE_TO_SITE_TYPE[role]
                if resname in ('SER', 'THR', 'TYR') and dih in ('alpha', 'beta'):
                    h_name = {'SER': 'HG', 'THR': 'HG1', 'TYR': 'HH'}[resname]
                    try:
                        p4 = residue[h_name].get_coord()
                    except KeyError:
                        continue
                    base_dihedral = _dihedral(p1, p2, p3, p4)
                    offset = 120 if dih == 'alpha' else 240
                    pos = _pos4(p1, p2, p3, dist, ang, base_dihedral + offset)
                else:
                    pos = _pos4(p1, p2, p3, dist, ang, dih)
                points.append((float(pos[0]), float(pos[1]), float(pos[2]), site_type))

            if resname not in ('HOH', 'WAT'):
                for a1, a2, a3, dist, ang, dih, role in _backbone_interactions():
                    if name != a3:
                        continue
                    try:
                        p1 = residue[a1].get_coord()
                        p2 = residue[a2].get_coord()
                        p3 = atom.get_coord()
                    except KeyError:
                        continue
                    pos = _pos4(p1, p2, p3, dist, ang, dih)
                    points.append((float(pos[0]), float(pos[1]), float(pos[2]),
                                   _ROLE_TO_SITE_TYPE[role]))

        for r_name, ring_atoms, arom_dist, par_disp in _aromatic_defs():
            if resname != r_name:
                continue
            try:
                coords = np.array([residue[a].get_coord() for a in ring_atoms], dtype=float)
            except KeyError:
                continue
            center = np.mean(coords, axis=0)
            normal = _normal_vector(coords[0], coords[1], coords[2], arom_dist)
            above, below = center + normal, center - normal
            for p in (above, below):
                points.append((float(p[0]), float(p[1]), float(p[2]), 'aro'))
            parallel = coords[0] - center
            par_mag = np.linalg.norm(parallel)
            for rot_ang in (0, 1.0472, 2.0944):  # 0°, 60°, 120° (por simetría +/- cubre 6 sentidos)
                rotated = _rotate_around_axis(normal, rot_ang, parallel) / par_mag * par_disp
                for p in (above, below):
                    for disp in (rotated, -rotated):
                        d = p + disp
                        points.append((float(d[0]), float(d[1]), float(d[2]), 'aro'))

    return points
