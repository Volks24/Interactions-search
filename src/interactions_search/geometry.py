"""Cálculos geométricos puros compartidos por el pipeline: centros de masa,
ángulos, centroides de anillos aromáticos, planaridad y volumen de envolvente
convexa. Sin dependencias de RDKit ni de las estructuras de datos del pipeline
(DF_Lig, DF_Active_Site, etc.) — solo coordenadas y vectores.
"""
from __future__ import annotations

import numpy as np
from Bio.PDB import Entity
from scipy.spatial import ConvexHull, QhullError

__all__ = [
    "center_of_mass",
    "angle_three_points",
    "center_aromatic_ring",
    "get_aromatic_coord",
    "convex_hull_volume",
    "_ring_planarity_rmsd",
]


def center_of_mass(entity, geometric=False):
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    """

    # Structure, Model, Chain, Residue
    if isinstance(entity, Entity.Entity):
        atom_list = entity.get_atoms()
    # List of Atoms
    elif hasattr(entity, '__iter__') and [x for x in entity if x.level == 'A']:
        atom_list = entity
    else: # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                            "Structure, Model, Chain, Residue, list of Atoms.")

    masses = []
    positions = [ [], [], [] ] # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]

    for atom in atom_list:
        masses.append(atom.mass)

        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError("Some Atoms don't have an element assigned.\n"
                         "Try adding them manually or calculate the geometrical center of mass instead.")

    if geometric:
        return [sum(coord_list)/len(masses) for coord_list in positions]
    else:
        w_pos = [ [], [], [] ]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index]*atom_mass)
            w_pos[1].append(positions[1][atom_index]*atom_mass)
            w_pos[2].append(positions[2][atom_index]*atom_mass)

        return [sum(coord_list)/sum(masses) for coord_list in w_pos]


def angle_three_points(Donor,Aceptor,Aceptor_Antecedent):

    Donor_coord = np.array(Donor)

    Aceptor_coord = np.array(Aceptor)

    Aceptor_Antecedent_coord = np.array(Aceptor_Antecedent)

    ba = Donor_coord-Aceptor_coord # normalization of vectors
    bc = Aceptor_Antecedent_coord-Aceptor_coord # normalization of vectors

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return (np.degrees(angle))  # calculated angle in radians to degree


def center_aromatic_ring(Aromatic_Ring):
    x,y,z = [],[],[]

    for j in range(0,len(Aromatic_Ring)):
        x.append(float(Aromatic_Ring[j][0]))
        y.append(float(Aromatic_Ring[j][1]))
        z.append(float(Aromatic_Ring[j][2]))

    CD1 = (x[1],y[1],z[1])
    CE1 = (x[3],y[3],z[3])

    vector_1 = (np.add(CD1, CE1))

    CD2 = (x[2],y[2],z[2])
    CE2 = (x[4],y[4],z[4])

    vector_2 = (np.add(CD2, CE2))

    center = (np.add(vector_2/2, vector_1/2))/2
    return(np.round(center,3))


def get_aromatic_coord(Res,AA):
    Aromatic_Ring = []
    if (Res == 'TYR') or (Res == 'PHE'):
            Coordenada = (AA.loc[AA['Atom'] == "CG", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD1", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE1", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
    elif ((Res) == 'TRP') :
            Coordenada = (AA.loc[AA['Atom'] == "CE3", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ3", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CH2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
    center = center_aromatic_ring(Aromatic_Ring)
    return(center[0],center[1],center[2])


def _ring_planarity_rmsd(pdb_coords, ring):
    """RMSD (Å) de los átomos del anillo respecto de su plano de mejor ajuste.
    Los anillos aromáticos son planos (RMSD ~0); un anillo saturado (ej. ciclohexano
    en silla) se desvía notablemente (~0.2-0.3 Å)."""
    coords = np.array([[pdb_coords[a][5], pdb_coords[a][6], pdb_coords[a][7]] for a in ring])
    centered = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centered)
    normal = vh[-1]
    return float(np.sqrt(np.mean((centered @ normal) ** 2)))


def convex_hull_volume(points):
    """Volumen (Å³) de la envolvente convexa de 'points' (array Nx3) vía
    scipy.spatial.ConvexHull. Retorna (nan, None) si hay menos de 4 puntos o si
    son coplanares/degenerados (QhullError) — no se puede definir un volumen 3D."""
    if len(points) < 4:
        return np.nan, None
    try:
        hull = ConvexHull(points)
        return round(float(hull.volume), 2), hull
    except QhullError:
        return np.nan, None
