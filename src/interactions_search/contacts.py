"""Búsqueda de contactos ligando-receptor por distancia (H-bond, hidrofóbico,
puente salino, π-catión) y validación por ángulo (H-bond, aromático)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem

__all__ = [
    "residuos_contacto",
    "Busqueda_Antecesor_Lig",
    "Interaccion_Aromatica",
    "aromatic_angle",
    "search_hydrophobic",
    "search_salt_bridges",
    "search_pi_cation",
]

_DF_COLS = ['Pos R', 'Res', 'Atom', 'Dist', 'Lig', 'Type', 'Angle', 'Interaction', 'LigID']


def residuos_contacto(Receptor_Caso,Lig_Caso,receptor_points,DF_Lig,DF_Interacciones,threshold_PH):

    Sub_Set_Receptor = receptor_points.query('Type == @Receptor_Caso')
    Matriz_receptor = np.array(Sub_Set_Receptor.iloc[:, [4, 5, 6]]).astype(float)

    Sub_Set_Ligando = DF_Lig.query('Caso == @Lig_Caso')

    for j in range(Sub_Set_Ligando.shape[0]):
        Coor_Lig = np.array(Sub_Set_Ligando.iloc[j, [1, 2, 3]]).astype(float)
        distances = np.linalg.norm(Matriz_receptor - Coor_Lig, axis=1)

        # Filtrar las distancias que son menores a 4.5
        within_distance_indices = np.where(distances < threshold_PH)[0]

        for idx in within_distance_indices:
            closest_data = Sub_Set_Receptor.iloc[idx]
            min_distance = distances[idx]

            # Agregar la información al DataFrame
            DF_Interacciones.loc[len(DF_Interacciones.index)] = [
                closest_data.iloc[1],  # X
                closest_data.iloc[2],  # Y
                closest_data.iloc[3],  # Z
                min_distance,           # Distancia
                Sub_Set_Ligando.iloc[j, 0],  # Nombre del átomo del ligando (solo display)
                Lig_Caso,               # Caso del ligando
                0.0, 'No',              # Angle (placeholder, se completa después), Interaction
                Sub_Set_Ligando.iloc[j, 5],  # Atom ID (serial único, para joins internos)
            ]
    return(DF_Interacciones)


def Busqueda_Antecesor_Lig(Atomo_ID,Lig_DF):
    # Atomo_ID es el serial único de átomo (columna 'Atom ID'), no el nombre:
    # el nombre puede repetirse en el ligando y matchear el átomo equivocado.
    punto_dado = np.array(Lig_DF.query('`Atom ID` == @Atomo_ID')[['X' , 'Y' , 'Z']])
    # Filtrar átomos que no sean H, excluyendo el propio átomo por serial exacto
    df_filtrado = Lig_DF[~Lig_DF['Element'].str.contains('H')]
    df_filtrado = df_filtrado[df_filtrado['Atom ID'] != Atomo_ID]
    # Calcular la distancia euclidiana
    df_filtrado['Distancia'] = np.sqrt((df_filtrado['X'] - punto_dado[0][0])**2 +
                                    (df_filtrado['Y'] - punto_dado[0][1])**2 +
                                    (df_filtrado['Z'] - punto_dado[0][2])**2)

    # Encontrar el índice del mínimo valor de distancia
    indice_min = df_filtrado['Distancia'].idxmin()

    # Obtener la fila con la distancia mínima
    coord = np.array(df_filtrado.loc[indice_min][['X' , 'Y' , 'Z']])


    return(coord)

def Interaccion_Aromatica(Anillo_Proteina,Anillo_Lig):
    # Anillo receptor #
    if (Anillo_Proteina.iloc[0,2] == 'TYR') or (Anillo_Proteina.iloc[0,2] == 'PHE'):
        Puntos_Interes = ['CG' , 'CD1' , 'CD2']
        Anillo_Name = 'CG-CD1-CD2'
        anillo_recept = np.array(Anillo_Proteina[Anillo_Proteina['Atom'].isin(Puntos_Interes)][['X','Y','Z']]).astype(float)
    elif (Anillo_Proteina.iloc[0,2]) == 'TRP':
        Puntos_Interes = ['CZ3' , 'CE3' , 'CH2']
        anillo_recept = np.array(Anillo_Proteina[Anillo_Proteina['Atom'].isin(Puntos_Interes)][['X','Y','Z']]).astype(float)
        Anillo_Name = 'CZ3-CE3-CH2'
    Anillo_Lig = (np.array(Anillo_Lig.iloc[0:3,1:4]))
    return(aromatic_angle(Anillo_Lig,anillo_recept))

def aromatic_angle(anillo_ligand, anillo_recept):
    # Encontrar los átomos comunes más cercanos
    atomos_comunes = [anillo_ligand[0], anillo_recept[1]]
    # Calcular los vectores normales a los planos aromáticos
    vector_normal1 = np.cross(anillo_ligand[1] - atomos_comunes[0], anillo_ligand[2] - atomos_comunes[0])
    vector_normal2 = np.cross(anillo_recept[2] - atomos_comunes[1], anillo_recept[0] - atomos_comunes[1])
    # Calcular el ángulo entre los vectores normales
    producto_punto = np.dot(vector_normal1, vector_normal2)
    norma_vector1 = np.linalg.norm(vector_normal1)
    norma_vector2 = np.linalg.norm(vector_normal2)
    # Calcular el ángulo en radianes y convertir a grados
    angulo_rad = np.arccos(producto_punto / (norma_vector1 * norma_vector2))
    angulo_deg = np.degrees(angulo_rad)
    # Asegurarse de que el ángulo esté en el rango de 0° a 90°
    if angulo_deg > 90:
        angulo_deg = 180 - angulo_deg

    return angulo_deg


# ──────────────────────────────────────────────────────────────────────────────
# Interacciones hidrofóbicas
# ──────────────────────────────────────────────────────────────────────────────

_HYDROPHOBIC_ATOMS = {
    'ALA': {'CB'},
    'VAL': {'CB', 'CG1', 'CG2'},
    'ILE': {'CB', 'CG1', 'CG2', 'CD1'},
    'LEU': {'CB', 'CG', 'CD1', 'CD2'},
    'MET': {'CB', 'CG', 'CE'},
    'PHE': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
    'TRP': {'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
    'PRO': {'CB', 'CG', 'CD'},
    'TYR': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2'},
}
_HPHO_LIG_SMARTS = '[c,C;!$([C,c]~[#7,#8,#16,#15,#9,#17,#35,#53])]'


def _collapse_same_residue_contacts(df):
    """Colapsa en una sola fila los contactos de un mismo átomo de ligando con
    varios átomos de un mismo residuo del receptor (ej: C20 contacta CG, CD1 y
    CD2 de una misma LEU): 1 contacto por residuo, distancia = promedio de todas
    las distancias átomo-átomo, 'Atom' lista los átomos involucrados."""
    if df.empty:
        return df
    # Agrupa por LigID (serial único del átomo del ligando), no por nombre: dos
    # átomos distintos pueden compartir nombre en ligandos mal nombrados.
    collapsed = df.groupby(['Pos R', 'LigID'], as_index=False).agg(
        Res=('Res', 'first'),
        Atom=('Atom', lambda s: ','.join(sorted(s.unique()))),
        Dist=('Dist', 'mean'),
        Lig=('Lig', 'first'),
        Type=('Type', 'first'),
        Angle=('Angle', 'first'),
        Interaction=('Interaction', 'first'),
    )
    collapsed['Dist'] = collapsed['Dist'].round(3)
    return collapsed[_DF_COLS]


def search_hydrophobic(mol, pdb_coords, DF_Active_Site, Distancia_Hidrofobica):
    """Contactos hidrofóbicos C-C entre ligando y residuos apolares del receptor."""
    pattern     = Chem.MolFromSmarts(_HPHO_LIG_SMARTS)
    hpho_idx    = {i for match in mol.GetSubstructMatches(pattern) for i in match}
    if not hpho_idx:
        return pd.DataFrame(columns=_DF_COLS)

    lig_pts    = [(pdb_coords[i][1], pdb_coords[i][5], pdb_coords[i][6], pdb_coords[i][7], pdb_coords[i][0])
                  for i in hpho_idx if i < len(pdb_coords)]
    lig_coords = np.array([[p[1], p[2], p[3]] for p in lig_pts], dtype=float)

    rec_mask   = DF_Active_Site.apply(
        lambda r: r['Atom'] in _HYDROPHOBIC_ATOMS.get(r['Residue'], set()), axis=1)
    rec_rows   = DF_Active_Site[rec_mask]
    if rec_rows.empty:
        return pd.DataFrame(columns=_DF_COLS)

    rec_coords = np.array(rec_rows[['X', 'Y', 'Z']], dtype=float)
    results = []
    for j, lig_pt in enumerate(lig_pts):
        dists = np.linalg.norm(rec_coords - lig_coords[j], axis=1)
        for k in np.where(dists < Distancia_Hidrofobica)[0]:
            r = rec_rows.iloc[k]
            results.append([int(r['Pos']), r['Residue'], r['Atom'],
                             round(dists[k], 3), lig_pt[0], 'hydrophobic', 0.0, 'Yes', lig_pt[4]])
    if not results:
        return pd.DataFrame(columns=_DF_COLS)
    return _collapse_same_residue_contacts(pd.DataFrame(results, columns=_DF_COLS))


# ──────────────────────────────────────────────────────────────────────────────
# Puentes salinos
# ──────────────────────────────────────────────────────────────────────────────

_SALT_POS_ATOMS = {'ARG': {'NH1', 'NH2', 'NE'}, 'LYS': {'NZ'}, 'HIP': {'ND1', 'NE2'}}
_SALT_NEG_ATOMS = {'ASP': {'OD1', 'OD2'}, 'GLU': {'OE1', 'OE2'}}
_SALT_DIST      = 4.0

_CATION_LIG_SMARTS = ['[N+;H3]', '[N+;H2]', '[N+;H1]', '[n+]', '[NH2]C(=[NH])[NH2]']
_ANION_LIG_SMARTS  = ['[O-]', '[$(C(=O)[OH])]', '[$(S(=O)(=O)[OH])]']


def search_salt_bridges(mol, pdb_coords, DF_Active_Site):
    """Detecta puentes salinos entre grupos cargados del ligando y del receptor."""
    def _lig_pts(smarts_list):
        idx = set()
        for s in smarts_list:
            pat = Chem.MolFromSmarts(s)
            if pat:
                for match in mol.GetSubstructMatches(pat):
                    idx.update(match)
        # Se conserva el serial (pdb_coords[i][0]) además del nombre: el nombre
        # puede repetirse en el ligando, y sin el serial no hay forma de volver
        # a ubicar la coordenada exacta del átomo (ej. para add_interaction_coords).
        return [(pdb_coords[i][0], pdb_coords[i][1], pdb_coords[i][5], pdb_coords[i][6], pdb_coords[i][7])
                for i in idx if i < len(pdb_coords)]

    cation_lig = _lig_pts(_CATION_LIG_SMARTS)
    anion_lig  = _lig_pts(_ANION_LIG_SMARTS)
    results = []

    for rec in DF_Active_Site.itertuples(index=False):
        rc = np.array([rec.X, rec.Y, rec.Z])
        if rec.Atom in _SALT_NEG_ATOMS.get(rec.Residue, set()):
            for serial, atom, x, y, z in cation_lig:
                d = np.linalg.norm(rc - np.array([x, y, z]))
                if d < _SALT_DIST:
                    results.append([rec.Pos, rec.Residue, rec.Atom, round(d,3),
                                     atom, 'salt_bridge', 0.0, 'Yes', serial])
        if rec.Atom in _SALT_POS_ATOMS.get(rec.Residue, set()):
            for serial, atom, x, y, z in anion_lig:
                d = np.linalg.norm(rc - np.array([x, y, z]))
                if d < _SALT_DIST:
                    results.append([rec.Pos, rec.Residue, rec.Atom, round(d,3),
                                     atom, 'salt_bridge', 0.0, 'Yes', serial])
    return pd.DataFrame(results, columns=_DF_COLS) if results else pd.DataFrame(columns=_DF_COLS)


# ──────────────────────────────────────────────────────────────────────────────
# Interacciones π-catión
# ──────────────────────────────────────────────────────────────────────────────

_PI_CATION_DIST     = 5.0
_CATION_REC_ATOMS   = {'ARG': {'NH1', 'NH2', 'NE'}, 'LYS': {'NZ'},
                        'HIS': {'ND1', 'NE2'}, 'HIP': {'ND1', 'NE2'}}


def search_pi_cation(DF_Active_Site, aromatic_lig_df):
    """Detecta interacciones π-catión: anillo aromático del ligando vs catión del receptor."""
    results = []
    for cas in aromatic_lig_df['Caso'].unique():
        ring_atoms  = aromatic_lig_df[aromatic_lig_df['Caso'] == cas]
        ring_center = np.mean(np.array(ring_atoms[['Coord X', 'Coord Y', 'Coord Z']]).astype(float), axis=0)
        for rec in DF_Active_Site.itertuples(index=False):
            if rec.Atom in _CATION_REC_ATOMS.get(rec.Residue, set()):
                d = np.linalg.norm(ring_center - np.array([rec.X, rec.Y, rec.Z]))
                if d < _PI_CATION_DIST:
                    results.append([rec.Pos, rec.Residue, rec.Atom, round(d,3),
                                     cas, 'pi_cation', 0.0, 'Yes', np.nan])
    return pd.DataFrame(results, columns=_DF_COLS) if results else pd.DataFrame(columns=_DF_COLS)
