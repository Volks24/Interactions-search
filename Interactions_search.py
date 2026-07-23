### Librerias ###

import sys
import yaml  # noqa: F401 — kept for backward compat; carga_variables now delegates to load_config
from pathlib import Path
from interactions_search.config import load_config
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import pandas as pd
import numpy as np
from Bio.PDB import *
import math
import argparse
import os
import shutil
from scipy.spatial import ConvexHull, QhullError
import matplotlib
matplotlib.use('Agg')  # headless: sin esto matplotlib puede requerir un $DISPLAY inexistente en servidores
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# Función para obtener las coordenadas 2D de los átomos
def get_atom_coords(mol, atom_idx):
    conf = mol.GetConformer()
    pos = conf.GetAtomPosition(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    return f"{atom.GetSymbol()} {atom_idx}: ({pos.x}, {pos.y}, {pos.z})"

# Leer el archivo PDB y extraer la información
def extract_coords_from_pdb(pdb_filename):
    coords = []
    CM = []
    with open(pdb_filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21]
                res_seq = line[22:26].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom_id = int(line[6:11].strip())
                CM.append([x,y,z])
                coords.append((atom_id, atom_name, res_name, chain_id, res_seq, x, y, z))
        CM_Coord_Set = np.array(CM)
        center_of_mass = np.mean(CM_Coord_Set, axis=0)

    return (coords,center_of_mass)


def _draw_mol_labeled(mol, highlight_atoms, atom_labels, filename, size=(600, 600)):
    """Dibuja la molécula con átomos resaltados y etiquetas atomNote."""
    mol_copy = Chem.RWMol(Chem.Mol(mol))
    rdDepictor.Compute2DCoords(mol_copy)
    for idx, lbl in atom_labels.items():
        mol_copy.GetAtomWithIdx(idx).SetProp('atomNote', str(lbl))
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    drawer.drawOptions().addAtomIndices = False
    colors = {idx: (0.9, 0.35, 0.35) for idx in highlight_atoms}
    radii  = {idx: 0.4 for idx in highlight_atoms}
    drawer.DrawMolecule(mol_copy,
                        highlightAtoms=list(highlight_atoms),
                        highlightAtomColors=colors,
                        highlightAtomRadii=radii)
    drawer.FinishDrawing()
    with open(filename, 'wb') as fh:
        fh.write(drawer.GetDrawingText())


def search_hot_points(Ligand_imput, mol, pdb_coords, ligand_plot, folder):

    acceptor_smarts = ['[O;H1]', '[O;H0]', '[N;H1]', '[N;H0]', '[n]', '[o]', '[N+]']
    donor_smarts    = ['[O;H]', '[N;H2]', '[N;H]', '[S;H]', '[nH]']

    acceptor_atoms, donor_atoms = [], []
    for smarts in acceptor_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        for match in mol.GetSubstructMatches(pattern):
            for atom_idx in match:
                acceptor_atoms.append(atom_idx)

    for smarts in donor_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        for match in mol.GetSubstructMatches(pattern):
            for atom_idx in match:
                donor_atoms.append(atom_idx)

    if ligand_plot == 'Yes':
        stem = Path(Ligand_imput).stem
        acc_labels = {idx: pdb_coords[idx][1] for idx in acceptor_atoms if idx < len(pdb_coords)}
        _draw_mol_labeled(mol, acceptor_atoms, acc_labels, f"{folder}/{stem}_acceptors.png")
        don_labels = {idx: pdb_coords[idx][1] for idx in donor_atoms if idx < len(pdb_coords)}
        _draw_mol_labeled(mol, donor_atoms, don_labels, f"{folder}/{stem}_donors.png")

    return acceptor_atoms, donor_atoms


def active_site_residues(structure, Ligando_Centro,cadena, centroid_distance , lig):
    model = structure[0][cadena]

    active_site = pd.DataFrame(columns=['Serial', 'Pos', 'Residue', 'Atom', 'X' , 'Y' ,'Z', 'CM X' , 'CM Y' , 'CM Z' ])

    Residuos_Interes = []

    for residue in model.get_residues():
        Residuo_Center = list(center_of_mass(residue))
        if (math.dist(Ligando_Centro, Residuo_Center)) < centroid_distance:
            Residuos_Interes.append([residue.get_resname(), residue.get_id()[1]])
            for atom in residue:
                Res_name = residue.get_resname()
                Res_id = residue.get_id()[1]
                atom_name = atom.get_name()
                Coor = list(atom.get_coord())
                Serial = atom.get_serial_number()
                #atoms.append([Serial, Res_id, Res_name, atom_name, Coor, Residuo_Center])
                active_site.loc[len(active_site.index)] = [Serial ,Res_id, Res_name, atom_name ,round(float(Coor[0]),3),round(float(Coor[1]),3),round(float(Coor[2]),3),round(Residuo_Center[0],3),round(Residuo_Center[1],3),round(Residuo_Center[2],3)]
    
    ### Elimino cosas que no sirven ###
    active_site = active_site[active_site['Residue'] != 'HOH']
    active_site = active_site[active_site['Residue'] != lig]
    return active_site

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



def Coordenadas_interes_receptor(Aceptores_Prot,Dadores_Prot,DF_Active_Site):
    ### Obtengo las coordenadas de los atomos de interes en el receptor
    receptor_points = pd.DataFrame(columns=['Type','Pos','Residue', 'Atom', 'X' , 'Y' , 'Z'])
    for pos in range(0,DF_Active_Site.shape[0]):
        Atomo = (DF_Active_Site.iloc[pos,2])
        Res = (DF_Active_Site.iloc[pos,3])
        listado = Aceptores_Prot.get(Atomo, [])
        if Res in listado:
            receptor_points.loc[len(receptor_points.index)] = 'Aceptor',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],DF_Active_Site.iloc[pos,3],DF_Active_Site.iloc[pos,4],DF_Active_Site.iloc[pos,5],DF_Active_Site.iloc[pos,6]
    for pos in range(0,DF_Active_Site.shape[0]):
        Atomo = (DF_Active_Site.iloc[pos,2])
        Res = (DF_Active_Site.iloc[pos,3])
        listado = Dadores_Prot.get(Atomo, [])
        if Res in listado:
            receptor_points.loc[len(receptor_points.index)] = 'Dador',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],DF_Active_Site.iloc[pos,3],DF_Active_Site.iloc[pos,4],DF_Active_Site.iloc[pos,5],DF_Active_Site.iloc[pos,6]
    aa_aro = ['TYR' , 'PHE' , 'TRP']
    for pos in range(0,DF_Active_Site.shape[0]):
        Atomo = (DF_Active_Site.iloc[pos,2])
        ID = (DF_Active_Site.iloc[pos,1])
        if Atomo in aa_aro:
            Sub_Set = DF_Active_Site.query('Pos == @ID')
            x,y,z = get_aromatic_coord(Atomo,Sub_Set)
            if ID not in receptor_points['Pos'].values:
                receptor_points.loc[len(receptor_points.index)] = 'aromatic',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],'center',x,y,z
            elif 'aromatic' not in (receptor_points.query('Pos == @ID')['Type'].tolist()) :# Solo posicion 
                receptor_points.loc[len(receptor_points.index)] = 'aromatic',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],'center',x,y,z
        
            
    return(receptor_points)



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

def generate_df_ligand(pdb_coords):

    columns = ['Atom ID', 'Element', 'Residue Name', 'Chain ID', 'Residue Number', 'X', 'Y', 'Z']
    df_ligand = pd.DataFrame(pdb_coords, columns=columns)

    return(df_ligand)


# ──────────────────────────────────────────────────────────────────────────────
# Bias probe file (.bpf) para GOLD
# ──────────────────────────────────────────────────────────────────────────────

# Vset/r fijos por tipo, tomados de receptor_fs(6).bpf (referencia de formato
# provista) — no varían por átomo, son la misma convención para todo don/acc/aro.
_BPF_PARAMS = {
    'don': {'Vset': -2.72, 'r': 1.20},
    'acc': {'Vset': -2.28, 'r': 0.80},
    'aro': {'Vset': -2.00, 'r': 2.00},
}


def _collect_bias_points(DF_Lig, DF_true=None):
    """Junta los puntos de bias (un punto por átomo aceptor/donor, un punto
    por anillo aromático en su centroide) en el mismo orden que usan
    export_bpf() y export_bpf_pdb(), para que ambos archivos queden
    alineados fila a fila.

    Si DF_true es None (default): usa TODOS los hot-points químicos del
    ligando (search_hot_points/search_rings), sin importar si esa posición
    llegó a formar una interacción real con el receptor en esta pose.

    Si se pasa DF_true (interacciones validadas, Interaction == 'Yes'): se
    restringe a los hot-points que sí participan en una interacción validada
    — átomos aceptor/donor por serial (LigID) y anillos por su etiqueta
    ('aromatic' o 'pi_cation', ambos indican que el anillo hace contacto
    real). 'Evidence-based' en vez de 'todo lo químicamente posible'."""
    donors    = DF_Lig[DF_Lig['Caso'] == 'donor']
    acceptors = DF_Lig[DF_Lig['Caso'] == 'acceptor']
    aromatic  = DF_Lig[DF_Lig['Caso'].astype(str).str.startswith('aromatic')]

    if DF_true is not None:
        validated_ligids = set(
            DF_true.loc[DF_true['Type'].isin(['acceptor', 'donor']), 'LigID'].dropna().astype(int))
        donors    = donors[donors['Atom ID'].astype(int).isin(validated_ligids)]
        acceptors = acceptors[acceptors['Atom ID'].astype(int).isin(validated_ligids)]
        validated_rings = set(DF_true.loc[DF_true['Type'].isin(['aromatic', 'pi_cation']), 'Lig'])
        aromatic  = aromatic[aromatic['Caso'].isin(validated_rings)]

    rows = []
    for _, r in donors.iterrows():
        rows.append((float(r['Coord X']), float(r['Coord Y']), float(r['Coord Z']), 'don'))
    for _, r in acceptors.iterrows():
        rows.append((float(r['Coord X']), float(r['Coord Y']), float(r['Coord Z']), 'acc'))
    for _, group in aromatic.groupby('Caso'):
        center = group[['Coord X', 'Coord Y', 'Coord Z']].astype(float).mean(axis=0)
        rows.append((float(center['Coord X']), float(center['Coord Y']), float(center['Coord Z']), 'aro'))
    return rows


def export_bpf(DF_Lig, filepath, DF_true=None):
    """Genera un archivo .bpf (bias probe file, formato GOLD: header 'x y z
    Vset r type') a partir de los hot-points del ligando (ver
    _collect_bias_points; DF_true filtra a solo los validados). Vset y r son
    fijos por tipo (_BPF_PARAMS)."""
    rows = _collect_bias_points(DF_Lig, DF_true)
    with open(filepath, 'w') as f:
        f.write('x\ty\tz\tVset\tr\ttype\n')
        for x, y, z, tipo in rows:
            p = _BPF_PARAMS[tipo]
            f.write(f"{x:.3f}\t{y:.3f}\t{z:.3f}\t{p['Vset']}\t{p['r']}\t{tipo}\n")


_BPF_RESNAME = {'don': 'DON', 'acc': 'ACC', 'aro': 'ARO'}


def export_bpf_pdb(DF_Lig, filepath, DF_true=None):
    """PDB 'dummy' con un átomo H por punto de bias (mismos puntos y mismo
    orden que export_bpf; DF_true filtra a solo los validados), para poder
    cargar y visualizar los puntos de bias en VMD junto al receptor/ligando.
    resname = DON/ACC/ARO según el tipo, chain 'X', un residuo por átomo
    (dummy, sin significado bioquímico)."""
    rows = _collect_bias_points(DF_Lig, DF_true)
    with open(filepath, 'w') as f:
        for i, (x, y, z, tipo) in enumerate(rows, start=1):
            resname = _BPF_RESNAME[tipo]
            f.write(
                f"ATOM  {i:>5} {'H':<4} {resname:>3} X{i:>4}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {'H':>2}\n"
            )
        f.write('END\n')



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

def angle_three_points(Donor,Aceptor,Aceptor_Antecedent):
    
    Donor_coord = np.array(Donor)

    Aceptor_coord = np.array(Aceptor)

    Aceptor_Antecedent_coord = np.array(Aceptor_Antecedent)

    ba = Donor_coord-Aceptor_coord # normalization of vectors
    bc = Aceptor_Antecedent_coord-Aceptor_coord # normalization of vectors

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return (np.degrees(angle))  # calculated angle in radians to degree 

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

def _ring_planarity_rmsd(pdb_coords, ring):
    """RMSD (Å) de los átomos del anillo respecto de su plano de mejor ajuste.
    Los anillos aromáticos son planos (RMSD ~0); un anillo saturado (ej. ciclohexano
    en silla) se desvía notablemente (~0.2-0.3 Å)."""
    coords = np.array([[pdb_coords[a][5], pdb_coords[a][6], pdb_coords[a][7]] for a in ring])
    centered = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centered)
    normal = vh[-1]
    return float(np.sqrt(np.mean((centered @ normal) ** 2)))


def search_rings(mol, pdb_coords, numero_anillo_aromatico, planarity_rmsd_max):
    """Identifica anillos aromáticos: tamaño > numero_anillo_aromatico y planos
    (RMSD respecto del plano de mejor ajuste por debajo de planarity_rmsd_max).
    No se usa mol.GetIsAromatic(): RDKit no perfila aromaticidad de forma fiable
    para moléculas leídas desde PDB (sin órdenes de enlace explícitos), así que
    la planaridad geométrica 3D es el criterio real de aromaticidad aquí."""
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()
    ring_data = []
    for ring in ring_atoms:
        if len(ring) > numero_anillo_aromatico and \
           _ring_planarity_rmsd(pdb_coords, ring) <= planarity_rmsd_max:
            ring_data.append({'Ring': len(ring_data) + 1, 'Atoms': ring, 'Ring Size': len(ring)})
    rings_data = []
    for ring in ring_data:
        label = f"aromatic {ring['Ring']} (#{ring['Ring Size']})"
        for atom in ring['Atoms']:
            rings_data.append([pdb_coords[atom][1], pdb_coords[atom][5],
                                pdb_coords[atom][6], pdb_coords[atom][7], label])
    return ring_data, rings_data

_RING_COLORS = [
    (0.9, 0.35, 0.35),   # R1 — rojo
    (0.25, 0.55, 0.9),   # R2 — azul
    (0.2,  0.78, 0.45),  # R3 — verde
    (0.95, 0.70, 0.15),  # R4 — amarillo
]

def visualize_rings(mol, ring_data, Ligand_imput, folder):
    mol_copy  = Chem.RWMol(Chem.Mol(mol))
    rdDepictor.Compute2DCoords(mol_copy)
    highlight, colors, radii = [], {}, {}
    for ring in ring_data:
        if ring['Ring Size'] > 5:
            rnum  = ring['Ring']
            color = _RING_COLORS[(rnum - 1) % len(_RING_COLORS)]
            atoms = ring['Atoms']
            for atom in atoms:
                highlight.append(atom)
                colors[atom] = color
                radii[atom]  = 0.4
            # Etiqueta en el átomo central del anillo
            mid = atoms[len(atoms) // 2]
            mol_copy.GetAtomWithIdx(mid).SetProp('atomNote', f'R{rnum}')
    drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
    drawer.drawOptions().addAtomIndices = False
    drawer.DrawMolecule(mol_copy,
                        highlightAtoms=highlight,
                        highlightAtomColors=colors,
                        highlightAtomRadii=radii)
    drawer.FinishDrawing()
    with open(f"{folder}/{Path(Ligand_imput).stem}_aromatic.png", 'wb') as fh:
        fh.write(drawer.GetDrawingText())
    

def _vmd_write_interaction(VDM_TCL, j, chain, resid, resname, coord1, coord2, color,
                            mol_receptor='$molReceptor', mol_graphics='$molLigand'):
    """Escribe en el .tcl la línea punteada + etiqueta de distancia entre coord1 (ligando)
    y coord2 (receptor) para una interacción. Común a aromatic/acceptor/donor.

    mol_receptor/mol_graphics son variables Tcl (seteadas por el caller con [mol new ...])
    en vez de molids fijos: si la sesión de VMD ya tenía moléculas cargadas antes de
    correr este script, 'mol new' no asigna 0/1 sino los siguientes ids libres, y
    hardcodear 0/'top' hace fallar 'atomselect'/'graphics' con 'invalid molecule'."""
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    VDM_TCL.write(f'graphics {mol_graphics} color {color}\n')
    VDM_TCL.write(f'graphics {mol_graphics} line {{{x1} {y1} {z1}}} {{{x2} {y2} {z2}}} width 5 style dashed\n')
    VDM_TCL.write(f'set Recptor{j} [atomselect {mol_receptor} "chain {chain} and resid {resid} and resname {resname} and name CZ"]\n')
    VDM_TCL.write(f'set x1 {{{x1}}}\n')
    VDM_TCL.write(f'set y1 {{{y1}}}\n')
    VDM_TCL.write(f'set z1 {{{z1}}}\n')
    VDM_TCL.write(f'set x2 {{{x2}}}\n')
    VDM_TCL.write(f'set y2 {{{y2}}}\n')
    VDM_TCL.write(f'set z2 {{{z2}}}\n')
    VDM_TCL.write('set dx [expr {$x1 - $x2}]\n')
    VDM_TCL.write('set dy [expr {$y1 - $y2}]\n')
    VDM_TCL.write('set dz [expr {$z1 - $z2}]\n')
    VDM_TCL.write('set distance [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]\n')
    VDM_TCL.write('set xm [expr {($x1 + $x2) / 2}]\n')
    VDM_TCL.write('set ym [expr {($y1 + $y2) / 2}]\n')
    VDM_TCL.write('set zm [expr {($z1 + $z2) / 2}]\n')
    VDM_TCL.write(f'graphics {mol_graphics} color white\n')
    VDM_TCL.write(f'graphics {mol_graphics} text [list $xm $ym $zm] [format "%.2f A" $distance]\n')


_VMD_COLORS = {'aromatic': 'white', 'acceptor': 'red', 'donor': 'yellow'}


def _write_hbond_aromatic_lines(VDM_TCL, DF_Interacciones, receptor_points, DF_Lig, chain):
    """Escribe las líneas punteadas + etiquetas de distancia de cada interacción
    H-bond/aromática de DF_Interacciones. Compartido por scripting_vmd() y
    scripting_vmd_combined() para no duplicar la lógica de lookup de coordenadas."""
    # Acceso por nombre de columna (no por posición): DF_Interacciones puede
    # traer X/Y/Z insertadas en medio del esquema original (ver
    # add_interaction_coords), lo que corre las posiciones fijas de columna.
    for j in range(0, DF_Interacciones.shape[0]):
        row     = DF_Interacciones.iloc[j]
        tipo    = row['Type']
        Recept  = row['Pos R']
        Resname = row['Res']

        if tipo == 'aromatic':
            Coord2 = np.array(receptor_points[(receptor_points['Pos'] == Recept) &
                                               (receptor_points['Atom'] == 'center')][['X','Y','Z']])[0]
            anillo = row['Lig']
            punto_dado = np.array(DF_Lig.query('Caso == @anillo')[['Coord X' , 'Coord Y' , 'Coord Z']])
            Coord1 = np.mean(punto_dado, axis=0)
        elif tipo in ('acceptor', 'donor'):
            Coord2 = np.array(receptor_points[(receptor_points['Pos'] == Recept) &
                                               (receptor_points['Atom'] == row['Atom'])][['X','Y','Z']])[0]
            atomo_id = row['LigID']
            Coord1 = np.array(DF_Lig[(DF_Lig['Atom ID'] == atomo_id)][['Coord X' , 'Coord Y' , 'Coord Z']])[0]
        else:
            continue

        _vmd_write_interaction(VDM_TCL, j, chain, Recept, Resname, Coord1, Coord2, _VMD_COLORS[tipo])


def scripting_vmd(DF_Interacciones,receptor_points,aromatic_lig_df,DF_Lig,Prot,chain,Lig,folder):
    receptor_name = Path(Prot).stem
    Lig_name = Path(Lig).stem
    # Nombres de archivo (no la ruta original): analyze_pair ya copió ambos PDB a
    # 'folder' antes de llamar a este script, y el .tcl queda guardado ahí mismo.
    # Referenciar la ruta original rompe con --complex (usa un tmp_dir que se
    # borra al terminar) y además hace el .tcl no portable si se mueve la carpeta.
    Prot_file = Path(Prot).name
    Lig_file  = Path(Lig).name

    Res_All = DF_Interacciones['Pos R'].tolist()
    residues = ' '.join(map(str, Res_All))

    with open(f'{folder}/vmd_{receptor_name}_{Lig_name}.tcl', 'w') as VDM_TCL:
        # Cargar el archivo PDB. Se captura el molid real en variables Tcl en vez de
        # asumir 0/1: si la sesión de VMD ya tenía moléculas cargadas, 'mol new' no
        # asigna esos ids y hardcodearlos rompe atomselect/graphics más abajo.
        VDM_TCL.write(f'display projection orthographic\n')
        VDM_TCL.write(f'set molReceptor [mol new "{Prot_file}"]\n')
        VDM_TCL.write(f'mol modselect 0 $molReceptor all\n')
        # Lines en vez de Tube/NewCartoon/Trace: en builds alpha de VMD (ej. 2.0.0a9)
        # esas representaciones basadas en spline por backbone truncan la geometría
        # a ~40 residuos sin avisar (bug confirmado: seleccionar explícitamente un
        # tramo lejano no dibuja nada), sea cual sea la selección o la molécula.
        # Lines no depende de ese cálculo (dibuja enlace por enlace) y siempre
        # muestra la proteína completa.
        VDM_TCL.write(f'mol modstyle 0 $molReceptor Lines 3\n')
        VDM_TCL.write(f'mol modcolor 0 $molReceptor ColorID 6\n')
        VDM_TCL.write(f'mol modmaterial 0 $molReceptor Opaque\n')
        VDM_TCL.write(f'mol addrep $molReceptor\n')
        VDM_TCL.write(f'mol modselect 1 $molReceptor resid {residues} and chain {chain}\n')
        VDM_TCL.write(f'mol modstyle 1 $molReceptor Licorice\n')
        VDM_TCL.write(f'set molLigand [mol new "{Lig_file}"]\n')
        # Crear una representación en Licorice para el ligando
        VDM_TCL.write(f'mol addrep $molLigand\n')
        VDM_TCL.write(f'mol modstyle 0 $molLigand Licorice\n')
        # Sin esto la cámara queda encuadrada según la última molécula cargada (el
        # ligando, mucho más chico) y no se reajusta tras cargar el receptor antes:
        # recorta tramos enteros de la proteína por los planos de clipping. Al cargar
        # manualmente por GUI, VMD hace resetview solo; en modo scripteado (-e) no.
        VDM_TCL.write(f'display resetview\n')

        ### Busco Interaccion
        _write_hbond_aromatic_lines(VDM_TCL, DF_Interacciones, receptor_points, DF_Lig, chain)


_VMD_HYDROPHOBIC_COLOR = 'orange'


def scripting_vmd_hydrophobic(DF_Interacciones, DF_Active_Site, DF_Lig_All, Prot, chain, Lig, folder):
    """Genera un .tcl de VMD exclusivo para los contactos hidrofóbicos (línea naranja).
    'Atom' puede traer varios átomos del mismo residuo colapsados (ej. 'CD1,CD2,CG');
    el punto del receptor es el centroide de esos átomos."""
    DF_Hpho = DF_Interacciones[DF_Interacciones['Type'] == 'hydrophobic']
    if DF_Hpho.empty:
        return

    receptor_name = Path(Prot).stem
    Lig_name = Path(Lig).stem
    Prot_file = Path(Prot).name
    Lig_file  = Path(Lig).name

    residues = ' '.join(map(str, DF_Hpho['Pos R'].tolist()))

    with open(f'{folder}/vmd_hydrophobic_{receptor_name}_{Lig_name}.tcl', 'w') as VDM_TCL:
        VDM_TCL.write(f'display projection orthographic\n')
        VDM_TCL.write(f'set molReceptor [mol new "{Prot_file}"]\n')
        VDM_TCL.write(f'mol modselect 0 $molReceptor all\n')
        VDM_TCL.write(f'mol modstyle 0 $molReceptor Lines 3\n')
        VDM_TCL.write(f'mol modcolor 0 $molReceptor ColorID 6\n')
        VDM_TCL.write(f'mol modmaterial 0 $molReceptor Opaque\n')
        VDM_TCL.write(f'mol addrep $molReceptor\n')
        VDM_TCL.write(f'mol modselect 1 $molReceptor resid {residues} and chain {chain}\n')
        VDM_TCL.write(f'mol modstyle 1 $molReceptor Licorice\n')
        VDM_TCL.write(f'set molLigand [mol new "{Lig_file}"]\n')
        VDM_TCL.write(f'mol addrep $molLigand\n')
        VDM_TCL.write(f'mol modstyle 0 $molLigand Licorice\n')
        VDM_TCL.write(f'display resetview\n')

        for j, row in enumerate(DF_Hpho.itertuples(index=False)):
            rec_atoms = row.Atom.split(',')
            rec_sub = DF_Active_Site[(DF_Active_Site['Pos'] == row._0) &
                                      (DF_Active_Site['Atom'].isin(rec_atoms))]
            Coord2 = rec_sub[['X', 'Y', 'Z']].astype(float).mean(axis=0).values
            lig_sub = DF_Lig_All[DF_Lig_All['Atom ID'] == row.LigID]
            Coord1 = lig_sub[['X', 'Y', 'Z']].astype(float).values[0]
            _vmd_write_interaction(VDM_TCL, j, chain, row._0, row.Res, Coord1, Coord2,
                                   _VMD_HYDROPHOBIC_COLOR)


def remove_bias(file_path, folder):
    old_file_path = Path(file_path).stem + '_old.pdb'
    shutil.copy(file_path, f'{folder}/{old_file_path}')

    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Filtrar las líneas que no contienen "CM"
    lines = [line for line in lines if ' CM ' not in line]

    # Guardar el archivo sin las líneas "CM"
    with open(file_path, 'w') as f:
        f.writelines(lines)


def split_pdb(pdb_path, output_dir='.', exclude_water=True, force_ligand_names=None):
    """
    Separa un PDB complejo en proteína (ATOM) y un archivo por cada grupo HETATM único.
    El agua (HOH, WAT, TIP3, SOL) se excluye por defecto.
    Los registros CONECT se distribuyen al archivo del grupo HETATM correspondiente.

    force_ligand_names : set de resnames (ej: {'TF3', '7FW'}) que se tratan como ligando
                         aunque estén escritos como ATOM en vez de HETATM en el PDB.

    Retorna:
        protein_path : Path al PDB de la proteína
        het_paths    : dict {resname: Path}  —  vacío si no hay HETATM
    """
    from collections import defaultdict

    WATER_NAMES = {'HOH', 'WAT', 'TIP', 'TIP3', 'SOL', 'DOD'}
    force_ligand_names = {n.upper() for n in force_ligand_names} if force_ligand_names else set()

    protein_lines = []
    het_lines     = defaultdict(list)   # resname -> [líneas HETATM]
    conect_lines  = []
    header_lines  = []

    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == 'ATOM':
                resname = line[17:20].strip()
                if resname in force_ligand_names:
                    het_lines[resname].append('HETATM' + line[6:])
                else:
                    protein_lines.append(line)
            elif rec == 'HETATM':
                resname = line[17:20].strip()
                if exclude_water and resname in WATER_NAMES:
                    continue
                het_lines[resname].append(line)
            elif rec == 'CONECT':
                conect_lines.append(line)
            elif rec in ('TER', 'REMARK', 'HEADER', 'TITLE', 'COMPND', 'SOURCE', 'SEQRES'):
                protein_lines.append(line)
                header_lines.append(line)

    out  = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = Path(pdb_path).stem

    # --- Proteína ---
    protein_path = out / f'{stem}_protein.pdb'
    with open(protein_path, 'w') as f:
        f.writelines(protein_lines)
        if not protein_lines or not protein_lines[-1].startswith('END'):
            f.write('END\n')

    # --- Grupos HETATM ---
    # Pre-indexar seriales de cada grupo para filtrar CONECT
    group_serials = {}
    for resname, lines in het_lines.items():
        serials = set()
        for l in lines:
            try:
                serials.add(int(l[6:11]))
            except ValueError:
                pass
        group_serials[resname] = serials

    het_paths = {}
    for resname, lines in het_lines.items():
        het_path = out / f'{stem}_{resname}.pdb'
        with open(het_path, 'w') as f:
            f.writelines(lines)
            # CONECT cuyos átomos pertenecen a este grupo
            my_serials = group_serials[resname]
            for cl in conect_lines:
                referenced = set()
                for i in range(6, min(len(cl.rstrip()), 31), 5):
                    tok = cl[i:i+5].strip()
                    if tok:
                        try:
                            referenced.add(int(tok))
                        except ValueError:
                            pass
                if referenced & my_serials:
                    f.write(cl)
            f.write('END\n')
        het_paths[resname] = het_path

    return protein_path, het_paths



# ──────────────────────────────────────────────────────────────────────────────
# Validación de inputs
# ──────────────────────────────────────────────────────────────────────────────

def validate_inputs(receptor_pdb, ligand_pdb, chain):
    """Valida archivos y cadena antes del análisis. Retorna lista de errores."""
    errors = []
    cwd = Path.cwd()
    pdbs_in_cwd = sorted(p.name for p in cwd.glob('*.pdb'))

    if not Path(receptor_pdb).exists():
        errors.append(f"Receptor not found: {receptor_pdb}")
        errors.append(f"  Current directory : {cwd}")
        errors.append(f"  Available PDBs    : {pdbs_in_cwd or '(none)'}")
        return errors
    if not Path(ligand_pdb).exists():
        errors.append(f"Ligand not found: {ligand_pdb}")
        errors.append(f"  Current directory : {cwd}")
        errors.append(f"  Available PDBs    : {pdbs_in_cwd or '(none)'}")
        return errors
    chains_found, has_atoms = set(), False
    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith('ATOM'):
                chains_found.add(line[21])
                has_atoms = True
    if not has_atoms:
        errors.append(f"Receptor has no ATOM records: {receptor_pdb}")
    elif chain not in chains_found:
        errors.append(f"Chain '{chain}' not found. Available: {sorted(chains_found)}")
    lig_atoms = sum(1 for line in open(ligand_pdb) if line.startswith(('ATOM', 'HETATM')))
    if lig_atoms == 0:
        errors.append(f"Ligand has no atoms: {ligand_pdb}")
    return errors


# ──────────────────────────────────────────────────────────────────────────────
# Volumen (envolvente convexa) — compartido entre sitio activo y pockets
# ──────────────────────────────────────────────────────────────────────────────

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


def plot_hull_volume(points, hull, title, filename):
    """PNG de dispersión 3D: todos los puntos (negro) y los vértices de su
    envolvente convexa (rojo), con el volumen en el título."""
    vertices = points[hull.vertices]
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], marker='.', color='black')
    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], marker='x', color='red')
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def plot_hull_surface(points, hull, title, filename):
    """PNG de superficie 3D sólida de la envolvente convexa: cada cara
    triangular del hull (hull.simplices) coloreada según su altura promedio
    (colormap viridis) — vista tipo malla sólida, complementaria a
    plot_hull_volume() (dispersión de puntos + vértices)."""
    faces = points[hull.simplices]  # (n_faces, 3, 3)
    face_z = faces[:, :, 2].mean(axis=1)
    norm = plt.Normalize(face_z.min(), face_z.max()) if face_z.max() > face_z.min() \
        else plt.Normalize(face_z.min() - 1, face_z.max() + 1)
    colors = plt.cm.viridis(norm(face_z))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    poly = Poly3DCollection(faces, facecolors=colors, edgecolor='k', linewidths=0.3, alpha=0.95)
    ax.add_collection3d(poly)

    # add_collection3d no autoescala los límites de los ejes: hay que fijarlos
    # a mano con el bounding box de los puntos, si no el plot queda vacío.
    mins, maxs = points.min(axis=0), points.max(axis=0)
    ax.set_xlim(mins[0], maxs[0])
    ax.set_ylim(mins[1], maxs[1])
    ax.set_zlim(mins[2], maxs[2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    fig.savefig(filename, dpi=200)
    plt.close(fig)


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
_DF_COLS = ['Pos R', 'Res', 'Atom', 'Dist', 'Lig', 'Type', 'Angle', 'Interaction', 'LigID']


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
# Pockets hidrofóbicos (agrupación multi-residuo por fragmento del ligando)
# ──────────────────────────────────────────────────────────────────────────────

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


def scripting_vmd_pockets(df_detail, Prot, chain, Lig, folder):
    """Genera un .tcl de VMD con una representación de superficie (Surf/MSMS)
    por cada pocket hidrofóbico validado, para visualizar la cavidad que
    envuelve al fragmento del ligando.

    NOTA VMD: el estilo de dibujo Wireframe/Solid Surface/Points de Surf y MSMS
    no es scripteable vía 'mol modstyle' (solo acepta probe radius y
    resolución); hay que cambiarlo a mano en Graphics > Representations >
    Draw style > Wireframe para cada representación agregada por este script."""
    if df_detail.empty:
        return

    receptor_name = Path(Prot).stem
    Lig_name = Path(Lig).stem
    Prot_file = Path(Prot).name
    Lig_file  = Path(Lig).name

    with open(f'{folder}/vmd_pockets_{receptor_name}_{Lig_name}.tcl', 'w') as VDM_TCL:
        # molid real en variable Tcl (ver nota en scripting_vmd): evita romper si la
        # sesión de VMD ya tenía moléculas cargadas antes de correr este script.
        VDM_TCL.write('display projection orthographic\n')
        VDM_TCL.write(f'set molReceptor [mol new "{Prot_file}"]\n')
        VDM_TCL.write('mol modselect 0 $molReceptor all\n')
        VDM_TCL.write('mol modstyle 0 $molReceptor Lines 3\n')
        VDM_TCL.write('mol modcolor 0 $molReceptor ColorID 6\n')
        VDM_TCL.write('mol modmaterial 0 $molReceptor Opaque\n')

        VDM_TCL.write('\n# --- Pockets hidrofobicos: superficie Surf por pocket ---\n')
        VDM_TCL.write('# Cambiar a mano "Draw style" -> Wireframe en Graphics > Representations\n')
        VDM_TCL.write('# para cada representacion Surf agregada abajo.\n')
        rep = 1
        for pocket_n in sorted(df_detail['Pocket'].unique()):
            sub = df_detail[df_detail['Pocket'] == pocket_n]
            residues = ' '.join(sorted({str(p) for p in sub['Pos R']}))
            colorid  = _POCKET_COLORIDS[(int(pocket_n) - 1) % len(_POCKET_COLORIDS)]
            VDM_TCL.write('mol addrep $molReceptor\n')
            VDM_TCL.write(f'mol modselect {rep} $molReceptor "resid {residues} and chain {chain}"\n')
            VDM_TCL.write(f'mol modstyle {rep} $molReceptor Surf 1.4 0\n')
            VDM_TCL.write(f'mol modcolor {rep} $molReceptor ColorID {colorid}\n')
            VDM_TCL.write(f'mol modmaterial {rep} $molReceptor Opaque\n')
            rep += 1

        VDM_TCL.write(f'\nset molLigand [mol new "{Lig_file}"]\n')
        VDM_TCL.write('mol addrep $molLigand\n')
        VDM_TCL.write('mol modstyle 0 $molLigand Licorice\n')
        VDM_TCL.write('display resetview\n')


def scripting_vmd_combined(DF_Interacciones, receptor_points, DF_Lig, df_pocket_detail,
                            Prot, chain, Lig, folder):
    """Una sola escena de VMD con todo junto: H-bonds y aromáticas (líneas
    punteadas, igual que scripting_vmd) + superficie Surf por pocket
    hidrofóbico validado (igual que scripting_vmd_pockets) — para no tener que
    cargar dos .tcl por separado y comparar a ojo. No incluye las líneas
    hidrofóbicas individuales de scripting_vmd_hydrophobic (la superficie del
    pocket ya representa esa región; agregarlas encima satura la escena)."""
    receptor_name = Path(Prot).stem
    Lig_name = Path(Lig).stem
    Prot_file = Path(Prot).name
    Lig_file  = Path(Lig).name

    Res_All = DF_Interacciones['Pos R'].tolist()
    residues = ' '.join(map(str, Res_All))

    with open(f'{folder}/vmd_combined_{receptor_name}_{Lig_name}.tcl', 'w') as VDM_TCL:
        VDM_TCL.write('display projection orthographic\n')
        VDM_TCL.write(f'set molReceptor [mol new "{Prot_file}"]\n')
        VDM_TCL.write('mol modselect 0 $molReceptor all\n')
        VDM_TCL.write('mol modstyle 0 $molReceptor Lines 3\n')
        VDM_TCL.write('mol modcolor 0 $molReceptor ColorID 6\n')
        VDM_TCL.write('mol modmaterial 0 $molReceptor Opaque\n')
        VDM_TCL.write('mol addrep $molReceptor\n')
        VDM_TCL.write(f'mol modselect 1 $molReceptor resid {residues} and chain {chain}\n')
        VDM_TCL.write('mol modstyle 1 $molReceptor Licorice\n')

        rep = 2
        if not df_pocket_detail.empty:
            VDM_TCL.write('\n# --- Pockets hidrofobicos: superficie Surf por pocket ---\n')
            VDM_TCL.write('# Cambiar a mano "Draw style" -> Wireframe en Graphics > Representations\n')
            VDM_TCL.write('# para cada representacion Surf agregada abajo.\n')
            for pocket_n in sorted(df_pocket_detail['Pocket'].unique()):
                sub = df_pocket_detail[df_pocket_detail['Pocket'] == pocket_n]
                pocket_residues = ' '.join(sorted({str(p) for p in sub['Pos R']}))
                colorid = _POCKET_COLORIDS[(int(pocket_n) - 1) % len(_POCKET_COLORIDS)]
                VDM_TCL.write('mol addrep $molReceptor\n')
                VDM_TCL.write(f'mol modselect {rep} $molReceptor "resid {pocket_residues} and chain {chain}"\n')
                VDM_TCL.write(f'mol modstyle {rep} $molReceptor Surf 1.4 0\n')
                VDM_TCL.write(f'mol modcolor {rep} $molReceptor ColorID {colorid}\n')
                VDM_TCL.write(f'mol modmaterial {rep} $molReceptor Opaque\n')
                rep += 1

        VDM_TCL.write(f'\nset molLigand [mol new "{Lig_file}"]\n')
        VDM_TCL.write('mol addrep $molLigand\n')
        VDM_TCL.write('mol modstyle 0 $molLigand Licorice\n')
        VDM_TCL.write('display resetview\n')

        VDM_TCL.write('\n# --- Interacciones H-bond / aromaticas ---\n')
        _write_hbond_aromatic_lines(VDM_TCL, DF_Interacciones, receptor_points, DF_Lig, chain)


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


# ──────────────────────────────────────────────────────────────────────────────
# Análisis de un par receptor-ligando
# ──────────────────────────────────────────────────────────────────────────────

_ALL_TYPES = ['acceptor', 'donor', 'aromatic', 'hydrophobic', 'salt_bridge', 'pi_cation']


def analyze_pair(receptor_pdb, Ligand_imput, chain_receptor, cfg):
    """Ejecuta el pipeline completo para un par receptor-ligando."""
    ligand_plot           = cfg['ligand_plot']
    vmd_output            = cfg['vmd_output']
    cumulative_output     = cfg['cumulative_output']
    Interaction_Coord_Source = cfg['Interaction_Coord_Source']
    Volume_Plot           = cfg['Volume_Plot']
    Bias                  = cfg['Bias']
    Bias_Validated_Only   = cfg['Bias_Validated_Only']
    Distances_Hidrogen_Bonds = cfg['Distances_Hidrogen_Bonds']
    Distances_Aromatic    = cfg['Distances_Aromatic']
    Distancia_Hidrofobica = cfg['Distancia_Hidrofobica']
    Distancia_Centro_Activo = cfg['Distancia_Centro_Activo']
    Angle_Hidrogen_Bonds_Min = cfg['Angle_Hidrogen_Bonds_Min']
    Angle_Hidrogen_Bonds_Max = cfg['Angle_Hidrogen_Bonds_Max']
    Ring_Planarity_RMSD_Max = cfg['Ring_Planarity_RMSD_Max']
    Pocket_Min_Residues   = cfg['Pocket_Min_Residues']
    Pocket_Coverage_Threshold = cfg['Pocket_Coverage_Threshold']
    Pocket_Density_Radius = cfg['Pocket_Density_Radius']
    Aceptores_Prot        = cfg['Aceptores_Prot']
    Dadores_Prot          = cfg['Dadores_Prot']
    Aceptot_antecedent    = cfg['Aceptot_antecedent']

    threshold_PH            = 4
    numero_anillo_aromatico = 5

    receptor = Path(receptor_pdb).stem
    ligand   = Path(Ligand_imput).stem
    folder   = f'{receptor}_{ligand}'
    Path(folder).mkdir(exist_ok=True)

    # ── Limpieza del ligando ──────────────────────────────────────
    remove_bias(Ligand_imput, folder)

    # ── Ligando: hot-points ───────────────────────────────────────
    mol = Chem.MolFromPDBFile(Ligand_imput, removeHs=False)
    if mol is None:
        print(f"  [WARN] RDKit could not read ligand: {Ligand_imput}")
        return

    pdb_coords, CM = extract_coords_from_pdb(Ligand_imput)
    acceptor_atoms, donor_atoms = search_hot_points(Ligand_imput, mol, pdb_coords, ligand_plot, folder)

    # ── Anillos aromáticos ────────────────────────────────────────
    aromatic_rings_data, rings_data = search_rings(mol, pdb_coords, numero_anillo_aromatico,
                                                    Ring_Planarity_RMSD_Max)
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
        distances       = np.linalg.norm(Matriz_receptor - ring_center, axis=1)
        for idx in np.where(distances < Distances_Aromatic)[0]:
            closest = Sub_Set_Receptor.iloc[idx]
            DF_Interacciones.loc[len(DF_Interacciones)] = [
                closest.iloc[1], closest.iloc[2], closest.iloc[3],
                distances[idx], Sub_Set_Ligando.iloc[0, 4], 'aromatic', 0.0, 0, np.nan]

    DF_Lig_All = generate_df_ligand(pdb_coords)
    DF_Interacciones = DF_Interacciones.drop_duplicates()

    # ── Nuevos tipos de interacción ───────────────────────────────
    df_hpho = search_hydrophobic(mol, pdb_coords, DF_Active_Site, Distancia_Hidrofobica)
    df_salt = search_salt_bridges(mol, pdb_coords, DF_Active_Site)
    df_pica = search_pi_cation(DF_Active_Site, aromatic_lig_df)
    df_pocket_summary, df_pocket_detail, pocket_hulls = search_hydrophobic_pockets(
        mol, pdb_coords, DF_Active_Site, Distancia_Hidrofobica,
        Pocket_Min_Residues, Pocket_Coverage_Threshold, Pocket_Density_Radius)
    # Frames vacíos (sin matches) no tienen dtypes declarados (columns=_DF_COLS sin data);
    # concatenarlos junto con DF_Interacciones (tipado) dispara el FutureWarning de pandas
    # sobre inferencia de dtype sobre columnas vacías/all-NA. Se excluyen antes de concatenar.
    frames = [f for f in (DF_Interacciones, df_hpho, df_salt, df_pica) if not f.empty]
    if frames:
        DF_Interacciones = pd.concat(frames, ignore_index=True).drop_duplicates()

    # ── Validación por ángulo ─────────────────────────────────────
    for j in range(DF_Interacciones.shape[0]):
        tipo = DF_Interacciones.iloc[j, 5]
        if tipo == 'acceptor':
            Aceptor_Antecedent = Busqueda_Antecesor_Lig(DF_Interacciones.iloc[j, 8], DF_Lig_All)
            Aceptor  = np.array(DF_Lig[DF_Lig['Atom ID'] == DF_Interacciones.iloc[j,8]].iloc[0, [1,2,3]])
            resultado = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0]) &
                                       (DF_Active_Site['Atom'] == DF_Interacciones.iloc[j,2])]
            Donor = np.array(resultado[['X','Y','Z']]).reshape(-1)
            DF_Interacciones.iloc[j, 6] = float(angle_three_points(Donor, Aceptor, Aceptor_Antecedent))
        elif tipo == 'donor':
            Donor    = np.array(DF_Lig[DF_Lig['Atom ID'] == DF_Interacciones.iloc[j,8]].iloc[0, [1,2,3]])
            resultado = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0]) &
                                       (DF_Active_Site['Atom'] == DF_Interacciones.iloc[j,2])]
            Aceptor  = np.array(resultado[['X','Y','Z']]).reshape(-1)
            try:
                Atomo = Aceptot_antecedent[DF_Interacciones.iloc[j,1]][DF_Interacciones.iloc[j,2]]
                resultado = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0]) &
                                           (DF_Active_Site['Atom'] == Atomo)]
                Aceptor_Antecedent = np.array(resultado[['X','Y','Z']]).reshape(-1)
            except KeyError:
                resultado = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0]) &
                                           (DF_Active_Site['Atom'] == 'C')]
                Aceptor_Antecedent = np.array(resultado[['X','Y','Z']]).reshape(-1)
            DF_Interacciones.iloc[j, 6] = float(angle_three_points(Donor, Aceptor, Aceptor_Antecedent))
        elif tipo == 'aromatic':
            Anillo_Proteina = DF_Active_Site[DF_Active_Site['Pos'] == DF_Interacciones.iloc[j, 0]]
            Anillo_Lig      = DF_Lig[DF_Lig['Caso'] == DF_Interacciones.iloc[j, 4]]
            DF_Interacciones.iloc[j, 6] = Interaccion_Aromatica(Anillo_Proteina, Anillo_Lig)

    # ── Clasificación final ───────────────────────────────────────
    for k in range(DF_Interacciones.shape[0]):
        tipo = DF_Interacciones.iloc[k, 5]
        dist = float(DF_Interacciones.iloc[k, 3])
        ang  = float(DF_Interacciones.iloc[k, 6])
        if tipo in ('hydrophobic', 'salt_bridge', 'pi_cation'):
            pass  # validadas en sus funciones con criterio de distancia
        elif tipo in ('acceptor', 'donor'):
            DF_Interacciones.iloc[k, 7] = (
                'Yes' if dist < Distances_Hidrogen_Bonds
                and Angle_Hidrogen_Bonds_Min < ang <= Angle_Hidrogen_Bonds_Max else 'No')
        elif tipo == 'aromatic':
            if dist < Distances_Aromatic:
                # parallel/sandwich: 0-30°  |  T-shaped: 60-90°
                DF_Interacciones.iloc[k, 7] = 'Yes' if (ang < 30 or ang > 60) else 'No'
            else:
                DF_Interacciones.iloc[k, 7] = 'No'

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
    DF_dist = DF_Interacciones[DF_Interacciones['Dist'] < Distances_Aromatic]
    DF_dist.drop(columns=['LigID']).to_csv(f'{folder}/Interaction_{receptor}_{ligand}_threshold.csv')
    DF_true = DF_Interacciones[DF_Interacciones['Interaction'] == 'Yes']
    DF_true.drop(columns=['LigID']).to_csv(f'{folder}/Interaction_{receptor}_{ligand}_true.csv')

    if Bias == 'Yes':
        bias_df_true = DF_true if Bias_Validated_Only == 'Yes' else None
        export_bpf(DF_Lig, f'{folder}/{receptor}_{ligand}.bpf', bias_df_true)
        export_bpf_pdb(DF_Lig, f'{folder}/{receptor}_{ligand}_bias.pdb', bias_df_true)

    df_pocket_summary.to_csv(f'{folder}/Pockets_{receptor}_{ligand}.csv', index=False)

    if Volume_Plot == 'Yes':
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


#### Busqueda de interacciones ####


def main():

    parser = argparse.ArgumentParser(
        description='Análisis de interacciones proteína-ligando.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            'Modos de uso:\n'
            '  Un ligando   : -r proteina.pdb -l ligando.pdb -c A\n'
            '  Batch        : -r proteina.pdb -l lig1.pdb lig2.pdb lig3.pdb -c A\n'
            '  PDB complejo : -x complejo.pdb -c A\n'
            '                 -x complejo.pdb -c A -n LIG\n'
            '  Sin HETATM   : -x complejo.pdb -c A -f TF3 (ligando guardado como ATOM)\n'
        )
    )
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument('-x', '--complex', dest='complex_pdb', metavar='COMPLEX.pdb',
                     help='PDB complejo. Se separa automáticamente.')
    grp.add_argument('-r', '--receptor_pdb', default=None,
                     help='PDB del receptor.')
    parser.add_argument('-l', '--ligand_input', nargs='+', default=None,
                        help='PDB(s) del ligando. Acepta múltiples para análisis batch.')
    parser.add_argument('-c', '--chain_receptor', required=True,
                        help='Cadena de la proteína.')
    parser.add_argument('-n', '--lig_name', default=None,
                        help='Nombre del HETATM a usar como ligando (con --complex).')
    parser.add_argument('-f', '--force_ligand', nargs='+', default=None,
                        help='Resname(s) a tratar como ligando aunque figuren como ATOM '
                             'en vez de HETATM en el PDB complejo (ej: -f TF3 7FW).')
    parser.add_argument('--config', default=None, metavar='CONFIG.yml',
                        help='Ruta al archivo YAML de configuración '
                             '(por defecto: Interacciones_variables.yml en la raíz del proyecto).')

    args = parser.parse_args()

    # ── Resolver lista de pares (receptor, ligando) ───────────────
    import tempfile
    pairs    = []
    tmp_dir  = None   # directorio temporal para --complex, se borra al final

    if args.complex_pdb:
        print(f"\n[Split] Splitting: {args.complex_pdb}")
        tmp_dir = tempfile.mkdtemp(prefix='interactions_split_')
        protein_path, het_paths = split_pdb(args.complex_pdb, output_dir=tmp_dir,
                                            force_ligand_names=args.force_ligand)
        print(f"  Protein  -> {protein_path}")
        if not het_paths:
            print("  No HETATM groups found (water excluded).")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            sys.exit(1)
        print(f"  HETATM   -> {list(het_paths.keys())}")
        if args.lig_name:
            if args.lig_name not in het_paths:
                print(f"  Error: '{args.lig_name}' not found. Available: {list(het_paths.keys())}")
                shutil.rmtree(tmp_dir, ignore_errors=True)
                sys.exit(1)
            pairs = [(str(protein_path), str(het_paths[args.lig_name]))]
        elif len(het_paths) == 1:
            resname, lig_path = next(iter(het_paths.items()))
            print(f"  Selected: {resname}")
            pairs = [(str(protein_path), str(lig_path))]
        else:
            print(f"  Multiple HETATM groups. Use -n to select: {list(het_paths.keys())}")
            shutil.rmtree(tmp_dir, ignore_errors=True)
            sys.exit(1)
    else:
        if not args.ligand_input:
            parser.error("Con -r debes proveer también -l/--ligand_input.")
        pairs = [(args.receptor_pdb, lig) for lig in args.ligand_input]

    # ── Cargar configuración una sola vez ─────────────────────────
    (ligand_plot, vmd_output, cumulative_output, Interaction_Coord_Source, Volume_Plot, Bias,
     Bias_Validated_Only, Distances_Hidrogen_Bonds, Distances_Aromatic,
     Distancia_Hidrofobica, Distancia_Centro_Activo, Angle_Hidrogen_Bonds_Min,
     Angle_Hidrogen_Bonds_Max, Ring_Planarity_RMSD_Max, Pocket_Min_Residues,
     Pocket_Coverage_Threshold, Pocket_Density_Radius, Aceptores_Prot, Dadores_Prot,
     Aceptot_antecedent, Special_case) = carga_variables(args.config)

    cfg = {
        'ligand_plot':              ligand_plot,
        'vmd_output':               vmd_output,
        'cumulative_output':        cumulative_output,
        'Interaction_Coord_Source': Interaction_Coord_Source,
        'Volume_Plot':              Volume_Plot,
        'Bias':                     Bias,
        'Bias_Validated_Only':      Bias_Validated_Only,
        'Distances_Hidrogen_Bonds': Distances_Hidrogen_Bonds,
        'Distances_Aromatic':       Distances_Aromatic,
        'Distancia_Hidrofobica':    Distancia_Hidrofobica,
        'Distancia_Centro_Activo':  Distancia_Centro_Activo,
        'Angle_Hidrogen_Bonds_Min': Angle_Hidrogen_Bonds_Min,
        'Angle_Hidrogen_Bonds_Max': Angle_Hidrogen_Bonds_Max,
        'Ring_Planarity_RMSD_Max':  Ring_Planarity_RMSD_Max,
        'Pocket_Min_Residues':      Pocket_Min_Residues,
        'Pocket_Coverage_Threshold': Pocket_Coverage_Threshold,
        'Pocket_Density_Radius':    Pocket_Density_Radius,
        'Aceptores_Prot':           Aceptores_Prot,
        'Dadores_Prot':             Dadores_Prot,
        'Aceptot_antecedent':       Aceptot_antecedent,
        'Special_case':             Special_case,
    }

    # ── Análisis (uno o batch) ────────────────────────────────────
    n_ok, n_skip = 0, 0
    for receptor_pdb, Ligand_imput in pairs:
        print(f"\n{'─'*60}")
        print(f"  Receptor : {receptor_pdb}")
        print(f"  Ligand   : {Ligand_imput}")
        print(f"  Chain    : {args.chain_receptor}")
        print(f"{'─'*60}")
        errors = validate_inputs(receptor_pdb, Ligand_imput, args.chain_receptor)
        if errors:
            for e in errors:
                print(f"  [ERROR] {e}")
            print("  Skipping this pair.")
            n_skip += 1
            continue
        analyze_pair(receptor_pdb, Ligand_imput, args.chain_receptor, cfg)
        n_ok += 1

    # ── Limpiar directorio temporal del split ─────────────────────
    if tmp_dir:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"\nAnalysis complete: {n_ok} pair(s) processed, {n_skip} skipped.")


if __name__ == '__main__':
    main()
