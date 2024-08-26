### Librerias ###

import yaml
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from Bio.PDB import *
import math
import argparse
import os
import shutil


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

# Función para obtener coordenadas por índice de átomo
def get_coord_by_atom_id(atom_id):
    for coord in pdb_coords:
        if coord[1] == atom_id:
            return f"{atom_id}:{coord[5]}:{coord[6]}:{coord[7]}"
    return None


def search_hot_points(Ligand_imput):

    # Definir los patrones SMARTS para cada caso
    acceptor_smarts = ['[O;H1]', '[O;H0]', '[N;H1]', '[N;H0]' , '[n]' , '[o]' , '[N+]']  # Aceptores
    donor_smarts = ['[O;H]', '[N;H]', '[S;H]' ,'[nH]']  # Donadores
    aromatic_smarts = 'a'  # Aromáticos

    # Diccionarios para almacenar las coordenadas
    acceptor_coords = []
    donor_coords = []
    aromatic_coords = []

    # Buscar y almacenar las coordenadas para aceptores de puente de hidrógeno
    acceptor_atoms = []
    for smarts in acceptor_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            for atom_idx in match:
                acceptor_atoms.append(atom_idx)
                acceptor_coords.append(get_atom_coords(mol, atom_idx))

    # Buscar y almacenar las coordenadas para donadores de puente de hidrógeno
    donor_atoms = []
    for smarts in donor_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            for atom_idx in match:
                donor_atoms.append(atom_idx)
                donor_coords.append(get_atom_coords(mol, atom_idx))

    # Crear gráficos separados para cada caso y guardarlos como imágenes
    # Gráfico para aceptores
    if ligand_plot == 'Yes':
        mol_copy = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(mol_copy)
        img_acceptors = Draw.MolToImage(mol_copy, highlightAtoms=acceptor_atoms, size=(600, 600))
        img_acceptors.save(f"{folder}/{Ligand_imput.split('.')[0]}_acceptors.png")

        # Gráfico para donadores
        mol_copy = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(mol_copy)
        img_donors = Draw.MolToImage(mol_copy, highlightAtoms=donor_atoms, size=(600, 600))
        img_donors.save(f"{folder}/{Ligand_imput.split('.')[0]}_donors.png")


    return(acceptor_atoms,donor_atoms)


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


def carga_variables():
    # Cargo Variables Generales #
    with open(r'Interacciones_variables.yml') as file:
        Interaciones = yaml.load(file, Loader=yaml.FullLoader)

    ligand_plot = str(Interaciones['options']['ligand_plot'])
    vmd_output = str(Interaciones['options']['vmd_output'])

    Distances_Hidrogen_Bonds =float(Interaciones['distancias']['Distances_Hidrogen_Bonds'])
    Distances_Aromatic = float(Interaciones['distancias']['Distances_Aromatic'])
    Distancia_Hidrofobica = float(Interaciones['distancias']['Distances_Hidrofobica'])
    
    Aceptores_Prot = Interaciones['acceptors']
    Dadores_Prot = Interaciones['donors']
    Aceptot_antecedent = Interaciones['acceptors_antecedent']
    Special_case = Interaciones['special']
    return(ligand_plot,vmd_output,Distances_Hidrogen_Bonds,Distances_Aromatic,Distancia_Hidrofobica,Aceptores_Prot,Dadores_Prot,Aceptot_antecedent,Special_case)



def Coordenadas_interes_receptor(Aceptores_Prot,Dadores_Prot,DF_Active_Site):
    ### Obtengo las coordenadas de los atomos de interes en el receptor
    receptor_points = pd.DataFrame(columns=['Type','Pos','Residue', 'Atom', 'X' , 'Y' , 'Z'])
    for pos in range(0,DF_Active_Site.shape[0]):
        Atomo = (DF_Active_Site.iloc[pos,2])
        Res = (DF_Active_Site.iloc[pos,3])
        listado = (Aceptores_Prot[Atomo])
        if Res in listado:
            receptor_points.loc[len(receptor_points.index)] = 'Aceptor',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],DF_Active_Site.iloc[pos,3],DF_Active_Site.iloc[pos,4],DF_Active_Site.iloc[pos,5],DF_Active_Site.iloc[pos,6]
    for pos in range(0,DF_Active_Site.shape[0]):
        Atomo = (DF_Active_Site.iloc[pos,2])
        Res = (DF_Active_Site.iloc[pos,3])
        listado = (Dadores_Prot[Atomo])
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

	for  j in range(0,len(Aromatic_Ring)):
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
                Sub_Set_Ligando.iloc[j, 0],  # Identificador del ligando
                Lig_Caso,               # Caso del ligando
                0, 0                    # Valores adicionales
            ]
    return(DF_Interacciones)

def generate_df_ligand(pdb_coords):
    
    columns = ['Atom ID', 'Element', 'Residue Name', 'Chain ID', 'Residue Number', 'X', 'Y', 'Z']
    df_ligand = pd.DataFrame(pdb_coords, columns=columns)

    return(df_ligand)



def Busqueda_Antecesor_Lig(Atomo_ID,Lig_DF):
    # Punto dado
    punto_dado = np.array(Lig_DF.query('Element == @Atomo_ID')[['X' , 'Y' , 'Z']])
    # Filtrar átomos que no sean H
    df_filtrado = Lig_DF[~Lig_DF['Element'].str.contains('H')]
    df_filtrado = df_filtrado[~df_filtrado['Element'].str.contains(Atomo_ID)]
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

def search_rings(mol):
    # Kekulizar la molécula para asegurar que los anillos se detecten correctamente
    Chem.Kekulize(mol, clearAromaticFlags=True)
    
    # Obtener información de los anillos
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()
    
    # Listar los anillos encontrados
    ring_data = []
    for i, ring in enumerate(ring_atoms):
        ring_data.append({
            'Ring': i + 1,
            'Atoms': ring,
            'Ring Size': len(ring)
        })

    rings_data = []
    # Imprimir información de los anillos
    for ring in ring_data:
        #### Filtro por el numero de anillos que busco como aromatico
        if ring['Ring Size'] > numero_anillo_aromatico:
            for atom in ring['Atoms']:
                rings_data.append([pdb_coords[atom][1] ,pdb_coords[atom][5] , pdb_coords[atom][6] ,pdb_coords[atom][7],'aromatic '+str(ring['Ring'])+' (#'+str(ring['Ring Size'])+')'])

   
    return ring_data,rings_data

def visualize_rings(mol, ring_data, Ligand_imput):
    # Dibujar la molécula y resaltar los átomos de los anillos
    ring_atoms = [atom for ring in ring_data for atom in ring['Atoms']]
    mol_copy = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mol_copy)
    img = Draw.MolToImage(mol_copy, highlightAtoms=ring_atoms, size=(600, 600), kekulize=True)
    img.save(f"{folder}/{Ligand_imput.split('.')[0]}_aromatic.png")
    

def scripting_vmd(DF_Interacciones,receptor_points,aromatic_lig_df,DF_Lig,Prot,chain,Lig):
    receptor_name = Prot.split('.')[0]
    Lig_name = Lig.split('.')[0]
    
    Res_All = DF_Interacciones['Pos R'].tolist()
    residues = ' '.join(map(str, Res_All))

    VDM_TCL = open(f'{folder}/vmd_{receptor_name}_{Lig_name}.tcl' , 'w')
    
    # Cargar el archivo PDB
    VDM_TCL.write(f'display projection orthographic\n')
    VDM_TCL.write(f'mol new "{Prot}"\n')
    VDM_TCL.write(f'mol modselect 0 0 all\n')
    VDM_TCL.write(f'mol modstyle 0 0 NewCartoon\n')
    VDM_TCL.write(f'mol modcolor 0 0 ColorID 6\n')
    VDM_TCL.write(f'mol modmaterial 0 0 Transparent\n')
    VDM_TCL.write(f'mol addrep top\n')
    VDM_TCL.write(f'mol modselect 1 top resid {residues} and chain {chain}\n')
    VDM_TCL.write(f'mol modstyle 1 top Licorice\n')
    VDM_TCL.write(f'mol new {Lig}\n')
    # Crear una representación en Licorice para el ligando
    VDM_TCL.write(f'mol addrep 1\n')
    VDM_TCL.write(f'mol modstyle 0 1 Licorice\n')
    
    ### Busco Interaccion
    for j in range(0,DF_Interacciones.shape[0]):
        if DF_Interacciones.iloc[j,5] == 'aromatic':
            Recept = DF_Interacciones.iloc[j,0]
            Coord2 = np.array((receptor_points[(receptor_points['Pos'] == Recept) & (receptor_points['Atom'] == 'center')])[['X','Y','Z']])[0]
            # ligando #
            anillo = DF_Interacciones.iloc[j,4]
            punto_dado = np.array(DF_Lig.query('Caso == @anillo')[['Coord X' , 'Coord Y' , 'Coord Z']])
            Coord1 = np.mean(punto_dado, axis=0) 
            # Receptor #
            VDM_TCL.write('graphics top color white\n')
            VDM_TCL.write('graphics top line {'+str(Coord1[0])+' '+str(Coord1[1])+' '+str(Coord1[2])+'} {'+str(Coord2[0])+' '+str(Coord2[1])+' '+str(Coord2[2])+'} width 5 style dashed\n')
            VDM_TCL.write(f'set Recptor{j} [atomselect 0 "chain {chain} and resid {DF_Interacciones.iloc[j,0]} and resname {DF_Interacciones.iloc[j,1]} and name CZ"]\n')
            VDM_TCL.write('set x1 {'+str(Coord1[0])+'}\n')
            VDM_TCL.write('set y1 {'+str(Coord1[1])+'}\n')
            VDM_TCL.write('set z1 {'+str(Coord1[2])+'}\n')
            VDM_TCL.write('set x2 {'+str(Coord2[0])+'}\n')
            VDM_TCL.write('set y2 {'+str(Coord2[1])+'}\n')
            VDM_TCL.write('set z2 {'+str(Coord2[2])+'}\n')
            VDM_TCL.write('set dx [expr {$x1 - $x2}]\n')
            VDM_TCL.write('set dy [expr {$y1 - $y2}]\n')
            VDM_TCL.write('set dz [expr {$z1 - $z2}]\n')
            VDM_TCL.write('set distance [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]\n')
            VDM_TCL.write('set xm [expr {($x1 + $x2) / 2}]\n')
            VDM_TCL.write('set ym [expr {($y1 + $y2) / 2}]\n')
            VDM_TCL.write('set zm [expr {($z1 + $z2) / 2}]\n')
            VDM_TCL.write('graphics top color white\n')
            VDM_TCL.write('graphics top text [list $xm $ym $zm] [format "%.2f A" $distance]\n')
        if DF_Interacciones.iloc[j,5] == 'acceptor':
            # Receptor #
            # busco donnor #
            Recept = DF_Interacciones.iloc[j,0]
            Coord2 = np.array(receptor_points[(receptor_points['Pos'] == Recept) & (receptor_points['Atom'] == DF_Interacciones.iloc[j,2])][['X','Y','Z']])[0]
            # Ligand #
            atomo = DF_Interacciones.iloc[j,4]
            Coord1 = np.array(DF_Lig[(DF_Lig['Átomo'] == atomo)][['Coord X' , 'Coord Y' , 'Coord Z']])[0]
            VDM_TCL.write('graphics top color red\n')
            VDM_TCL.write('graphics top line {'+str(Coord1[0])+' '+str(Coord1[1])+' '+str(Coord1[2])+'} {'+str(Coord2[0])+' '+str(Coord2[1])+' '+str(Coord2[2])+'} width 5 style dashed\n')
            VDM_TCL.write(f'set Recptor{j} [atomselect 0 "chain {chain} and resid {DF_Interacciones.iloc[j,0]} and resname {DF_Interacciones.iloc[j,1]} and name CZ"]\n')
            VDM_TCL.write('set x1 {'+str(Coord1[0])+'}\n')
            VDM_TCL.write('set y1 {'+str(Coord1[1])+'}\n')
            VDM_TCL.write('set z1 {'+str(Coord1[2])+'}\n')
            VDM_TCL.write('set x2 {'+str(Coord2[0])+'}\n')
            VDM_TCL.write('set y2 {'+str(Coord2[1])+'}\n')
            VDM_TCL.write('set z2 {'+str(Coord2[2])+'}\n')
            VDM_TCL.write('set dx [expr {$x1 - $x2}]\n')
            VDM_TCL.write('set dy [expr {$y1 - $y2}]\n')
            VDM_TCL.write('set dz [expr {$z1 - $z2}]\n')
            VDM_TCL.write('set distance [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]\n')
            VDM_TCL.write('set xm [expr {($x1 + $x2) / 2}]\n')
            VDM_TCL.write('set ym [expr {($y1 + $y2) / 2}]\n')
            VDM_TCL.write('set zm [expr {($z1 + $z2) / 2}]\n')
            VDM_TCL.write('graphics top color white\n')
            VDM_TCL.write('graphics top text [list $xm $ym $zm] [format "%.2f A" $distance]\n')
        if DF_Interacciones.iloc[j,5] == 'donnor':
            # Receptor #
            # busco donnor #
            Recept = DF_Interacciones.iloc[j,0]
            Coord2 = np.array(receptor_points[(receptor_points['Pos'] == Recept) & (receptor_points['Atom'] == DF_Interacciones.iloc[j,2])][['X','Y','Z']])[0]
            # Ligand #
            atomo = DF_Interacciones.iloc[j,4]
            Coord1 = np.array(DF_Lig[(DF_Lig['Átomo'] == atomo)][['Coord X' , 'Coord Y' , 'Coord Z']])[0]
            VDM_TCL.write('graphics top color yellow\n')
            VDM_TCL.write('graphics top line {'+str(Coord1[0])+' '+str(Coord1[1])+' '+str(Coord1[2])+'} {'+str(Coord2[0])+' '+str(Coord2[1])+' '+str(Coord2[2])+'} width 5 style dashed\n')
            VDM_TCL.write(f'set Recptor{j} [atomselect 0 "chain {chain} and resid {DF_Interacciones.iloc[j,0]} and resname {DF_Interacciones.iloc[j,1]} and name CZ"]\n')
            VDM_TCL.write('set x1 {'+str(Coord1[0])+'}\n')
            VDM_TCL.write('set y1 {'+str(Coord1[1])+'}\n')
            VDM_TCL.write('set z1 {'+str(Coord1[2])+'}\n')
            VDM_TCL.write('set x2 {'+str(Coord2[0])+'}\n')
            VDM_TCL.write('set y2 {'+str(Coord2[1])+'}\n')
            VDM_TCL.write('set z2 {'+str(Coord2[2])+'}\n')
            VDM_TCL.write('set dx [expr {$x1 - $x2}]\n')
            VDM_TCL.write('set dy [expr {$y1 - $y2}]\n')
            VDM_TCL.write('set dz [expr {$z1 - $z2}]\n')
            VDM_TCL.write('set distance [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]\n')
            VDM_TCL.write('set xm [expr {($x1 + $x2) / 2}]\n')
            VDM_TCL.write('set ym [expr {($y1 + $y2) / 2}]\n')
            VDM_TCL.write('set zm [expr {($z1 + $z2) / 2}]\n')
            VDM_TCL.write('graphics top color white\n')
            VDM_TCL.write('graphics top text [list $xm $ym $zm] [format "%.2f A" $distance]\n')
    VDM_TCL.close()


def remove_bias(file_path):
    old_file_path = file_path.replace('.pdb', '_old.pdb')
    shutil.copy(file_path, folder+'/'+old_file_path)

    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Filtrar las líneas que no contienen "CM"
    lines = [line for line in lines if ' CM ' not in line]

    # Guardar el archivo sin las líneas "CM"
    with open(file_path, 'w') as f:
        f.writelines(lines)





#### Busqueda de interacciones ####


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script para evaluar interacciones.')
    
    parser.add_argument('-r', '--receptor_pdb', required=True, help='Archivo PDB del receptor.')
    parser.add_argument('-c', '--chain_receptor', required=True, help='Archivo PDB del receptor.')
    parser.add_argument('-l', '--ligand_input', required=True, help='Archivo PDB del ligando.')

    args = parser.parse_args()

    # Ligand
    Ligand_imput = args.ligand_input
    # Receptor
    receptor_pdb = args.receptor_pdb
    chain_receptor= args.chain_receptor
    # SA
    Distancia_Centro_Activo = 12
    threshold_inter_aro = 5.5
    threshold_PH = 4
    numero_anillo_aromatico = 5

    ligand_plot,vmd_output,Distances_Hidrogen_Bonds,Distances_Aromatic,Distancia_Hidrofobica,Aceptores_Prot,Dadores_Prot,Aceptot_antecedent,Special_case = carga_variables()
    Caso = ['acceptors','donors']

    ## Crear Carpetas ##
    try:
        os.mkdir('{}_{}'.format(receptor_pdb.split('.')[0],Ligand_imput.split('.')[0]))
    except FileExistsError:
        pass

    folder = '{}_{}'.format(receptor_pdb.split('.')[0],Ligand_imput.split('.')[0])

    ### Elimino bias si el ligando lo tiene ###

    remove_bias(Ligand_imput)

    ### Busqueda de hot point en Ligando ###

    # Cargar la molécula desde un archivo PDB
    mol = Chem.MolFromPDBFile(f'{Ligand_imput}', removeHs=False)

    
    acceptor_atoms,donor_atoms = search_hot_points(Ligand_imput)
    ### Busco Coordenas 3D ###

    # Obtener coordenadas desde el archivo PDB
    # Y CM del ligado para establecer SA
    pdb_coords , CM = extract_coords_from_pdb(Ligand_imput)
    ###### Busqueda Aromaticos ##### 
    

    aromatic_rings_data , rings_data = search_rings(mol)

    if ligand_plot == 'Yes':
        visualize_rings(mol, aromatic_rings_data,Ligand_imput)

    
    DF_Aro = pd.DataFrame(rings_data,columns=['Átomo', 'Coord X', 'Coord Y', 'Coord Z', 'Caso'])
    
    ###### Busqueda Aceptores #####
    coordenadas = []
    
    for j in range(0,len(acceptor_atoms)):
        Coord = (get_coord_by_atom_id(pdb_coords[acceptor_atoms[j]][1]))
        temp = Coord.split(':')
        coordenadas.append([temp[0],float(temp[1]),float(temp[2]),float(temp[3]),'acceptor'])

    ###### Busqueda Dadores #####
    for j in range(0,len(donor_atoms)):
        Coord = (get_coord_by_atom_id(pdb_coords[donor_atoms[j]][1]))
        temp = Coord.split(':')
        coordenadas.append([temp[0],float(temp[1]),float(temp[2]),float(temp[3]),'donnor'])

    DF_Lig = pd.DataFrame(coordenadas,columns=['Átomo', 'Coord X', 'Coord Y', 'Coord Z', 'Caso'])
    
    
    DF_Lig = pd.concat([DF_Lig, DF_Aro], ignore_index=True)

    ### Fin Ligando ####

    ### Inicio Receptor ###

    ## Load the PDB file of receptor ##

    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('pdb', receptor_pdb)

    DF_Active_Site = active_site_residues(structure, CM,chain_receptor, Distancia_Centro_Activo , Ligand_imput.split('.')[0])

    #### Interacciones #####
    
    receptor_points = Coordenadas_interes_receptor(Aceptores_Prot,Dadores_Prot,DF_Active_Site)

    DF_Interacciones = pd.DataFrame( columns=['Pos R','Res','Atom' ,'Dist','Lig','Type', 'Angle' , 'Interaction'])

    ### Caso Ligando Aceptor - Receptor Dador
    
    
    Receptor_Caso ='Dador'
    Lig_Caso = 'acceptor'

    DF_Interacciones = residuos_contacto(Receptor_Caso,Lig_Caso,receptor_points,DF_Lig,DF_Interacciones,threshold_PH)

    ### Caso Ligando Dador - Receptor Aceptor

    Receptor_Caso ='Aceptor'
    Lig_Caso = 'donnor'

    DF_Interacciones = residuos_contacto(Receptor_Caso,Lig_Caso,receptor_points,DF_Lig,DF_Interacciones,threshold_PH)
    


    ## Aromaticas ##

    aromatic_lig_df = DF_Lig.query('Caso.str.contains("aromatic")', engine='python')
    Receptor_Caso = 'aromatic'
    Lig_Caso = 'aromatic'
    Casos = (aromatic_lig_df['Caso'].unique())

    Sub_Set_Receptor = receptor_points.query('Type == @Receptor_Caso')

    for cas in Casos:
        Sub_Set_Ligando = aromatic_lig_df.query('Caso == @cas')
        center_of_mass = np.mean(np.array(Sub_Set_Ligando.iloc[:,[1,2,3]]),axis=0)
        Matriz_receptor = (np.array(Sub_Set_Receptor.iloc[:,[4,5,6]]))
        distances = np.linalg.norm(Matriz_receptor - center_of_mass, axis=1)
        
        # Filtrar las distancias que son menores a 4.5
        within_distance_indices = np.where(distances < 15)[0]
        
        for idx in within_distance_indices:
            closest_data = Sub_Set_Receptor.iloc[idx]
            min_distance = distances[idx]
            
            # Agregar la información al DataFrame
            DF_Interacciones.loc[len(DF_Interacciones.index)] = [
                closest_data.iloc[1],  # X
                closest_data.iloc[2],  # Y
                closest_data.iloc[3],  # Z
                min_distance,           # Distancia
                Sub_Set_Ligando.iloc[0, 4],  # Identificador del ligando
                Lig_Caso,               # Caso del ligando
                0, 0                    # Valores adicionales
            ]

    DF_Lig_All = generate_df_ligand(pdb_coords)

    DF_Interacciones = DF_Interacciones.drop_duplicates()

   
    ### Validacion de las interacciones ###
    for j in range(0,DF_Interacciones.shape[0]):
        if DF_Interacciones.iloc[j,5] == 'acceptor':
            Aceptor_Antecedent = Busqueda_Antecesor_Lig(DF_Interacciones.iloc[j,4],DF_Lig_All) # Ligando
            Aceptor =  DF_Lig[DF_Lig.iloc[:,0] == DF_Interacciones.iloc[j,4]]
            Aceptor = np.array(Aceptor.iloc[0,[1,2,3]])
            # Filtrar el DataFrame
            resultado = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0]) & (DF_Active_Site['Atom'] == DF_Interacciones.iloc[j,2])]
            Donor = (np.array(resultado[['X' , 'Y' , 'Z']])).reshape(-1)
            DF_Interacciones.iloc[j,6] = angle_three_points(Donor,Aceptor,Aceptor_Antecedent)
        if DF_Interacciones.iloc[j,5] == 'donnor':
            # Ligando
            Donor =  DF_Lig[DF_Lig.iloc[:,0] == DF_Interacciones.iloc[j,4]] 
            Donor = np.array(Donor.iloc[0,[1,2,3]])
            # Receptor
            resultado = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0]) & (DF_Active_Site['Atom'] == DF_Interacciones.iloc[j,2])]
            Aceptor = (np.array(resultado[['X' , 'Y' , 'Z']])).reshape(-1)
            # Ante Receptor
            try:
                atomo_ante = Aceptot_antecedent[DF_Interacciones.iloc[j,1]]
                Atomo = (atomo_ante[DF_Interacciones.iloc[j,2]])
                resultado = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0]) & (DF_Active_Site['Atom'] == Atomo)]
                Aceptor_Antecedent = (np.array(resultado[['X' , 'Y' , 'Z']])).reshape(-1)
            except KeyError:
                resultado = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0]) & (DF_Active_Site['Atom'] == 'C')]
                Aceptor_Antecedent = (np.array(resultado[['X' , 'Y' , 'Z']])).reshape(-1)
            DF_Interacciones.iloc[j,6] = angle_three_points(Donor,Aceptor,Aceptor_Antecedent)
        if DF_Interacciones.iloc[j,5] == 'aromatic':
            Anillo_Proteina = DF_Active_Site[(DF_Active_Site['Pos'] == DF_Interacciones.iloc[j,0])]
            Anillo_Lig = DF_Lig[(DF_Lig['Caso'] == DF_Interacciones.iloc[j,4])]
            DF_Interacciones.iloc[j,6] = Interaccion_Aromatica(Anillo_Proteina,Anillo_Lig)
    
    
    for k in range(0,DF_Interacciones.shape[0]):
        if DF_Interacciones.iloc[k,5] == 'acceptor':
            if (float(DF_Interacciones.iloc[k,3]) < Distances_Hidrogen_Bonds):
                if (float(DF_Interacciones.iloc[k,6]) > 100) and ((float(DF_Interacciones.iloc[k,6]) < 200)): ### Angulos Aceptor
                    DF_Interacciones.iloc[k,7] = 'Yes'
                else:
                    DF_Interacciones.iloc[k,7] = 'No'
            else:
                DF_Interacciones.iloc[k,7] = 'No'
        if DF_Interacciones.iloc[k,5] == 'donnor':
            if (float(DF_Interacciones.iloc[k,3]) < Distances_Hidrogen_Bonds):
                if (float(DF_Interacciones.iloc[k,6]) > 100) and ((float(DF_Interacciones.iloc[k,6]) < 200)): ### Angulos Dador
                    DF_Interacciones.iloc[k,7] = 'Yes'
                else:
                    DF_Interacciones.iloc[k,7] = 'No'
            else:
                DF_Interacciones.iloc[k,7] = 'No'
        if DF_Interacciones.iloc[k,5] == 'aromatic':
            if (float(DF_Interacciones.iloc[k,3]) < Distances_Aromatic):
                if (float(DF_Interacciones.iloc[k,6]) > 0) and (float(DF_Interacciones.iloc[k,6]) < 30):
                    DF_Interacciones.iloc[k,7] = 'Yes'
                elif (float(DF_Interacciones.iloc[k,6]) > 85) and (float(DF_Interacciones.iloc[k,6]) < 95):
                    DF_Interacciones.iloc[k,7] = 'Yes'
                else:
                    DF_Interacciones.iloc[k,7] = 'Yes' ### Hasta ver bien los angulos 
            else:
                DF_Interacciones.iloc[k,7] = 'No'

    DF_Interacciones = DF_Interacciones.drop_duplicates()
    

    ## Calculos Final ##
    dat_final = pd.DataFrame(columns=['Receptor' , 'Ligand' , 'Total' , 'Distance','Ac' ,'Do', 'Ar' , '#Inter','acceptor' ,'donnor', 'aromatic'])
    dat_list = []
    receptor = receptor_pdb.split('.')[0]
    ligand = Ligand_imput.split('.')[0]
    dat_list.append(receptor)
    dat_list.append(ligand)
    DF_Interacciones.to_csv(f'{folder}/Interaction_{receptor}_{ligand}_all.csv')
    dat_list.append(DF_Interacciones.shape[0])
    ## Filtro Distancia##
    DF_Interacciones = DF_Interacciones[DF_Interacciones['Dist'] < threshold_inter_aro]
    DF_Interacciones.to_csv(f'{folder}/Interaction_{receptor}_{ligand}_threshold.csv')

    dat_list.append(DF_Interacciones.shape[0])
    Valores = (dict(DF_Interacciones['Type'].value_counts()))
    try:
        dat_list.append(Valores['acceptor']) # PH
    except KeyError:
        dat_list.append(0) # PH
    try:
        dat_list.append(Valores['donnor']) # PH
    except KeyError:
        dat_list.append(0) # PH        
    try:
        dat_list.append(Valores['aromatic']) # PH
    except KeyError:
        dat_list.append(0) # PH        
    ## Filtro True##
    if not DF_Interacciones.empty:
        DF_Interacciones = DF_Interacciones[DF_Interacciones['Interaction'] == 'Yes']
        DF_Interacciones.to_csv(f'{folder}/Interaction_{receptor}_{ligand}_true.csv')
        dat_list.append(DF_Interacciones.shape[0])
    else:
        dat_list.append(0)
    Valores = (dict(DF_Interacciones['Type'].value_counts()))
    try:
        dat_list.append(Valores['acceptor']) # PH
    except KeyError:
        dat_list.append(0) # PH
    try:
        dat_list.append(Valores['donnor']) # PH
    except KeyError:
        dat_list.append(0) # PH        
    try:
        dat_list.append(Valores['aromatic']) # PH
    except KeyError:
        dat_list.append(0) # PH        
    
    ## VMD_Ploting ##

    if vmd_output == 'Yes':
        scripting_vmd(DF_Interacciones,receptor_points,aromatic_lig_df,DF_Lig,receptor_pdb,chain_receptor,Ligand_imput)
    shutil.copy(Ligand_imput, folder+'/'+Ligand_imput)
    shutil.copy(receptor_pdb, folder+'/'+receptor_pdb)

    ## Resumen data ##
    out_put_file = 'Interactions_close.csv'
    fila = pd.DataFrame([dat_list], columns=dat_final.columns)
    dat_final = pd.concat([dat_final, fila], ignore_index=True)
    dat_final = dat_final.drop_duplicates()
    if not os.path.isfile(out_put_file):
        # Escribir el DataFrame con encabezado si el archivo no existe
         dat_final.to_csv(out_put_file, mode='w', header=True, index=False)
    else:
        # Appendear los datos al archivo existente sin el encabezado
         dat_final.to_csv(out_put_file,  mode='a', header=False, index=False)

    ## Resumen data ##
    out_put_file = 'Interactions_all_count.csv'
    subset_df = dat_final[['Receptor' , 'Ligand' , 'Total', '#Inter','acceptor' ,'donnor', 'aromatic']]
    if not os.path.isfile(out_put_file):
        # Escribir el DataFrame con encabezado si el archivo no existe
         subset_df.to_csv(out_put_file, mode='w', header=True, index=False)
    else:
        # Appendear los datos al archivo existente sin el encabezado
         subset_df.to_csv(out_put_file,  mode='a', header=False, index=False)

    
    
    