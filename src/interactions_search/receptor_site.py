"""Sitio activo del receptor: residuos dentro del radio de búsqueda alrededor
del centro de masa del ligando, y sus puntos de interés (aceptores, donores,
centroides de anillos aromáticos) para la búsqueda de contactos."""
from __future__ import annotations

import math

import pandas as pd

from interactions_search.geometry import center_of_mass, get_aromatic_coord

__all__ = [
    "active_site_residues",
    "Coordenadas_interes_receptor",
]


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
