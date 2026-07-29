"""Generación de scripts .tcl para visualización en VMD: escena base
(H-bond/aromático), contactos hidrofóbicos, superficies de pockets, y una
escena combinada con todo junto."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from interactions_search.pockets import _POCKET_COLORIDS

__all__ = [
    "scripting_vmd",
    "scripting_vmd_hydrophobic",
    "scripting_vmd_pockets",
    "scripting_vmd_combined",
]


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
