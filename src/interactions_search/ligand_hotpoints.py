"""Hot-points químicos del ligando: aceptores/donores de puentes de H
(RDKit + SMARTS) y anillos aromáticos (geometría + planaridad), más los PNGs
de visualización 2D correspondientes."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from interactions_search.geometry import _ring_planarity_rmsd

__all__ = [
    "search_hot_points",
    "generate_df_ligand",
    "search_rings",
    "visualize_rings",
]


def _openbabel_h_flags(pdb_path, n_atoms):
    """Por cada atomo (mismo orden/indice 0-based que RDKit al leer el mismo
    PDB), True si OpenBabel -- que SI percibe correctamente enlaces dobles
    C=O/C=N a partir de geometria 3D sin CONECT records -- le asigna al menos
    un H (implicito o explicito).

    Por que hace falta esto: Chem.MolFromPDBFile no percibe ordenes de enlace
    de forma fiable cuando el PDB no trae CONECT (el caso normal para
    ligandos de cristalografia); en la practica termina asignando enlace
    simple a la mayoria de los oxigenos terminales (carbonilo, sulfonilo,
    carboxilato) y rellena la valencia con un H implicito, marcando como
    "donor" atomos que en realidad no tienen ningun H real. OpenBabel
    (`PerceiveBondOrders`) resuelve esto correctamente para los mismos casos
    (ya se usa con el mismo proposito, para deteccion de anillos aromaticos,
    en el companion script prepare_bias.py de AutoDock Bias).

    Devuelve una lista de bool de largo n_atoms; si no se puede leer el PDB
    con OpenBabel, devuelve None (el llamador debe tratarlo como "no
    disponible" y no filtrar nada, para no romper el comportamiento si
    OpenBabel no esta instalado)."""
    try:
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("pdb", "pdb")
        obmol = ob.OBMol()
        if not conv.ReadFile(obmol, str(pdb_path)):
            return None
        obmol.PerceiveBondOrders()
    except Exception:
        return None

    flags = [None] * n_atoms
    for atom in ob.OBMolAtomIter(obmol):
        idx = atom.GetIdx() - 1  # OpenBabel es 1-indexed, RDKit 0-indexed
        if 0 <= idx < n_atoms:
            flags[idx] = (atom.GetImplicitHCount() + atom.ExplicitHydrogenCount()) > 0
    if any(f is None for f in flags):
        # desalineacion de indices entre RDKit y OpenBabel (atomos distintos,
        # PDB no estandar, etc.) -- mas seguro no filtrar que filtrar mal
        return None
    return flags


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

    # [N;H1] se saco de acceptor_smarts: en RDKit [N;H] (donor_smarts) es
    # exactamente [N;H1], asi que un N con un solo H (amina secundaria,
    # N-H de amida/sulfonamida) matcheaba los dos patrones a la vez.
    # Quimicamente ese N dona via su H (el par libre suele estar
    # deslocalizado hacia el grupo vecino, ej. sulfonamida/amida) -- no tiene
    # sentido tratarlo tambien como aceptor. Un N sin H (piridina, amina
    # terciaria, amida sin H) sigue siendo aceptor via [N;H0].
    acceptor_smarts = ['[O;H1]', '[O;H0]', '[N;H0]', '[n]', '[o]', '[N+]']
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

    # Validacion cruzada con OpenBabel: RDKit (Chem.MolFromPDBFile) no
    # percibe bien los ordenes de enlace sin CONECT records y suele marcar
    # oxigenos/nitrogenos terminales como si tuvieran H aunque en realidad
    # sean carbonilo/sulfonilo/carboxilato (ver docstring de
    # _openbabel_h_flags). Se descartan los "donor" que OpenBabel dice que
    # no tienen ningun H real.
    ob_has_h = _openbabel_h_flags(Ligand_imput, mol.GetNumAtoms())
    if ob_has_h is not None:
        donor_atoms = [idx for idx in donor_atoms if ob_has_h[idx]]

    if ligand_plot == 'Yes':
        stem = Path(Ligand_imput).stem
        acc_labels = {idx: pdb_coords[idx][1] for idx in acceptor_atoms if idx < len(pdb_coords)}
        _draw_mol_labeled(mol, acceptor_atoms, acc_labels, f"{folder}/{stem}_acceptors.png")
        don_labels = {idx: pdb_coords[idx][1] for idx in donor_atoms if idx < len(pdb_coords)}
        _draw_mol_labeled(mol, donor_atoms, don_labels, f"{folder}/{stem}_donors.png")

    return acceptor_atoms, donor_atoms


def generate_df_ligand(pdb_coords):

    columns = ['Atom ID', 'Element', 'Residue Name', 'Chain ID', 'Residue Number', 'X', 'Y', 'Z']
    df_ligand = pd.DataFrame(pdb_coords, columns=columns)

    return(df_ligand)


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
