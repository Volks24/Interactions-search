"""Hot-points químicos del ligando: aceptores/donores de puentes de H
(RDKit + SMARTS) y anillos aromáticos (geometría + planaridad), más los PNGs
de visualización 2D correspondientes."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
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
