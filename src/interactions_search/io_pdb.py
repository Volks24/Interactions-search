"""Lectura, separación y validación de archivos PDB usados por el pipeline
de interacciones (previo al análisis de hot-points / contactos)."""
from __future__ import annotations

import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np

__all__ = [
    "get_atom_coords",
    "extract_coords_from_pdb",
    "remove_bias",
    "split_pdb",
    "validate_inputs",
]


def get_atom_coords(mol, atom_idx):
    conf = mol.GetConformer()
    pos = conf.GetAtomPosition(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    return f"{atom.GetSymbol()} {atom_idx}: ({pos.x}, {pos.y}, {pos.z})"


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

    Los grupos se separan por (cadena, resname): si un mismo resname aparece en más
    de una cadena (ej. varias copias del complejo en la unidad asimétrica), cada
    copia se escribe en su propio archivo — nunca se fusionan copias de cadenas
    distintas, porque promediar sus coordenadas para el centro de masa del ligando
    arruina la búsqueda de sitio activo (que sí queda acotada a una única cadena).

    Retorna:
        protein_path : Path al PDB de la proteína
        het_paths    : dict {label: Path}  —  vacío si no hay HETATM.
                       label es el resname si es única su cadena, o
                       '{resname}_{cadena}' si el resname aparece en varias cadenas.
    """
    WATER_NAMES = {'HOH', 'WAT', 'TIP', 'TIP3', 'SOL', 'DOD'}
    force_ligand_names = {n.upper() for n in force_ligand_names} if force_ligand_names else set()

    protein_lines = []
    het_lines     = defaultdict(list)   # (chain, resname) -> [líneas HETATM]
    conect_lines  = []
    header_lines  = []

    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == 'ATOM':
                resname = line[17:20].strip()
                if resname in force_ligand_names:
                    het_lines[(line[21], resname)].append('HETATM' + line[6:])
                else:
                    protein_lines.append(line)
            elif rec == 'HETATM':
                resname = line[17:20].strip()
                if exclude_water and resname in WATER_NAMES:
                    continue
                het_lines[(line[21], resname)].append(line)
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
    # Etiqueta cada grupo: 'RESNAME' si su cadena es única, 'RESNAME_CHAIN' si el
    # resname se repite en más de una cadena (evita fusionar copias distintas).
    resname_chains = defaultdict(set)
    for chain, resname in het_lines:
        resname_chains[resname].add(chain)

    def _label(chain, resname):
        return resname if len(resname_chains[resname]) == 1 else f'{resname}_{chain}'

    # Pre-indexar seriales de cada grupo para filtrar CONECT
    group_serials = {}
    for key, lines in het_lines.items():
        serials = set()
        for l in lines:
            try:
                serials.add(int(l[6:11]))
            except ValueError:
                pass
        group_serials[key] = serials

    het_paths = {}
    for (chain, resname), lines in het_lines.items():
        label = _label(chain, resname)
        het_path = out / f'{stem}_{label}.pdb'
        with open(het_path, 'w') as f:
            f.writelines(lines)
            # CONECT cuyos átomos pertenecen a este grupo
            my_serials = group_serials[(chain, resname)]
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
        het_paths[label] = het_path

    return protein_path, het_paths


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
