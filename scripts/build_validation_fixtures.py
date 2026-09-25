"""Extract fixed crystallographic pockets from locally downloaded RCSB files.

No production detector is imported. This script never creates expected contacts.
Source files: 1IEP.pdb, 1STP.pdb, STI.cif and BTN.cif in --source-dir.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

CASES = [('1IEP', 'STI', 201), ('1STP', 'BTN', 300)]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xyz(line):
    return tuple(float(line[i:i+8]) for i in (30, 38, 46))


def build(source_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {'retrieved': '2026-09-25', 'selection_radius_A': 8.0, 'cases': []}
    for pdb_id, ligand_id, ligand_resid in CASES:
        source = source_dir / f'{pdb_id}.pdb'
        ccd_path = source_dir / f'{ligand_id}.cif'
        lines = source.read_text().splitlines(keepends=True)
        # Both selected structures have a single model. Preserve author atom
        # names, serials, coordinates and residue IDs; retain blank/A altlocs.
        selected = [line for line in lines if line.startswith(('ATOM  ', 'HETATM'))
                    and line[21] == 'A' and line[16] in (' ', 'A')]
        ligand = [line for line in selected if line[17:20] == ligand_id
                  and int(line[22:26]) == ligand_resid]
        assert ligand, pdb_id
        ligand_xyz = [xyz(line) for line in ligand]
        protein = [line for line in selected if line.startswith('ATOM  ')]
        residues = {line[22:27] for line in protein
                    if any(math.dist(xyz(line), point) <= 8.0 for point in ligand_xyz)}
        protein = [line for line in protein if line[22:27] in residues]
        serials = {int(line[6:11]) for line in ligand}
        conect = []
        for line in lines:
            if not line.startswith('CONECT'):
                continue
            ids = [int(line[i:i+5]) for i in range(6, len(line.rstrip()), 5)
                   if line[i:i+5].strip()]
            if ids and ids[0] in serials:
                ids = [ids[0], *[i for i in ids[1:] if i in serials]]
                if len(ids) > 1:
                    conect.append('CONECT' + ''.join(f'{i:5d}' for i in ids) + '\n')
        ccd = MMCIF2Dict(str(ccd_path))
        smiles = next(value for kind, program, value in zip(
            ccd['_pdbx_chem_comp_descriptor.type'], ccd['_pdbx_chem_comp_descriptor.program'],
            ccd['_pdbx_chem_comp_descriptor.descriptor'])
            if kind == 'SMILES_CANONICAL' and program == 'CACTVS')
        prefix = pdb_id.lower()
        outputs = {
            f'{prefix}_receptor.pdb': ''.join(protein) + 'END\n',
            f'{prefix}_ligand.pdb': ''.join(ligand + conect) + 'END\n',
            f'{ligand_id}.smi': smiles + '\n',
        }
        for name, content in outputs.items():
            (output_dir / name).write_text(content)
        manifest['cases'].append({
            'pdb_id': pdb_id, 'chain': 'A', 'ligand': ligand_id, 'ligand_resid': ligand_resid,
            'source_url': f'https://files.rcsb.org/download/{pdb_id}.pdb',
            'ccd_url': f'https://files.rcsb.org/ligands/download/{ligand_id}.cif',
            'source_sha256': digest(source), 'ccd_sha256': digest(ccd_path),
            'smiles_source': 'CCD SMILES_CANONICAL / CACTVS; deposited neutral component',
            'receptor_residues': sorted(int(residue[:4]) for residue in residues),
            'receptor_atoms': len(protein), 'ligand_atoms': len(ligand),
            'files': {name: digest(output_dir / name) for name in outputs},
        })
    (output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=Path('tests/fixtures/real'))
    args = parser.parse_args()
    build(args.source_dir, args.output_dir)
