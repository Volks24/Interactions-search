#!/usr/bin/env python3
"""Export all built hotspot pocket grids as two site-selectable PDB files."""
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

_SURFACE_SAMPLE_SPACING = 0.5
_TYPE_CODE = {'acceptor': 'A', 'donor': 'D', 'hydrophobic': 'H', 'none': 'N'}
_BASE36 = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'


def _read_grid(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def _read_dx(path):
    text = path.read_text()
    header, values_text = text.split('data follows', 1)
    counts = tuple(map(int, re.search(
        r'gridpositions counts\s+(\d+)\s+(\d+)\s+(\d+)', header).groups()))
    origin = np.array(re.search(
        r'origin\s+(\S+)\s+(\S+)\s+(\S+)', header).groups(), dtype=float)
    deltas = [np.array(row, dtype=float) for row in re.findall(
        r'delta\s+(\S+)\s+(\S+)\s+(\S+)', header)]
    spacing = np.array([deltas[0][0], deltas[1][1], deltas[2][2]])
    values_text = values_text.split('attribute', 1)[0]
    values = np.fromstring(values_text, sep=' ').reshape(counts)
    return origin, spacing, values


def _surface_rows(site_dir, rows):
    """Sample marching-cubes isosurface vertices and retain nearest grid metadata."""
    try:
        from skimage.measure import marching_cubes
    except ImportError as error:
        raise RuntimeError(
            'Smooth PDB surface export requires scikit-image (pip install scikit-image).') from error
    origin, spacing, mask = _read_dx(site_dir / 'pocket_mask.dx')
    vertices, _, _, _ = marching_cubes(mask.astype(np.float32), level=0.5,
                                        spacing=tuple(spacing))
    vertices += origin
    vertices = np.unique(np.round(vertices, 3), axis=0)
    samples = {}
    for point in vertices:
        cell = tuple(np.floor(point / _SURFACE_SAMPLE_SPACING).astype(int))
        center = (np.array(cell) + 0.5) * _SURFACE_SAMPLE_SPACING
        distance_sq = float(np.sum((point - center) ** 2))
        if cell not in samples or distance_sq < samples[cell][0]:
            samples[cell] = distance_sq, point
    vertices = np.array([entry[1] for entry in samples.values()])
    grid_xyz = np.array([[float(row[axis]) for axis in ('X', 'Y', 'Z')] for row in rows])
    _, nearest = cKDTree(grid_xyz).query(vertices)
    sampled = []
    for point, row_index in zip(vertices, nearest):
        row = rows[int(row_index)].copy()
        row.update(X=str(point[0]), Y=str(point[1]), Z=str(point[2]))
        sampled.append(row)
    return sampled


def _remarks(sites, style, metadata):
    description = 'all voxel centers' if style == 'spheres' else 'sampled interpolated isosurface'
    lines = [f'REMARK 900 HOTSPOT POCKETS: {description}; SITE SELECTED BY CHAIN',
             'REMARK 900 RESNAME PCK; OCCUPANCY = BURIEDNESS; B FACTOR = BEST_DG',
             'REMARK 900 ATOM NAME FIRST LETTER: A=acceptor D=donor H=hydrophobic N=none']
    for site_index, (site_id, rows) in enumerate(sites):
        lines.append(f'REMARK 901 CHAIN {_CHAIN_IDS[site_index]} = SITE {site_id}; {len(rows)} points')
        info = metadata.get(site_id)
        if info:
            lines.append(
                f"REMARK 902 SITE {site_id}: HOTSPOTS {info['N_Hotspots']}; "
                f"RESIDUES {info['N_Residues']}; VOLUME_A3 {info['Volume_A3']}")
            lines.append(
                f"REMARK 903 SITE {site_id}: DG_SUM {info['DG_Sum']}; "
                f"DG_BEST {info['DG_Best']}; CENTER {info['Center_X']},"
                f"{info['Center_Y']},{info['Center_Z']}")
    return lines


def _write_pdb(path, sites, style, metadata):
    count = 0
    with path.open('w') as stream:
        for line in _remarks(sites, style, metadata):
            stream.write(line + '\n')
        for site_index, (site_id, rows) in enumerate(sites):
            site_chain = _CHAIN_IDS[site_index]
            type_counts = {name: 0 for name in _TYPE_CODE}
            for row in rows:
                count += 1
                best_type = row['Best_Type'] if row['Best_Type'] in _TYPE_CODE else 'none'
                type_counts[best_type] += 1
                suffix = type_counts[best_type] - 1
                encoded_suffix = ''
                while True:
                    encoded_suffix = _BASE36[suffix % 36] + encoded_suffix
                    suffix //= 36
                    if not suffix:
                        break
                atom_name = _TYPE_CODE[best_type] + encoded_suffix.rjust(3, '0')
                x, y, z = (float(row[axis]) for axis in ('X', 'Y', 'Z'))
                occupancy = float(row['Buriedness'])
                beta_value = float(row['Best_DG']) if row['Best_DG'] else 0.0
                beta = beta_value if math.isfinite(beta_value) else 0.0
                stream.write(
                    f'HETATM{count:5d} {atom_name:4s} PCK {site_chain}   1    '
                    f'{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{beta:6.2f}'
                    f'          C\n')
        stream.write('END\n')
    return count


def _combine_with_protein(path, protein_path, chain_id, input_dir, built_sites,
                          surface_path):
    if len(built_sites) > len(_CHAIN_IDS):
        raise ValueError('PDB supports at most 62 one-character site chain IDs.')
    protein_source = Path(protein_path).read_text().splitlines()
    chain_to_site = {_CHAIN_IDS[i]: site_id for i, (site_id, _) in enumerate(built_sites)}
    surface_by_site = {}
    for line in Path(surface_path).read_text().splitlines():
        if line.startswith('HETATM') and line[21:22] in chain_to_site:
            surface_by_site.setdefault(chain_to_site[line[21:22]], []).append(line)

    with path.open('w') as stream:
        stream.write('REMARK 900 EACH CHAIN CONTAINS ONE SITE PROTEIN FRAGMENT AND ITS POCKET\n')
        stream.write('REMARK 900 SELECT A CHAIN TO SHOW ITS PROTEIN FRAGMENT AND POCKET\n')
        stream.write(f'REMARK 900 SOURCE RECEPTOR: {Path(protein_path).name}\n')
        stream.write(f'REMARK 900 SOURCE POCKET DATA: {Path(input_dir).name}\n')
        serial = 0
        protein_atom_count = 0
        summary_path = Path(input_dir) / 'sites_summary.csv'
        summary = {int(row['Site']): row for row in _read_grid(summary_path)} \
            if summary_path.is_file() else {}
        for site_index, (site_id, site_dir) in enumerate(built_sites):
            site_chain = _CHAIN_IDS[site_index]
            residue_ids = {
                (int(row['Pos']), row['Residue'])
                for row in _read_grid(site_dir / 'residues.csv')
            }
            stream.write(f'REMARK 901 CHAIN {site_chain} = SITE {site_id}; PROTEIN RESIDUES\n')
            info = summary.get(site_id)
            if info:
                stream.write(
                    f"REMARK 902 SITE {site_id}: HOTSPOTS {info['N_Hotspots']}; "
                    f"DG_SUM {info['DG_Sum']}; DG_BEST {info['DG_Best']}\n")
            protein_count = 0
            used_resids = set()
            for source_line in protein_source:
                if (source_line.startswith('ATOM  ') and source_line[21:22] == chain_id
                        and (int(source_line[22:26]), source_line[17:20].strip()) in residue_ids):
                    serial += 1
                    line = (source_line[:6] + f'{serial:5d}' + source_line[11:21]
                            + site_chain + source_line[22:])
                    stream.write(line + '\n')
                    protein_count += 1
                    used_resids.add(int(source_line[22:26]))
            if not protein_count:
                raise ValueError(f'No pocket residues from chain {chain_id!r} found for site {site_id}')
            protein_atom_count += protein_count
            stream.write('TER\n')
            grid_resid = next(candidate for candidate in range(1, 10000)
                              if candidate not in used_resids)
            for source_line in surface_by_site.get(site_id, []):
                serial += 1
                line = (source_line[:6] + f'{serial:5d}' + source_line[11:21] + site_chain
                        + f'{grid_resid:4d}' + source_line[26:])
                stream.write(line + '\n')
            stream.write('TER\n')
        stream.write('END\n')
    return protein_atom_count


def _write_vmd_scenes(output_dir):
    scenes = {
        'view_pockets_spheres.tcl': (
            'pockets_global_spheres.pdb',
            'mol modstyle 0 $m VDW 0.18 12\nmol modcolor 0 $m Beta\n'),
        'view_pockets_surface.tcl': (
            'pockets_global_surface.pdb',
            'mol modstyle 0 $m Surf 1.4 0\nmol modcolor 0 $m Beta\n'),
    }
    combined_name = 'pockets_global_surface_with_protein.pdb'
    combined = [
        'mol modselect 0 $m "not resname PCK"',
        'mol modstyle 0 $m Licorice 0.2',
        'mol modcolor 0 $m ColorID 8',
    ]
    for chain in _CHAIN_IDS[:5]:
        rep = _CHAIN_IDS.index(chain) + 1
        combined.extend([
            'mol addrep $m',
            f'mol modselect {rep} $m "chain {chain} and resname PCK"',
            f'mol modstyle {rep} $m Surf 1.4 0',
            f'mol modcolor {rep} $m Beta',
        ])
    scenes['view_pockets_with_protein.tcl'] = (combined_name, '\n'.join(combined) + '\n')

    for filename, (pdb_name, representation) in scenes.items():
        path = output_dir / filename
        path.write_text(
            'set scriptDir [file dirname [info script]]\n'
            'cd $scriptDir\n'
            f'set m [mol new "{pdb_name}" type pdb]\n'
            'mol modselect 0 $m all\n'
            + representation
            + 'display projection orthographic\ndisplay resetview\n')


def export(input_dir, output_dir=None, protein_pdb=None, chain_id='B'):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir) if output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    sites = []
    built_sites = []
    for directory in sorted(input_dir.glob('site_*'), key=lambda path: int(path.name.split('_')[1])):
        grid_path = directory / 'grid.csv'
        if not grid_path.is_file():
            continue
        rows = _read_grid(grid_path)
        if rows:
            site_id = int(directory.name.split('_')[1])
            sites.append((site_id, rows))
            built_sites.append((site_id, directory))
    if not sites:
        raise ValueError(f'No built site grids found in {input_dir}')
    metadata = {}
    summary_path = input_dir / 'sites_summary.csv'
    if summary_path.is_file():
        for row in _read_grid(summary_path):
            try:
                site_id = int(row['Site'])
                if row['Built'] == 'Yes':
                    metadata[site_id] = row
            except (KeyError, ValueError):
                continue

    all_points = [(site_id, rows) for site_id, rows in sites]
    smooth_surface = [(site_id, _surface_rows(site_dir, rows))
                      for (site_id, rows), (_, site_dir) in zip(sites, built_sites)]
    sphere_path = output_dir / 'pockets_global_spheres.pdb'
    surface_path = output_dir / 'pockets_global_surface.pdb'
    sphere_count = _write_pdb(sphere_path, all_points, 'spheres', metadata)
    surface_count = _write_pdb(surface_path, smooth_surface, 'surface', metadata)
    combined_path = None
    protein_atom_count = None
    if protein_pdb:
        combined_path = output_dir / 'pockets_global_surface_with_protein.pdb'
        protein_atom_count = _combine_with_protein(
            combined_path, protein_pdb, chain_id, input_dir, built_sites, surface_path)
    _write_vmd_scenes(output_dir)
    return (sphere_path, sphere_count, surface_path, surface_count,
            combined_path, protein_atom_count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_dir', help='hotspot_pockets output directory')
    parser.add_argument('--output-dir', default=None,
                        help='destination (defaults to input_dir)')
    parser.add_argument('--protein-pdb', default=None,
                        help='include residues used by the pockets from this PDB')
    parser.add_argument('--chain', default='B', help='protein chain (default: B)')
    args = parser.parse_args()
    sphere_path, sphere_count, surface_path, surface_count, combined_path, protein_count = export(
        args.input_dir, args.output_dir, args.protein_pdb, args.chain)
    print(f'Spheres: {sphere_path} ({sphere_count} points)')
    print(f'Surface shell: {surface_path} ({surface_count} sampled surface points)')
    if combined_path:
        print(f'Protein + surface shell: {combined_path} ({protein_count} protein atoms)')
    print('VMD scenes: view_pockets_spheres.tcl, view_pockets_surface.tcl'
          + (', view_pockets_with_protein.tcl' if combined_path else ''))
    print('Select pockets by chain: A=site 1, B=site 2, etc.')


if __name__ == '__main__':
    main()
