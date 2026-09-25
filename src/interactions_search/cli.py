"""Interfaz de línea de comandos: resuelve los pares (receptor, ligando) a
analizar (uno, batch, o split automático de un PDB complejo) y llama a
analyze_pair() por cada uno."""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

from interactions_search.config import load_config
from interactions_search.hotspot_pocket import analyze_hotspot_pockets
from interactions_search.io_pdb import split_pdb, validate_inputs
from interactions_search.ligand_chemistry import load_reference
from interactions_search.pipeline import (
    analyze_pair,
    analyze_probe,
    analyze_site_bias,
)
from interactions_search.probe import PROBE_TYPES, read_probe_file

__all__ = ["main"]


def _build_cfg(config_path):
    """Dict de configuración que consumen analyze_pair()/analyze_probe()."""
    config = load_config(config_path)
    return {
        **config.distancias.model_dump(),
        **config.angulos.model_dump(),
        'ligand_plot': config.options.ligand_plot,
        'vmd_output': config.options.vmd_output,
        'cumulative_output': config.options.cumulative_output,
        'Interaction_Coord_Source': config.options.interaction_coord,
        'Volume_Plot': config.options.volume_plot,
        'Bias': config.options.bias,
        'Bias_Validated_Only': config.options.bias_validated_only,
        'Distancia_Hidrofobica': config.distancias.Distances_Hidrofobica,
        'Distancia_Centro_Activo': config.distancias.centroid_distance,
        'Ring_Planarity_RMSD_Max': config.aromaticidad.Ring_Planarity_RMSD_Max,
        'Pocket_Min_Residues': config.pockets.min_residues,
        'Pocket_Coverage_Threshold': config.pockets.coverage_threshold,
        'Pocket_Density_Radius': config.pockets.density_radius,
        'Aceptores_Prot': config.acceptors,
        'Dadores_Prot': config.donors,
        'Aceptot_antecedent': config.acceptors_antecedent,
        'Special_case': config.special,
    }


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
            '  Site bias    : -r proteina.pdb -c A --site-point X Y Z [--site-radius 10]\n'
            '                 [--site-method atom|ideal]\n'
            '                 (sin -l: puntos de bias del receptor alrededor de una\n'
            '                 coordenada arbitraria, sin necesidad de una pose de ligando.\n'
            '                 "ideal" requiere el receptor protonado)\n'
            '  Sondeo       : -r proteina.pdb -c A --probe X Y Z [--probe X Y Z ...]\n'
            '                 -r proteina.pdb -c A --probe-file puntos.bpf|.pdb|.csv\n'
            '                 [--probe-type all|acceptor|donor|aromatic|hydrophobic|cation|anion]\n'
            '                 (sin -l: simula qué interacciones haría un átomo/grupo del\n'
            '                 ligando ubicado en cada coordenada)\n'
            '  Hotspots MD  : -r proteina.pdb -c A --hotspots results_global/ [--exclude-res STI]\n'
            '                 (pockets + grilla anotada a partir de clusters de hotspots\n'
            '                 aceptor/donor/hidrofóbico de una dinámica con cosolvente)\n'
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
    parser.add_argument('--site-point', nargs=3, type=float, default=None,
                        metavar=('X', 'Y', 'Z'),
                        help='Modo especial (requiere -r, no -l ni -x): en vez de analizar un '
                             'ligando, busca los residuos del receptor dentro de --site-radius '
                             'de esta coordenada y exporta sus puntos de bias '
                             '(aceptor/donor/aromático) como .bpf + PDB dummy.')
    parser.add_argument('--site-radius', type=float, default=10.0,
                        help='Radio (Å) de búsqueda alrededor de --site-point (default: 10.0).')
    parser.add_argument('--site-method', choices=['atom', 'ideal'], default='atom',
                        help="'atom' (default): un punto de bias por átomo/anillo del receptor, "
                             "en su propia coordenada. 'ideal': abanico de puntos ideales de "
                             "H-bond/stacking (distancia+ángulo+diedro) hacia donde debería caer "
                             "el átomo complementario del ligando; requiere el receptor protonado "
                             "(nomenclatura Maestro/Amber) para los grupos donores.")
    parser.add_argument('--probe', nargs=3, type=float, action='append', default=None,
                        metavar=('X', 'Y', 'Z'),
                        help='Modo sondeo (requiere -r, no -l ni -x): simula las interacciones '
                             'que haría un átomo/grupo del ligando en esta coordenada. '
                             'Repetible para varios puntos.')
    parser.add_argument('--probe-file', default=None, metavar='POINTS',
                        help='Puntos de sondeo desde archivo: .bpf (ej. salida de '
                             '--site-point), .pdb (un punto por átomo; resname DON/ACC/ARO '
                             'fija el tipo) o texto/CSV "x y z [tipo]".')
    parser.add_argument('--probe-type', nargs='+', choices=['all', *PROBE_TYPES],
                        default=['all'],
                        help='Rol(es) con que se sondea cada punto sin tipo propio '
                             '(default: all).')
    parser.add_argument('--hotspots', default=None, metavar='DIR',
                        help='Modo hotspots (requiere -r): directorio con acceptors/, donors/ '
                             'y/o hydrophobics/ (clusters.csv, cluster_points.pdb, *.dx). '
                             'Agrupa los hotspots en sitios y arma el pocket + grilla de cada uno '
                             '(parámetros en la sección hotspot_pocket del YAML).')
    parser.add_argument('--exclude-res', nargs='+', default=[], metavar='RESNAME',
                        help='Resnames a excluir del receptor en el modo --hotspots (ej. '
                             'ligandos guardados como ATOM: --exclude-res STI MS7).')

    reference_group = parser.add_mutually_exclusive_group()
    reference_group.add_argument('--ligand-smiles', help='SMILES de referencia química opcional.')
    reference_group.add_argument('--ligand-sdf', help='SDF de referencia con una sola molécula.')
    parser.add_argument('--legacy-rings', action='store_true',
                        help='Conservar el filtro histórico de anillos de más de cinco átomos.')
    args = parser.parse_args()
    has_reference = args.ligand_smiles is not None or args.ligand_sdf is not None
    if (has_reference or args.legacy_rings) and (
            args.hotspots or args.probe or args.probe_file or args.site_point is not None):
        parser.error('Las opciones de química del ligando requieren un análisis con ligando.')
    if has_reference and args.ligand_input and len(args.ligand_input) != 1:
        parser.error('La referencia química requiere un solo ligando por ejecución.')
    try:
        ligand_reference = load_reference(smiles=args.ligand_smiles, sdf=args.ligand_sdf)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))

    run_context = {
        'cli_arguments': vars(args),
        'source_files': {key: value for key, value in {
            'source_complex': args.complex_pdb, 'source_config': args.config,
            'source_ligand_sdf': args.ligand_sdf, 'source_probe_file': args.probe_file,
        }.items() if value is not None},
    }

    # ── Modo especial: pockets a partir de hotspots de dinámica ────
    if args.hotspots:
        if (args.complex_pdb or args.ligand_input or args.site_point is not None
                or args.probe or args.probe_file):
            parser.error("--hotspots requiere -r y no es compatible con -x, -l, "
                         "--site-point ni --probe.")
        analyze_hotspot_pockets(args.receptor_pdb, args.chain_receptor, args.hotspots,
                                _build_cfg(args.config), load_config(args.config).hotspot_pocket,
                                exclude_resnames=args.exclude_res, run_context=run_context)
        return

    # ── Modo especial: sondeo de coordenadas (sin ligando) ─────────
    if args.probe or args.probe_file:
        if args.complex_pdb or args.ligand_input or args.site_point is not None:
            parser.error("--probe/--probe-file requiere -r y no es compatible con "
                         "-x, -l ni --site-point.")
        default_types = list(PROBE_TYPES) if 'all' in args.probe_type else args.probe_type
        points = [(x, y, z, None) for x, y, z in (args.probe or [])]
        if args.probe_file:
            file_points = read_probe_file(args.probe_file)
            if not file_points:
                parser.error(f"No se encontraron puntos en {args.probe_file}.")
            points += file_points
        points = [(x, y, z, [t] if t else default_types) for x, y, z, t in points]

        if args.probe_file:
            label = Path(args.probe_file).stem
        elif len(points) == 1:
            x, y, z, _ = points[0]
            label = f'{x:.1f}_{y:.1f}_{z:.1f}'
        else:
            label = f'{len(points)}pts'
        analyze_probe(args.receptor_pdb, args.chain_receptor, points, _build_cfg(args.config),
                      label=label, run_context=run_context)
        return

    # ── Modo especial: site bias (sin ligando) ─────────────────────
    if args.site_point is not None:
        if args.complex_pdb:
            parser.error("--site-point no es compatible con -x/--complex; usá -r/--receptor_pdb.")
        if not args.receptor_pdb:
            parser.error("--site-point requiere -r/--receptor_pdb.")
        analyze_site_bias(args.receptor_pdb, args.chain_receptor, tuple(args.site_point),
                          args.site_radius, _build_cfg(args.config),
                          method=args.site_method, run_context=run_context)
        return

    # ── Resolver lista de pares (receptor, ligando) ───────────────
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
            lig_key = args.lig_name
            if lig_key not in het_paths:
                # Nombre "pelado" (ej. NAI) con copias en varias cadenas (NAI_B, NAI_C, ...):
                # se resuelve automáticamente a la copia de la cadena pasada en -c.
                by_chain_key = f'{lig_key}_{args.chain_receptor}'
                if by_chain_key in het_paths:
                    lig_key = by_chain_key
                else:
                    print(f"  Error: '{args.lig_name}' not found. Available: {list(het_paths.keys())}")
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                    sys.exit(1)
            pairs = [(str(protein_path), str(het_paths[lig_key]))]
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
    cfg = _build_cfg(args.config)

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
        try:
            analyze_pair(receptor_pdb, Ligand_imput, args.chain_receptor, cfg,
                         ligand_reference=ligand_reference, legacy_rings=args.legacy_rings,
                         run_context=run_context)
        except ValueError as exc:
            if tmp_dir:
                shutil.rmtree(tmp_dir, ignore_errors=True)
            parser.error(str(exc))
        n_ok += 1

    # ── Limpiar directorio temporal del split ─────────────────────
    if tmp_dir:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"\nAnalysis complete: {n_ok} pair(s) processed, {n_skip} skipped.")


if __name__ == '__main__':
    main()
