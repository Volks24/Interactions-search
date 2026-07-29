"""Interfaz de línea de comandos: resuelve los pares (receptor, ligando) a
analizar (uno, batch, o split automático de un PDB complejo) y llama a
analyze_pair() por cada uno."""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile

from interactions_search.io_pdb import split_pdb, validate_inputs
from interactions_search.pipeline import analyze_pair, carga_variables

__all__ = ["main"]


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
