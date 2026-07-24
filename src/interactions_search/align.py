"""align.py: Dado 2 PDB, alinea "protein" sobre "reference"."""

from __future__ import annotations

import argparse

from Bio import PDB
from Bio.PDB import Superimposer

__all__ = ["calcular_alineamiento_rmsd", "main"]


def calcular_alineamiento_rmsd(archivo_referencia, archivo_proteina, cadena_id='A', cadenas_a_mantener=None):
    # Cargar las estructuras desde los archivos PDB
    estructura1 = PDB.PDBParser().get_structure("estructura1", archivo_referencia)
    estructura2 = PDB.PDBParser().get_structure("estructura2", archivo_proteina)

    # Obtener el modelo y la cadena (ajustar según tus necesidades)
    for estructura, archivo in ((estructura1, archivo_referencia), (estructura2, archivo_proteina)):
        cadenas_disponibles = [c.id for c in estructura[0]]
        if cadena_id not in cadenas_disponibles:
            raise ValueError(
                f"La cadena '{cadena_id}' no existe en {archivo}. "
                f"Cadenas disponibles: {cadenas_disponibles}"
            )
    modelo1 = estructura1[0][cadena_id]
    modelo2 = estructura2[0][cadena_id]

    # Eliminar cadenas no deseadas de la proteína a alinear (se conserva la
    # cadena usada para el fit aunque no esté en la lista)
    if cadenas_a_mantener is not None:
        cadenas_a_mantener = set(cadenas_a_mantener) | {cadena_id}
        for chain in list(estructura2[0]):
            if chain.id not in cadenas_a_mantener:
                estructura2[0].detach_child(chain.id)

    # Obtener los residuos de cada cadena
    residuos_A = list(modelo1)
    residuos_B = list(modelo2)

    # Obtener el número del primer y último aminoácido
    primer_aminoacido_A = residuos_A[0].id[1]
    ultimo_aminoacido_A = residuos_A[-1].id[1]
    primer_aminoacido_B = residuos_B[0].id[1]
    ultimo_aminoacido_B = residuos_B[-1].id[1]

    # Calcular el inicio y fin del tramo
    inicio_tramo = max(primer_aminoacido_A, primer_aminoacido_B)
    fin_tramo = min(ultimo_aminoacido_A, ultimo_aminoacido_B)

    # Seleccionar el átomo CA de cada residuo del tramo, por número de residuo
    # (no por índice de lista: puede haber huecos o residuos faltantes)
    seleccion1 = [res['CA'] for res in residuos_A
                  if inicio_tramo <= res.id[1] <= fin_tramo and 'CA' in res]
    seleccion2 = [res['CA'] for res in residuos_B
                  if inicio_tramo <= res.id[1] <= fin_tramo and 'CA' in res]

    if len(seleccion1) != len(seleccion2):
        raise ValueError(
            f"Número de residuos con CA distinto en el tramo {inicio_tramo}-{fin_tramo}: "
            f"{archivo_referencia} tiene {len(seleccion1)}, {archivo_proteina} tiene {len(seleccion2)}."
        )
    if not seleccion1:
        raise ValueError(
            f"No se encontraron átomos CA en el tramo {inicio_tramo}-{fin_tramo} para la cadena '{cadena_id}'."
        )

    # Crear el objeto Superimposer y calcular la transformación de alineamiento
    superimposer = Superimposer()
    superimposer.set_atoms(seleccion1, seleccion2)
    # Aplicar la transformación a toda la estructura (todas las cadenas
    # restantes), no solo a la cadena usada para el fit
    superimposer.apply(list(estructura2.get_atoms()))

    # Calcular el RMSD
    rmsd = superimposer.rms

    # Guardar la estructura alineada en un nuevo archivo PDB
    io = PDB.PDBIO()
    io.set_structure(estructura2)
    output_pdb = f"{archivo_proteina.split('.')[0]}_alig.pdb"
    io.save(output_pdb)

    return inicio_tramo, fin_tramo, rmsd


def main():
    parser = argparse.ArgumentParser(description='Realizar alineamiento estructural y calcular RMSD.')
    parser.add_argument('-R', '--reference', required=True, help='Archivo PDB de referencia')
    parser.add_argument('-P', '--protein', required=True, help='Archivo PDB de la proteína a alinear')
    parser.add_argument('-C', '--Chain', default='A', help='Identificador de cadena (por defecto: A)')
    parser.add_argument('-K', '--keep-chains', nargs='+', default=None, dest='keep_chains',
                        help='Cadenas a conservar en el PDB de salida (elimina el resto). '
                             'La cadena usada para el alineamiento (-C) siempre se conserva.')

    args = parser.parse_args()

    inicio_tramo, fin_tramo, rmsd = calcular_alineamiento_rmsd(
        args.reference, args.protein, args.Chain, args.keep_chains
    )

    # Guardar resultados en un archivo de texto
    output_file = f"{args.protein.split('.')[0]}_resultados.txt"
    with open(output_file, 'w') as f:
        f.write(f"Alineamiento entre: {args.reference.split('.')[0]} y {args.protein.split('.')[0]}\n")
        f.write(f"Inicio Alineamiento: {inicio_tramo}\n")
        f.write(f"Fin Alineamiento: {fin_tramo}\n")
        f.write(f"RMSD: {rmsd:.4f} Å\n")

    print(f"Resultados guardados en {output_file}")


if __name__ == "__main__":
    main()
