"""Bias probe file (.bpf) para GOLD: puntos de los hot-points químicos del
ligando (aceptor/donor/aromático), exportados como .bpf y como PDB dummy
para visualización en VMD.

También soporta el caso inverso: puntos de bias del *receptor* (aceptor/
donor/aromático de los residuos del sitio activo, ver receptor_site.py) en
vez de los del ligando -- usado por el modo "site bias" (analyze_site_bias en
pipeline.py) para generar un .bpf alrededor de una coordenada arbitraria sin
necesidad de una pose de ligando."""
from __future__ import annotations

__all__ = ["export_bpf", "export_bpf_pdb", "export_bpf_receptor", "export_bpf_pdb_receptor",
           "export_bpf_points", "export_bpf_pdb_points"]

# Vset/r fijos por tipo, tomados de receptor_fs(6).bpf (referencia de formato
# provista) — no varían por átomo, son la misma convención para todo don/acc/aro.
_BPF_PARAMS = {
    'don': {'Vset': -2.72, 'r': 1.20},
    'acc': {'Vset': -2.28, 'r': 0.80},
    'aro': {'Vset': -2.00, 'r': 2.00},
}


def _collect_bias_points(DF_Lig, DF_true=None):
    """Junta los puntos de bias (un punto por átomo aceptor/donor, un punto
    por anillo aromático en su centroide) en el mismo orden que usan
    export_bpf() y export_bpf_pdb(), para que ambos archivos queden
    alineados fila a fila.

    Si DF_true es None (default): usa TODOS los hot-points químicos del
    ligando (search_hot_points/search_rings), sin importar si esa posición
    llegó a formar una interacción real con el receptor en esta pose.

    Si se pasa DF_true (interacciones validadas, Interaction == 'Yes'): se
    restringe a los hot-points que sí participan en una interacción validada
    — átomos aceptor/donor por serial (LigID) y anillos por su etiqueta
    ('aromatic' o 'pi_cation', ambos indican que el anillo hace contacto
    real). 'Evidence-based' en vez de 'todo lo químicamente posible'."""
    donors    = DF_Lig[DF_Lig['Caso'] == 'donor']
    acceptors = DF_Lig[DF_Lig['Caso'] == 'acceptor']
    aromatic  = DF_Lig[DF_Lig['Caso'].astype(str).str.startswith('aromatic')]

    if DF_true is not None:
        validated_ligids = set(
            DF_true.loc[DF_true['Type'].isin(['acceptor', 'donor']), 'LigID'].dropna().astype(int))
        donors    = donors[donors['Atom ID'].astype(int).isin(validated_ligids)]
        acceptors = acceptors[acceptors['Atom ID'].astype(int).isin(validated_ligids)]
        validated_rings = set(DF_true.loc[DF_true['Type'].isin(['aromatic', 'pi_cation']), 'Lig'])
        aromatic  = aromatic[aromatic['Caso'].isin(validated_rings)]

    rows = []
    for _, r in donors.iterrows():
        rows.append((float(r['Coord X']), float(r['Coord Y']), float(r['Coord Z']), 'don'))
    for _, r in acceptors.iterrows():
        rows.append((float(r['Coord X']), float(r['Coord Y']), float(r['Coord Z']), 'acc'))
    for _, group in aromatic.groupby('Caso'):
        center = group[['Coord X', 'Coord Y', 'Coord Z']].astype(float).mean(axis=0)
        rows.append((float(center['Coord X']), float(center['Coord Y']), float(center['Coord Z']), 'aro'))
    return rows


def _write_bpf(rows, filepath):
    with open(filepath, 'w') as f:
        f.write('x\ty\tz\tVset\tr\ttype\n')
        for x, y, z, tipo in rows:
            p = _BPF_PARAMS[tipo]
            f.write(f"{x:.3f}\t{y:.3f}\t{z:.3f}\t{p['Vset']}\t{p['r']}\t{tipo}\n")


_BPF_RESNAME = {'don': 'DON', 'acc': 'ACC', 'aro': 'ARO'}


def _write_bpf_pdb(rows, filepath):
    with open(filepath, 'w') as f:
        for i, (x, y, z, tipo) in enumerate(rows, start=1):
            resname = _BPF_RESNAME[tipo]
            f.write(
                f"ATOM  {i:>5} {'H':<4} {resname:>3} X{i:>4}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {'H':>2}\n"
            )
        f.write('END\n')


def export_bpf(DF_Lig, filepath, DF_true=None):
    """Genera un archivo .bpf (bias probe file, formato GOLD: header 'x y z
    Vset r type') a partir de los hot-points del ligando (ver
    _collect_bias_points; DF_true filtra a solo los validados). Vset y r son
    fijos por tipo (_BPF_PARAMS)."""
    _write_bpf(_collect_bias_points(DF_Lig, DF_true), filepath)


def export_bpf_pdb(DF_Lig, filepath, DF_true=None):
    """PDB 'dummy' con un átomo H por punto de bias (mismos puntos y mismo
    orden que export_bpf; DF_true filtra a solo los validados), para poder
    cargar y visualizar los puntos de bias en VMD junto al receptor/ligando.
    resname = DON/ACC/ARO según el tipo, chain 'X', un residuo por átomo
    (dummy, sin significado bioquímico)."""
    _write_bpf_pdb(_collect_bias_points(DF_Lig, DF_true), filepath)


_RECEPTOR_TYPE_MAP = {'Aceptor': 'acc', 'Dador': 'don', 'aromatic': 'aro'}


def _collect_bias_points_receptor(receptor_points):
    """Misma idea que _collect_bias_points pero a partir de receptor_points
    (salida de Coordenadas_interes_receptor en receptor_site.py: columnas
    Type/Pos/Residue/Atom/X/Y/Z). Un punto por átomo aceptor/donor del
    receptor y uno por centroide de anillo aromático (ya viene un solo row
    por anillo, Atom == 'center')."""
    rows = []
    for _, r in receptor_points.iterrows():
        tipo = _RECEPTOR_TYPE_MAP.get(r['Type'])
        if tipo is None:
            continue
        rows.append((float(r['X']), float(r['Y']), float(r['Z']), tipo))
    return rows


def export_bpf_receptor(receptor_points, filepath):
    """Igual que export_bpf pero para puntos de bias del *receptor* (ver
    _collect_bias_points_receptor) -- usado por el modo "site bias" para
    generar un .bpf alrededor de una coordenada arbitraria, sin ligando."""
    _write_bpf(_collect_bias_points_receptor(receptor_points), filepath)


def export_bpf_pdb_receptor(receptor_points, filepath):
    """Equivalente a export_bpf_pdb para puntos de bias del receptor."""
    _write_bpf_pdb(_collect_bias_points_receptor(receptor_points), filepath)


def export_bpf_points(rows, filepath):
    """.bpf a partir de una lista ya armada de (x, y, z, tipo) -- usada por
    el método 'ideal' de analyze_site_bias (ideal_sites.py), que arma las
    filas con la geometría completa (abanico) en vez de un punto por átomo."""
    _write_bpf(rows, filepath)


def export_bpf_pdb_points(rows, filepath):
    """Equivalente a export_bpf_points pero en formato PDB dummy."""
    _write_bpf_pdb(rows, filepath)
