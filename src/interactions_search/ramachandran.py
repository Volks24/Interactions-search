"""Ángulos de backbone phi/psi (Ramachandran) de los residuos del sitio
activo. A diferencia de los chi (chi_angles.py), que solo necesitan los
átomos del propio residuo, phi/psi requieren el C del residuo anterior y el N
del residuo siguiente en la cadena — por eso se buscan en la estructura
completa (`structure`), no en DF_Active_Site, que puede no incluir al vecino
si cae justo fuera del radio del sitio activo."""
from __future__ import annotations

import pandas as pd

from interactions_search.geometry import dihedral_angle

__all__ = ["compute_active_site_phi_psi"]

_COLS = ["Pos", "Residue", "phi", "psi"]


def compute_active_site_phi_psi(structure, chain_receptor, DF_Active_Site):
    """phi/psi (grados, redondeado; None si falta el vecino, hay un gap en la
    numeración de residuos, o falta algún átomo de backbone) para cada
    residuo de DF_Active_Site, en el orden en que aparecen en la cadena
    `chain_receptor` de `structure`."""
    if DF_Active_Site.empty:
        return pd.DataFrame(columns=_COLS)

    chain = structure[0][chain_receptor]
    residues = sorted((r for r in chain.get_residues() if r.id[0] == ' '),
                       key=lambda r: r.id[1])
    pos_to_index = {r.id[1]: i for i, r in enumerate(residues)}

    def dihedral_or_none(prev_ok, atoms):
        if not prev_ok or any(a is None for a in atoms):
            return None
        return round(dihedral_angle(*atoms), 2)

    rows = []
    targets = DF_Active_Site[["Pos", "Residue"]].drop_duplicates().sort_values("Pos")
    for _, row in targets.iterrows():
        pos, resname = int(row["Pos"]), row["Residue"]
        idx = pos_to_index.get(pos)
        if idx is None:
            rows.append([pos, resname, None, None])
            continue

        res = residues[idx]
        prev_res = residues[idx - 1] if idx > 0 else None
        next_res = residues[idx + 1] if idx < len(residues) - 1 else None
        # Solo calcula si el vecino es realmente el residuo contiguo en la
        # secuencia (id == pos ± 1): si hay un gap (residuo faltante en el
        # cristal), no hay enlace peptídico real que definir.
        has_prev = prev_res is not None and prev_res.id[1] == pos - 1
        has_next = next_res is not None and next_res.id[1] == pos + 1

        phi = dihedral_or_none(has_prev, [
            prev_res['C'].coord if has_prev and 'C' in prev_res else None,
            res['N'].coord if 'N' in res else None,
            res['CA'].coord if 'CA' in res else None,
            res['C'].coord if 'C' in res else None,
        ])
        psi = dihedral_or_none(has_next, [
            res['N'].coord if 'N' in res else None,
            res['CA'].coord if 'CA' in res else None,
            res['C'].coord if 'C' in res else None,
            next_res['N'].coord if has_next and 'N' in next_res else None,
        ])
        rows.append([pos, resname, phi, psi])

    return pd.DataFrame(rows, columns=_COLS)
