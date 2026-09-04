"""Ángulos chi (torsión de cadena lateral, chi1-chi5) de los residuos que
participan en un pocket hidrofóbico validado. Los cuatro átomos que definen
cada chi por residuo vienen del recurso empaquetado data/chi_angles.json
(`_CHI_ATOMS`) — es una tabla de referencia IUPAC estática, no configuración
de proyecto, por eso vive dentro del paquete (a diferencia de
Interacciones_variables.yml, que el usuario sí edita por proyecto) y se
carga con importlib.resources para que funcione también desde un wheel
instalado, no solo en un clone editable del repo. El diedro se calcula con
geometry.dihedral_angle sobre las coordenadas ya presentes en
DF_Active_Site (active_site_residues() guarda todos los átomos del residuo,
no solo los puntos de interés, así que no hace falta releer el PDB)."""
from __future__ import annotations

import json
import re
from importlib import resources

import numpy as np
import pandas as pd

from interactions_search.geometry import dihedral_angle

__all__ = ["compute_residue_chi_angles", "compute_pocket_chi_angles"]

_CHI_NAMES = ["chi1", "chi2", "chi3", "chi4", "chi5"]

# Variantes de protonación de histidina (usadas en Interacciones_variables.yml)
# comparten la definición de chi de HIS — chi_angles.json solo lista HIS.
_RESNAME_ALIASES = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS"}


def _load_chi_atoms() -> dict[str, dict[str, tuple[str, str, str, str]]]:
    resource = resources.files("interactions_search.data").joinpath("chi_angles.json")
    raw = resource.read_text(encoding="utf-8")
    data = json.loads(raw)["chi_angles"]
    by_residue: dict[str, dict[str, tuple[str, str, str, str]]] = {}
    for chi_name, per_residue in data.items():
        for resname, spec in per_residue.items():
            by_residue.setdefault(resname, {})[chi_name] = tuple(spec.split("-"))
    return by_residue


_CHI_ATOMS = _load_chi_atoms()


def compute_residue_chi_angles(DF_Active_Site, pos, resname):
    """Calcula chi1-chi5 (grados, redondeado; None si el residuo no tiene ese
    chi o falta alguno de los 4 átomos en DF_Active_Site, ej. PDB incompleto)
    para el residuo `resname` en la posición `pos`."""
    chi_defs = _CHI_ATOMS.get(_RESNAME_ALIASES.get(resname, resname), {})
    if not chi_defs:
        return {chi: None for chi in _CHI_NAMES}

    atoms = DF_Active_Site[DF_Active_Site["Pos"] == pos]
    coords = {row["Atom"]: np.array([row["X"], row["Y"], row["Z"]], dtype=float)
              for _, row in atoms.iterrows()}

    result = {}
    for chi in _CHI_NAMES:
        spec = chi_defs.get(chi)
        if spec is None or any(a not in coords for a in spec):
            result[chi] = None
            continue
        p0, p1, p2, p3 = (coords[a] for a in spec)
        result[chi] = round(dihedral_angle(p0, p1, p2, p3), 2)
    return result


_RESIDUE_TOKEN_RE = re.compile(r"^([A-Z]+)(\d+)$")


def compute_pocket_chi_angles(df_pocket_summary, DF_Active_Site):
    """Una fila por (Pocket, residuo) con chi1..chi5, solo para los pockets
    validados (Is_Pocket == 'Yes'). Parsea la columna 'Residues' del summary
    (ej. 'LEU63,VAL67,TYR129') en vez de recibir los residuos ya agrupados,
    para no acoplar pockets.py a este cálculo."""
    cols = ["Pocket", "Pos", "Residue", *_CHI_NAMES]
    if df_pocket_summary.empty:
        return pd.DataFrame(columns=cols)

    rows = []
    qualifying = df_pocket_summary[df_pocket_summary["Is_Pocket"] == "Yes"]
    for _, prow in qualifying.iterrows():
        for token in prow["Residues"].split(","):
            m = _RESIDUE_TOKEN_RE.match(token)
            if not m:
                continue
            resname, pos = m.group(1), int(m.group(2))
            chis = compute_residue_chi_angles(DF_Active_Site, pos, resname)
            rows.append([prow["Pocket"], pos, resname, *(chis[c] for c in _CHI_NAMES)])

    return pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)
