from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, Field, ValidationError, model_validator

__all__ = [
    "YesNo",
    "Options",
    "Distances",
    "Angles",
    "Aromaticity",
    "Pockets",
    "InteractionConfig",
    "load_config",
    "ValidationError",
]

YesNo = Literal["Yes", "No"]  # flag: convert to bool in Phase 14

_DEFAULT_CONFIG_PATH = Path(__file__).parent.parent.parent / "Interacciones_variables.yml"


class Options(BaseModel):
    ligand_plot: YesNo = Field(default="Yes", description="Generate 2D ligand PNG plot")
    vmd_output: YesNo = Field(default="Yes", description="Generate VMD TCL visualization script")
    cumulative_output: YesNo = Field(default="Yes", description="Append results to cumulative CSV files")


class Distances(BaseModel):
    Distances_Hidrogen_Bonds: float = Field(default=3.2, gt=0, description="H-bond distance cutoff (Å)")
    Distances_Aromatic: float = Field(default=5.5, gt=0, description="Aromatic interaction distance cutoff (Å)")
    Distances_Hidrofobica: float = Field(default=4.0, gt=0, description="Hydrophobic interaction distance cutoff (Å)")
    centroid_distance: float = Field(default=12.0, gt=0, description="Active-site search radius (Å)")
    Distances_C_Simple: float = Field(default=1.54, gt=0, description="C–C single bond distance cutoff (Å)")
    Distances_C_Doble: float = Field(default=2.56, gt=0, description="C=C double bond distance cutoff (Å)")

    @model_validator(mode="before")
    @classmethod
    def _strip_key_whitespace(cls, data: object) -> object:
        # YAML source has trailing spaces on some keys (e.g. "Distances_C_Simple ")
        if isinstance(data, dict):
            return {k.strip(): v for k, v in data.items()}
        return data


class Angles(BaseModel):
    Angle_Hidrogen_Bonds_Min: float = Field(default=100.0, ge=0, le=360, description="Minimum H-bond angle (°)")
    Angle_Hidrogen_Bonds_Max: float = Field(default=180.0, ge=0, le=360, description="Maximum H-bond angle (°)")


class Aromaticity(BaseModel):
    Ring_Planarity_RMSD_Max: float = Field(
        default=0.15, gt=0, description="Maximum ring-plane RMSD to classify ring as aromatic (Å)"
    )


class Pockets(BaseModel):
    min_residues: int = Field(default=3, gt=0, description="Minimum distinct residues contacting a ligand fragment to qualify as a pocket")
    coverage_threshold: float = Field(default=0.5, ge=0, le=1, description="Maximum coverage R (0=fully surrounded, 1=one-sided) to qualify as a pocket")


class InteractionConfig(BaseModel):
    options: Options = Field(default_factory=Options)
    distancias: Distances = Field(default_factory=Distances)
    angulos: Angles = Field(default_factory=Angles)
    aromaticidad: Aromaticity = Field(default_factory=Aromaticity)
    pockets: Pockets = Field(default_factory=Pockets)
    acceptors: dict[str, list[str]]
    donors: dict[str, list[str]]
    acceptors_antecedent: dict[str, dict[str, str]]
    special: dict[str, list] = Field(default_factory=dict)


def load_config(path: Path | str | None = None) -> InteractionConfig:
    """Read a YAML config file and return a validated InteractionConfig.

    Defaults to Interacciones_variables.yml at the project root when *path* is omitted.
    Raises pydantic.ValidationError on invalid values.
    """
    resolved = Path(path) if path is not None else _DEFAULT_CONFIG_PATH
    with resolved.open(encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return InteractionConfig.model_validate(data)
