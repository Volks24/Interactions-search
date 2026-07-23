"""Tests for src/interactions_search/config.py — Step 1.1 acceptance criteria."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from interactions_search.config import (
    Distances,
    InteractionConfig,
    Options,
    load_config,
)

_YAML_PATH = Path(__file__).parent.parent.parent / "Interacciones_variables.yml"


def test_valid_yaml_loads():
    config = load_config(_YAML_PATH)
    assert isinstance(config, InteractionConfig)
    assert config.distancias.Distances_Hidrogen_Bonds == pytest.approx(3.2)
    assert config.options.ligand_plot == "Yes"


def test_negative_distance_raises():
    with pytest.raises(ValidationError) as exc_info:
        Distances(
            Distances_Hidrogen_Bonds=-1.0,
            Distances_Aromatic=5.5,
            Distances_Hidrofobica=4.0,
            centroid_distance=12.0,
            Distances_C_Simple=1.54,
            Distances_C_Doble=2.56,
        )
    assert "Distances_Hidrogen_Bonds" in str(exc_info.value)


def test_zero_distance_raises():
    with pytest.raises(ValidationError):
        Distances(
            Distances_Hidrogen_Bonds=0.0,
            Distances_Aromatic=5.5,
            Distances_Hidrofobica=4.0,
            centroid_distance=12.0,
            Distances_C_Simple=1.54,
            Distances_C_Doble=2.56,
        )


def test_invalid_yes_no_raises():
    with pytest.raises(ValidationError) as exc_info:
        Options(
            ligand_plot="maybe",
            vmd_output="Yes",
            cumulative_output="No",
        )
    assert "ligand_plot" in str(exc_info.value)


def test_angle_out_of_range_raises():
    from interactions_search.config import Angles

    with pytest.raises(ValidationError):
        Angles(Angle_Hidrogen_Bonds_Min=-5.0, Angle_Hidrogen_Bonds_Max=180.0)

    with pytest.raises(ValidationError):
        Angles(Angle_Hidrogen_Bonds_Min=0.0, Angle_Hidrogen_Bonds_Max=400.0)


def test_load_config_default_path():
    """load_config() with no argument resolves to the project-root YAML."""
    config = load_config()
    assert isinstance(config, InteractionConfig)


def test_yaml_key_trailing_spaces_tolerated():
    """Distances_C_Simple and Distances_C_Doble have trailing spaces in the YAML source."""
    config = load_config(_YAML_PATH)
    assert config.distancias.Distances_C_Simple == pytest.approx(1.54)
    assert config.distancias.Distances_C_Doble == pytest.approx(2.56)
