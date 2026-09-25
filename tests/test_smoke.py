"""Smoke tests — freeze current behavior of Interactions_search.py.

Each test uses a module-scoped fixture that runs the script once against the
minimal fixture PDBs and cleans up all output afterward.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).parent.parent
PYTHON = sys.executable
SCRIPT = REPO_ROOT / "Interactions_search.py"
RECEPTOR = REPO_ROOT / "tests" / "fixtures" / "receptor_mini.pdb"
LIGAND = REPO_ROOT / "tests" / "fixtures" / "ligand_mini.pdb"

RECEPTOR_NOINT = REPO_ROOT / "tests" / "fixtures" / "receptor_noint.pdb"
LIGAND_NOINT = REPO_ROOT / "tests" / "fixtures" / "ligand_noint.pdb"
OUTPUT_FOLDER = REPO_ROOT / "receptor_mini_ligand_mini"
OUTPUT_FOLDER_NOINT = REPO_ROOT / "receptor_noint_ligand_noint"

EXPECTED_COLUMNS = ["Pos R", "Res", "Atom", "Dist", "Lig", "Type", "Angle",
                    "X", "Y", "Z", "Interaction", "Reason"]


@pytest.fixture(scope="module", autouse=True)
def clean_module(tmp_path_factory):
    """Run against temporary input copies; never delete user outputs or edit fixtures."""
    global RECEPTOR, LIGAND, RECEPTOR_NOINT, LIGAND_NOINT
    global OUTPUT_FOLDER, OUTPUT_FOLDER_NOINT
    workspace = tmp_path_factory.mktemp('smoke')
    copied = []
    for source in (RECEPTOR, LIGAND, RECEPTOR_NOINT, LIGAND_NOINT):
        target = workspace / source.name
        shutil.copy(source, target)
        copied.append(target)
    RECEPTOR, LIGAND, RECEPTOR_NOINT, LIGAND_NOINT = copied
    OUTPUT_FOLDER = workspace / OUTPUT_FOLDER.name
    OUTPUT_FOLDER_NOINT = workspace / OUTPUT_FOLDER_NOINT.name


@pytest.fixture(scope="module")
def run_output():
    return subprocess.run(
        [
            str(PYTHON), str(SCRIPT),
            "-r", str(RECEPTOR),
            "-l", str(LIGAND),
            "-c", "A",
        ],
        capture_output=True,
        text=True,
        cwd=str(OUTPUT_FOLDER.parent),
    )
    

@pytest.fixture(scope="module")
def run_output_noint():
    return subprocess.run(
        [
            str(PYTHON), str(SCRIPT),
            "-r", str(RECEPTOR_NOINT),
            "-l", str(LIGAND_NOINT),
            "-c", "A",
        ],
        capture_output=True,
        text=True,
        cwd=str(OUTPUT_FOLDER_NOINT.parent),
    )


def test_exit_code_zero(run_output, run_output_noint):
    assert run_output.returncode == 0, (
        f"Script exited {run_output.returncode}\n--- stderr ---\n{run_output.stderr}"
    )
    assert run_output_noint.returncode == 0, (
        f"Script exited {run_output_noint.returncode}\n--- stderr ---\n{run_output_noint.stderr}"
    )


def test_output_folder_and_all_csv_exist(run_output, run_output_noint):
    assert OUTPUT_FOLDER.is_dir(), "Output folder was not created"
    csv_files = list(OUTPUT_FOLDER.glob("Interaction_*_all.csv"))
    assert csv_files, "No Interaction_*_all.csv found in output folder"

    assert OUTPUT_FOLDER_NOINT.is_dir(), "Output folder for no interaction was not created"
    csv_files_noint = list(OUTPUT_FOLDER_NOINT.glob("Interaction_*_all.csv"))
    assert csv_files_noint, "No Interaction_*_all.csv found in no interaction output folder"


def test_all_csv_has_expected_columns(run_output, run_output_noint):
    csv = next(OUTPUT_FOLDER.glob("Interaction_*_all.csv"))
    df = pd.read_csv(csv, index_col=0)
    assert list(df.columns) == EXPECTED_COLUMNS, (
        f"Column mismatch.\n  got : {list(df.columns)}\n  want: {EXPECTED_COLUMNS}"
    )

    csv_noint = next(OUTPUT_FOLDER_NOINT.glob("Interaction_*_all.csv"))
    df_noint = pd.read_csv(csv_noint, index_col=0)
    assert list(df_noint.columns) == EXPECTED_COLUMNS, (
        f"Column mismatch.\n  got : {list(df_noint.columns)}\n  want: {EXPECTED_COLUMNS}"
    )


def test_true_csv_interactions_all_yes(run_output, run_output_noint):
    true_csv = next(OUTPUT_FOLDER.glob("Interaction_*_true.csv"))
    df = pd.read_csv(true_csv, index_col=0)
    assert not df.empty, "_true.csv is empty — expected at least one validated interaction"
    bad = df[df["Interaction"] != "Yes"]
    assert bad.empty, (
        f"{len(bad)} row(s) in _true.csv have Interaction != 'Yes':\n{bad}"
    )

    true_csv_noint = next(OUTPUT_FOLDER_NOINT.glob("Interaction_*_true.csv"))
    df_noint = pd.read_csv(true_csv_noint, index_col=0)
    assert df_noint.empty, "_true.csv is not empty — expected no validated interactions"
