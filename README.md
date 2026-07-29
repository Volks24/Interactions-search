# Interactions Search

Detects and classifies non-covalent interactions between a ligand and a protein from PDB files. Identifies hydrogen bonds, aromatic interactions (π-π and T-shaped), hydrophobic contacts, salt bridges, and π-cation interactions. Each contact is validated by distance and angle criteria. Outputs CSV files and TCL scripts for visualization in VMD.

---

## Files

| File | Description |
|---|---|
| `Interactions_search.py` | Compatibility shim at the repo root — `python Interactions_search.py ...` still works exactly as before. All logic lives in the `interactions_search` package (below); this file just re-exports the CLI entry point. |
| `Interacciones_variables.yml` | Distance thresholds, acceptors and donors per residue |
| `src/interactions_search/align.py` | Structural alignment utility (`align-protein` CLI), see [Structural Alignment](#structural-alignment-alignpy) below |

### Package layout (`src/interactions_search/`)

The pipeline described below is implemented as one module per stage:

| Module | Responsibility |
|---|---|
| `config.py` | Loads and validates `Interacciones_variables.yml` (pydantic) |
| `geometry.py` | Pure geometric helpers: center of mass, angles, ring centroids, planarity, convex-hull volume |
| `io_pdb.py` | PDB reading, `split_pdb`, input validation |
| `ligand_hotpoints.py` | Ligand acceptor/donor/aromatic hot-points (RDKit + SMARTS) and their 2D PNGs |
| `receptor_site.py` | Receptor active-site residues and their points of interest |
| `contacts.py` | Distance-based contact search (H-bond, hydrophobic, salt bridge, π-cation) and angle validation |
| `bias.py` | GOLD bias probe file (`.bpf`) export |
| `pockets.py` | Hydrophobic pocket detection (`search_hydrophobic_pockets`) |
| `plotting.py` | Convex-hull PNGs (scatter + solid surface) |
| `vmd.py` | VMD `.tcl` script generation |
| `pipeline.py` | Orchestrates all of the above into `analyze_pair()` |
| `cli.py` | Argument parsing and batch/complex-PDB pair resolution (`main()`) |
| `align.py` | Structural alignment utility (`align-protein` CLI) |

This is a structural split only — no behavior changed; `tests/test_smoke.py` runs the full pipeline end-to-end against fixture PDBs to guard against regressions.

---

## Installation

### Dev install (recommended)

```bash
git clone <repo-url>
cd Interactions-search
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

Installs the package in editable mode plus dev tools (pytest, ruff, mypy).

### Script-only (no package install)

If you only need to run the script without the package:

```bash
pip install biopython rdkit pandas numpy pyyaml scipy matplotlib
python Interactions_search.py -r receptor.pdb -l ligand.pdb -c A
```

---

## Usage

### Mode 1 — Separate PDB files

```bash
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A
```

### Mode 2 — Batch (multiple ligands, one receptor)

```bash
python Interactions_search.py -r protein.pdb -l lig1.pdb lig2.pdb lig3.pdb -c A
```

Each pair generates its own output folder. If `cumulative_output: 'Yes'` in the config, `Interactions_close.csv` and `CM_all.csv` (in the current directory) get one row appended per pair, tagged by `Receptor`/`Ligand`, so successive or batch runs accumulate without overwriting each other.

### Mode 3 — Complex PDB

A single PDB containing protein + ligand(s). The script splits it automatically.

```bash
# Single HETATM group (selected automatically)
python Interactions_search.py -x complex.pdb -c A

# Multiple HETATM groups (select with -n)
python Interactions_search.py -x complex.pdb -c A -n LIG
```

If multiple HETATM groups are present and `-n` is not specified, the script lists the available names and exits without analysing.

If the complex PDB has no `HETATM` records at all (ligand saved as `ATOM`, e.g. some CHARMM/AMBER-prepped structures), use `-f` to tell `split_pdb` which residue name(s) to treat as ligand:

```bash
# Ligand TF3 is stored as ATOM records, not HETATM
python Interactions_search.py -x complex.pdb -c A -f TF3

# Multiple forced ligand names; still use -n to pick one for this run
python Interactions_search.py -x complex.pdb -c A -f TF3 7FW -n TF3
```

### Arguments

| Argument | Description |
|---|---|
| `-x / --complex` | Complex PDB. Alternative to `-r`. |
| `-r / --receptor_pdb` | Receptor PDB (already separated). |
| `-l / --ligand_input` | Ligand PDB(s). Accepts one or several (batch). |
| `-c / --chain_receptor` | Protein chain (e.g. `A`). |
| `-n / --lig_name` | HETATM name when multiple groups exist in `--complex`. |
| `-f / --force_ligand` | Residue name(s) to treat as ligand even if stored as `ATOM` instead of `HETATM` in `--complex`. |

---

## Python API

> **Status: in development.** Programmatic access is being extracted into the `interactions_search` package. The script interface above is the stable way to run analyses for now.

Once available, the intended usage will be:

```python
from interactions_search import analyze_pair, load_config

cfg = load_config("Interacciones_variables.yml")
analyze_pair("receptor.pdb", "ligand.pdb", chain="A", cfg=cfg)
```

---

## Analysis Pipeline

```
Complex PDB (optional)
        │
        ▼
[1] split_pdb()
    ├── <stem>_protein.pdb      ← ATOM records
    └── <stem>_RESNAME.pdb      ← HETATM records per residue (water excluded)
        │
        ▼
[2] Ligand cleanup
    └── remove_bias()           ← removes CM atoms from ligand PDB
        │
        ▼
[3] Ligand hot-points  (RDKit + SMARTS)
    ├── H-bond acceptors:  [O;H1], [O;H0], [N;H1], [N;H0], [n], [o], [N+]
    ├── H-bond donors:     [O;H], [N;H2], [N;H], [S;H], [nH]
    └── Aromatic rings:    detected via ring_info, filtered by size > 5 and by
                           planarity (best-fit-plane RMSD ≤ Ring_Planarity_RMSD_Max)
        │
        ▼
[4] Receptor active site  (BioPython)
    └── active_site_residues()
        Residues whose centre of mass is within 12 Å of the ligand CM
        HOH and the ligand itself are excluded
        │
        ▼
[5] Receptor points of interest
    └── Coordenadas_interes_receptor()
        ├── Receptor acceptors  (from YAML table)
        ├── Receptor donors     (from YAML table)
        └── Aromatic ring centroids (TYR, PHE, TRP)
        │
        ▼
[6] Contact search  (numpy, vectorised distances)
    ├── H-bond:        lig acceptor ↔ rec donor      (threshold: Distances_Hidrogen_Bonds)
    │                  lig donor    ↔ rec acceptor    (threshold: Distances_Hidrogen_Bonds)
    ├── Aromatic:      centroid ↔ centroid             (threshold: Distances_Aromatic)
    ├── Hydrophobic:   apolar C lig ↔ apolar C rec     (threshold: Distances_Hidrofobica)
    ├── Salt bridge:   ± group lig ↔ ∓ group rec       (threshold: 4.0 Å)
    └── π-cation:      lig ring ↔ ARG/LYS/HIS rec      (threshold: 5.0 Å)
        │
        ▼
[7] Angle validation
    ├── H-bond:    D-A···Antecedent angle between 100° and 200°
    ├── Aromatic:  angle between ring planes
    │               0°–30°  → π-π (parallel / sandwich)
    │               60°–90° → T-shaped (perpendicular)
    └── Hyd. / salt / π-cat:  validated by distance only
        │
        ▼
[8] Hydrophobic pocket detection  (search_hydrophobic_pockets, independent of step 7)
    ├── Group ligand hydrophobic atoms into fragments by bond connectivity
    ├── Per fragment, collect distinct contacting receptor residues (≥ Pocket_Min_Residues)
    ├── Score spatial coverage around the fragment (Coverage_R, see "Hydrophobic Pockets" below)
    ├── Convex-hull volume of the contacting residues' atoms (Volume_A3, scipy.spatial.ConvexHull)
    └── Local hydrophobic contact density, fpocket-style (Density_Score)
        │
        ▼
[9] Outputs
    ├── CSV with all raw interactions
    ├── CSV filtered by distance
    ├── CSV filtered by distance + angle  (validated interactions)
    ├── CSV of hydrophobic pocket candidates (Pockets_<rec>_<lig>.csv)
    ├── Console summary  (table of validated interactions + pocket count)
    ├── Cumulative summary Interactions_close.csv
    ├── Centre of mass    CM_all.csv
    └── TCL scripts for VMD (if vmd_output: Yes)
```

---

## PDB Splitting (`split_pdb`)

When `--complex` is used, `split_pdb()` pre-processes the PDB before analysis:

- Separates `ATOM` records (protein) from `HETATM` records (ligands/cofactors)
- Groups HETATM records by residue name (`resName`)
- **Excludes water** automatically: HOH, WAT, TIP, TIP3, SOL, DOD
- Distributes `CONECT` records to the corresponding HETATM file (by atom serial)
- Saves files to a temporary directory that is deleted after the run

Example with a PDB containing protein + ligand LIG + HEM group:

```
<tmpdir>/
├── complex_protein.pdb   ← full protein
├── complex_LIG.pdb       ← organic ligand + its CONECT records
└── complex_HEM.pdb       ← haem group
```

---

## Structural Alignment (`align.py`)

Installed as the package's `interactions_search.align` module, with an `align-protein`
console script (available after `pip install -e .`). Given a reference PDB and a protein
PDB, aligns the protein onto the reference by CA-atom superposition (BioPython
`Superimposer`) over the overlapping residue range of a chosen chain, then writes out the
aligned structure and reports the RMSD. Useful as a pre-processing step before running
`Interactions_search.py` on structures that need to be in the same frame (e.g. comparing
poses/ligands across multiple crystal structures of the same protein).

```bash
align-protein -R reference.pdb -P protein.pdb -C A

# keep only chains A and B in the aligned output (the alignment chain -C is always kept)
align-protein -R reference.pdb -P protein.pdb -C A -K A B
```

### Arguments

| Argument | Description |
|---|---|
| `-R / --reference` | Reference PDB (the structure being aligned onto). |
| `-P / --protein` | Protein PDB to align. |
| `-C / --Chain` | Chain identifier used for alignment (default: `A`). Must exist in both PDBs. |
| `-K / --keep-chains` | Chain(s) to keep in the aligned output PDB (drops the rest). The alignment chain (`-C`) is always kept even if omitted here. |

### How it works

1. Loads both PDBs and selects chain `-C` from each.
2. Restricts the fit to the residue-number range common to both chains (`max(first)` to
   `min(last)`), matching residues by number (not list index) so gaps or missing residues
   don't misalign the selection.
3. Superimposes using each residue's `CA` atom in that range; raises if the atom counts
   don't match or no CA atoms are found in the range.
4. Applies the resulting transformation to **all** atoms of the protein structure (every
   chain, not just `-C`), then optionally drops chains not listed in `-K`.

### Outputs

| File | Content |
|---|---|
| `<protein>_alig.pdb` | Aligned structure (all kept chains, transformed) |
| `<protein>_resultados.txt` | Reference/protein names, aligned residue range, and RMSD (Å) |

---

## Configuration (`Interacciones_variables.yml`)

```yaml
options:
  ligand_plot: 'Yes'         # generates PNG images of ligand acceptors, donors and rings
  vmd_output:  'Yes'         # generates TCL script for VMD visualisation
  cumulative_output: 'Yes'   # appends each pair to Interactions_close.csv / CM_all.csv
  interaction_coord: 'center'  # X,Y,Z column added to Interaction_*.csv: 'receptor', 'ligand' or 'center' (midpoint)
  volume_plot: 'Yes'         # generates a 3D scatter PNG (convex-hull volume) for the whole
                              # active site and for each qualifying hydrophobic pocket (needs matplotlib)
  bias: 'No'                 # generates <rec>_<lig>.bpf (GOLD bias probe file) and
                              # <rec>_<lig>_bias.pdb (same points, dummy atoms for VMD)
  bias_validated_only: 'No'  # 'No': all of the ligand's chemical hot-points. 'Yes': only the
                              # ones that form a validated interaction (Interaction == 'Yes')

distancias:
  Distances_Hidrogen_Bonds: 3.2   # Å — H-bond threshold
  Distances_Aromatic:       5.5   # Å — centre-to-centre aromatic threshold
  Distances_Hidrofobica:    4.0   # Å — hydrophobic threshold (aligned with PLIP)
  centroid_distance:       12.0   # Å — active site search radius
  Distances_C_Simple:       1.54  # Å — C-C single bond (reference)
  Distances_C_Doble:        2.56  # Å — C=C double bond (reference)

angulos:
  Angle_Hidrogen_Bonds_Min: 100    # ° — minimum Donor-Acceptor-Antecedent angle
  Angle_Hidrogen_Bonds_Max: 180    # ° — maximum angle (180° is the geometric ceiling)

aromaticidad:
  Ring_Planarity_RMSD_Max:  0.15   # Å — max RMSD to the ring's best-fit plane to
                                   # be considered aromatic (real rings ~0.01, chair ~0.25)

acceptors:             # acceptor atoms per residue
  ALA: [O]
  TYR: [O, OH]
  ASP: [O, OD1, OD2]
  ...

donors:                # donor atoms per residue
  ALA: [N]
  ARG: [N, HNE, HH11, HH12, HH21, HH22, HE, NE, NH1, NH2]
  ...

acceptors_antecedent:  # antecedent atom of each acceptor (for angle calculation)
  TYR: {OH: CZ}
  ASP: {OD1: CG, OD2: CG}
  ...

special:               # special cases (e.g. haem group)
  HEM: [FE, 1.59]

pockets:
  min_residues:        3     # minimum distinct residues contacting the same ligand fragment
  coverage_threshold:  0.5   # max Coverage_R (0-1) to qualify as an enclosing pocket
  density_radius:      5.0   # Å — neighbourhood radius for Density_Score (fpocket-style local density)
```

---

## Hydrophobic Pockets

A single hydrophobic contact (one residue, one ligand atom) doesn't tell you whether the
ligand sits in a real, enclosing binding pocket, or just brushes past a residue on one
side. `search_hydrophobic_pockets()` (in `src/interactions_search/pockets.py`) answers that question
with two independent criteria, both must hold:

1. **Multiple residues on the same ligand fragment.** Ligand hydrophobic atoms (matched by
   the same SMARTS used for `hydrophobic` contacts, `_HPHO_LIG_SMARTS`) are grouped into
   *fragments* by **bond connectivity** (RDKit's bond graph), not spatial proximity — a ring
   or a contiguous aliphatic chain that is contacted counts as one fragment. For each
   fragment, every receptor residue with at least one apolar atom within
   `Distances_Hidrofobica` Å of any atom in that fragment is collected. The fragment
   qualifies only if it has **≥ `min_residues`** distinct contacting residues (default 3).

2. **Spatial coverage around the fragment.** Having 3+ residues touching the same fragment
   isn't enough on its own — they could all be sitting on the same face of the ligand
   (a flat, superficial contact) rather than wrapping around it. `Coverage_R` measures this:
   for each contacting residue, take the unit vector from the fragment's centroid to that
   residue's centroid, then compute the magnitude of the **average of those unit vectors**.

   ```
   Coverage_R = | Σ unit_vectors | / n_residues        (0 ≤ Coverage_R ≤ 1)
   ```

   - **Coverage_R ≈ 0** — the vectors point in different directions and cancel out:
     residues surround the fragment from multiple sides → a real, enclosing pocket,
     consistent with well-defined binding sites seen in crystal structures.
   - **Coverage_R ≈ 1** — the vectors mostly point the same way: all residues are on
     the same side → a superficial contact, not an enclosing pocket, even with 3+ residues.

   A fragment is marked `Is_Pocket = Yes` only if `n_residues ≥ min_residues` **and**
   `Coverage_R < coverage_threshold` (default 0.5).

**Volume.** For every candidate fragment (pass or fail), `Volume_A3` is the volume (Å³) of the
convex hull (`scipy.spatial.ConvexHull`) of all atoms belonging to the contacting residues
(the same atom set `vmd_pockets_*.tcl`'s `Surf` representation selects) — `NaN` if there
are fewer than 4 atoms or the geometry is degenerate (coplanar points). This is separate
from the whole active site's volume (`ActiveSite_Volume_A3` in `summary.csv`), which is the
convex hull of every atom in `DF_Active_Site` (all residues within `centroid_distance` Å of
the ligand), independent of hydrophobic pockets. If `options.volume_plot: 'Yes'`, two PNGs
are generated for the active site and for each qualifying pocket (`Is_Pocket == Yes`): a
scatter view (all atoms in black, convex-hull vertices in red) and a `_solid` view (the hull
rendered as a solid triangulated surface, faces coloured by height on the viridis colormap).

**Density.** `Density_Score` is a local hydrophobic density measure (fpocket-style "mean
local hydrophobic density"), distinct from `Coverage_R` (residue direction) and `Volume_A3`
(cavity size): each atom(ligand)–atom(receptor) contact in the fragment is represented by
its midpoint; for every midpoint, count how many *other* midpoints of the same fragment fall
within `density_radius` Å (default 5.0), then average those counts. A high score means the
individual contacts are tightly clustered (a snug hydrophobic fit); a low score means they're
spread out even within the same enclosing cavity (a loose fit). `0.0` if the fragment has
fewer than 2 contacts.

This runs independently of the per-contact `Interaction == Yes` validation in step 7 — a
fragment can have several individually-validated hydrophobic contacts and still fail the
pocket criteria (e.g. only 2 residues), or vice versa.

---

## Bias Probe File

If `options.bias: 'Yes'`, two files with the ligand's H-bond/aromatic hot-points are
generated per pair: `<rec>_<lig>.bpf` (GOLD's bias probe file format — header
`x y z Vset r type`, tab-separated) and `<rec>_<lig>_bias.pdb` (the same points as dummy `H`
atoms, `resname` `DON`/`ACC`/`ARO`, chain `X`, so they can be loaded and visualised in VMD
alongside the receptor/ligand).

One point per donor atom (`don`), one per acceptor atom (`acc`), and one per aromatic ring
**centroid** (`aro`, not one point per ring atom — that would produce several near-duplicate
points for the same ring). `Vset`/`r` are fixed per type, not per atom: `don` → `-2.72/1.20`,
`acc` → `-2.28/0.80`, `aro` → `-2.00/2.00`.

By default (`bias_validated_only: 'No'`) this uses **every** acceptor/donor/aromatic
hot-point the ligand has chemically (from `search_hot_points`/`search_rings`), regardless of
whether that group actually contacts the receptor in this pose — e.g. a solvent-facing
hydroxyl still gets a bias point. Set `bias_validated_only: 'Yes'` to restrict it to the
hot-points that form a validated interaction (`Interaction == 'Yes'` in `_true.csv`):
acceptor/donor atoms are matched by serial (`LigID`), and a ring counts as validated if it
appears as `aromatic` or `pi_cation` in the validated interactions.

---

## Outputs

### CSV files

| File | Content |
|---|---|
| `<folder>/Interaction_<rec>_<lig>_all.csv` | All interactions found (no filters) |
| `<folder>/Interaction_<rec>_<lig>_threshold.csv` | Filtered by distance |
| `<folder>/Interaction_<rec>_<lig>_true.csv` | Validated by distance and angle |
| `<folder>/Pockets_<rec>_<lig>.csv` | Hydrophobic pocket candidates (see "Hydrophobic Pockets" above), one row per ligand fragment |
| `<folder>/<rec>_<lig>.bpf` | GOLD bias probe file (see "Bias Probe File" above); requires `bias: 'Yes'` |
| `<folder>/<rec>_<lig>_bias.pdb` | Same bias points as dummy PDB atoms, for visualising in VMD; requires `bias: 'Yes'` |
| `<folder>/ActiveSite_<rec>_<lig>_volume.png` | 3D scatter of the whole active site's atoms + convex-hull vertices (if `volume_plot: 'Yes'`) |
| `<folder>/ActiveSite_<rec>_<lig>_volume_solid.png` | Same active-site hull, rendered as a solid triangulated surface coloured by height (Z, viridis colormap) instead of a point scatter |
| `<folder>/Pocket_<n>_<rec>_<lig>_volume.png` | 3D scatter of pocket `<n>`'s contacting atoms + convex-hull vertices (if `volume_plot: 'Yes'`, one per qualifying pocket) |
| `<folder>/Pocket_<n>_<rec>_<lig>_volume_solid.png` | Same pocket hull, rendered as a solid triangulated surface coloured by height (Z, viridis colormap) instead of a point scatter |
| `Interactions_close.csv` | Cumulative run summary, one row per pair (same content as `summary.csv`, including per-type counts) — requires `cumulative_output: 'Yes'` |
| `CM_all.csv` | Ligand centre of mass, one row per pair — requires `cumulative_output: 'Yes'` |

`Pockets_<rec>_<lig>.csv` columns:

| Column | Description |
|---|---|
| `Pocket` | Fragment id (arbitrary, stable within the run) |
| `Fragment_Atoms` | Ligand atom names in the fragment, comma-separated |
| `N_Ligand_Atoms` | Number of ligand atoms in the fragment actually in contact |
| `Residues` | Contacting receptor residues, e.g. `LEU63,VAL67,TYR129` |
| `N_Residues` | Distinct contacting residue count |
| `Coverage_R` | Spatial coverage score, 0–1 (see above); lower = more enclosing |
| `Volume_A3` | Convex-hull volume (Å³) of the contacting residues' atoms (see "Hydrophobic Pockets" above); `NaN` if not computable |
| `Density_Score` | Local hydrophobic contact density, fpocket-style (see "Hydrophobic Pockets" above); higher = tighter fit |
| `Is_Pocket` | `Yes` / `No` — whether both criteria (`N_Residues` and `Coverage_R`) are met |
| `X`, `Y`, `Z` | Centroid of the ligand fragment atoms actually in contact (the same centroid `Coverage_R` is computed around) |

Interaction CSV columns (same schema in `_all`/`_threshold`/`_true`):

| Column | Description |
|---|---|
| `Pos R` | Receptor residue number |
| `Res` | Residue name (e.g. SER, TYR) |
| `Atom` | Receptor atom involved. For `hydrophobic`, when the same ligand atom contacts several atoms of the same residue, they are collapsed into one row and listed comma-separated (e.g. `CD1,CD2,CG`), with `Dist` averaged across them |
| `Dist` | Distance in Å |
| `Lig` | Ligand atom or ring involved |
| `Type` | Type: `acceptor`, `donor`, `aromatic`, `hydrophobic`, `salt_bridge`, `pi_cation` |
| `Angle` | Validation angle in degrees |
| `X`, `Y`, `Z` | 3D coordinate of the interaction, selected by `options.interaction_coord`: `'receptor'` (receptor atom or aromatic-ring centroid; mean of the collapsed atoms for `hydrophobic`), `'ligand'` (ligand atom or ring centroid), or `'center'` (midpoint between both, default). `NaN` if the atom/ring couldn't be resolved |
| `Interaction` | `Yes` / `No` — whether distance and angle criteria are met |

### VMD scripts (if `vmd_output: 'Yes'`)

Four `.tcl` scripts are generated per pair, each self-contained (they load the PDB copies
already saved in the same output folder, so the folder can be moved or run on a different
machine without editing paths):

| File | Content |
|---|---|
| `vmd_<rec>_<lig>.tcl` | Full protein + active site residues (Licorice) + ligand (Licorice); dashed lines with distance labels for each validated H-bond/aromatic interaction — **white** aromatic, **red** ligand-acceptor, **yellow** ligand-donor |
| `vmd_hydrophobic_<rec>_<lig>.tcl` | Same base scene; dashed **orange** lines for each validated hydrophobic contact |
| `vmd_pockets_<rec>_<lig>.tcl` | Same base scene; one `Surf` (MSMS) representation per qualifying pocket (`Is_Pocket == Yes`), colour-rotated per pocket, covering that pocket's contacting residues |
| `vmd_combined_<rec>_<lig>.tcl` | `vmd_<rec>_<lig>.tcl` and `vmd_pockets_<rec>_<lig>.tcl` merged into one scene: H-bond/aromatic dashed lines **and** the pocket `Surf` representation(s) together, so both can be inspected at once without loading two scripts. Does not include the individual hydrophobic dashed lines (the pocket surface already covers that region; adding them on top clutters the view) |

The full protein is rendered with **Lines** (not `NewCartoon`/`Tube`/`Trace`): some VMD
builds — notably early `2.0.0` alpha releases — silently truncate spline-based backbone
representations to the first ~40 residues regardless of selection, a confirmed VMD bug
unrelated to the input PDB. `Lines` draws bond-by-bond and is unaffected, so it's used as
the reliable default; switch to `NewCartoon` manually in VMD's *Graphics > Representations*
if your VMD build renders it correctly.

`Surf`/`MSMS` in VMD doesn't expose a scriptable "Wireframe" draw style (only probe radius
and resolution are settable via `mol modstyle`) — to see the pocket surface as a mesh
instead of solid, change it manually: *Graphics > Representations* → select the pocket's
`Surf` rep → *Draw style* → Wireframe/Points.

### Ligand PNG images (if `ligand_plot: Yes`)

| File | Content |
|---|---|
| `<lig>_acceptors.png` | Ligand with acceptor atoms highlighted and labelled by PDB atom name |
| `<lig>_donors.png` | Ligand with donor atoms highlighted and labelled by PDB atom name |
| `<lig>_aromatic.png` | Ligand with aromatic rings highlighted; each ring has a distinct colour and is labelled R1, R2, … |

The atom names in the PNG images match the `Lig` column in the CSV files directly.

---

## Output folder structure

Everything is stored inside a single folder per pair `<receptor>_<ligand>/`:

```
<receptor>_<ligand>/
├── <receptor>.pdb             ← copy of the receptor PDB
├── <ligand>.pdb               ← copy of the ligand PDB
├── <ligand>_old.pdb           ← pre-cleanup copy (remove_bias)
├── Interaction_*_all.csv      ← all interactions, no filter
├── Interaction_*_threshold.csv← filtered by distance
├── Interaction_*_true.csv     ← validated by distance + angle
├── Pockets_*.csv              ← hydrophobic pocket candidates
├── <rec>_<lig>.bpf            ← GOLD bias probe file (if bias: Yes)
├── <rec>_<lig>_bias.pdb       ← same bias points as dummy atoms, for VMD (if bias: Yes)
├── ActiveSite_*_volume.png    ← 3D scatter + convex hull of the whole active site (if volume_plot: Yes)
├── ActiveSite_*_volume_solid.png ← same hull as a solid coloured surface (if volume_plot: Yes)
├── Pocket_<n>_*_volume.png    ← 3D scatter + convex hull per qualifying pocket (if volume_plot: Yes)
├── Pocket_<n>_*_volume_solid.png ← same pocket hull as a solid coloured surface (if volume_plot: Yes)
├── summary.csv                ← interaction count by type
├── CM.csv                     ← ligand centre of mass
├── vmd_*.tcl                  ← main VMD script (H-bonds/aromatic) (if vmd_output: Yes)
├── vmd_hydrophobic_*.tcl      ← hydrophobic contacts VMD script (if vmd_output: Yes)
├── vmd_pockets_*.tcl          ← hydrophobic pockets VMD script (if vmd_output: Yes)
├── vmd_combined_*.tcl         ← H-bonds/aromatic + pocket surfaces in one scene (if vmd_output: Yes)
├── *_acceptors.png            ← ligand with acceptors highlighted
├── *_donors.png               ← ligand with donors highlighted
└── *_aromatic.png             ← ligand with aromatic rings highlighted
```

In batch mode each pair generates its own independent folder.

---

## Development

```bash
# Run tests
pytest tests/ -v

# Lint
ruff check .

# Type-check (once annotations are added)
mypy src/
```

Smoke tests in [`tests/test_smoke.py`](tests/test_smoke.py) run the full pipeline against minimal fixture PDBs in [`tests/fixtures/`](tests/fixtures/) and verify exit code, output folder, CSV columns, and that all validated rows have `Interaction == 'Yes'`.

---

## Notes

- The script should be run from the directory containing the PDB files, or use absolute paths.
- For batch analysis of multiple ligands, the script can be called in a shell loop; with `cumulative_output: 'Yes'`, `Interactions_close.csv` and `CM_all.csv` are appended automatically across runs (set to `'No'` to disable).
- Non-standard residues not listed in `acceptors` / `donors` in the YAML are silently skipped.
- Aromatic rings in the ligand must contain more than 5 atoms and be planar (RMSD to the best-fit plane ≤ `Ring_Planarity_RMSD_Max`) to be considered. Planarity, not RDKit's aromaticity flag, is used because `Chem.MolFromPDBFile` does not reliably perceive aromaticity from PDB files without explicit bond orders.

---

## See also

- [`docs/BASELINE.md`](docs/BASELINE.md) — full behavioral snapshot: pipeline details, config keys, CSV schema, regression check command
- [`Interacciones_variables.yml`](Interacciones_variables.yml) — live configuration file with all thresholds and per-residue donor/acceptor tables
