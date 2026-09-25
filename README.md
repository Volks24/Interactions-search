# Interactions Search

Detects and classifies non-covalent interactions between a ligand and a protein from PDB files. Identifies hydrogen bonds, aromatic interactions (π-π and T-shaped), hydrophobic contacts, salt bridges, and π-cation interactions. Each contact is validated by distance and angle criteria. Outputs CSV files and TCL scripts for visualization in VMD.

---

## Files

| File | Description |
|---|---|
| `Interactions_search.py` | Compatibility shim at the repo root — `python Interactions_search.py ...` still works exactly as before. All logic lives in the `interactions_search` package (below); this file just re-exports the CLI entry point. |
| `Interacciones_variables.yml` | Distance thresholds, acceptors and donors per residue |
| `src/interactions_search/data/chi_angles.json` | Side-chain chi1-chi5 dihedral atom definitions per residue (static IUPAC reference data, packaged with the module — not project config), used by `chi_angles.py` to compute pocket residue rotamer angles |
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
| `interaction_rules.py` | Shared neighbor search, configurable cutoffs, H-bond/aromatic validation and receptor acceptor geometry for ligand and probe modes |
| `bias.py` | GOLD bias probe file (`.bpf`) export |
| `pockets.py` | Hydrophobic pocket detection (`search_hydrophobic_pockets`) |
| `chi_angles.py` | Side-chain chi1-chi5 angles of pocket-fragment residues (`compute_pocket_chi_angles`) and of every active-site residue (`compute_active_site_chi_angles`), from `chi_angles.json` |
| `ramachandran.py` | Backbone phi/psi angles of every active-site residue (`compute_active_site_phi_psi`) |
| `probe.py` | Probe mode: simulated interactions of arbitrary coordinates (`probe_interactions`, `read_probe_file`) |
| `hotspot_pocket.py` | Hotspot mode: pockets + annotated grid from MD cosolvent hotspot clusters (`analyze_hotspot_pockets`) |
| `plotting.py` | Convex-hull PNGs (scatter + solid surface) and the Ramachandran scatter PNG |
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

### Optional ligand chemistry reference

PDB-only commands remain supported; no SDF or SMILES is required. Optionally,
provide the chemical structure of the same ligand to assign bond orders, formal
charges and aromaticity while retaining PDB coordinates, atom names and order:

```bash
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A --ligand-smiles 'CC(=O)N'
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A --ligand-sdf ligand.sdf
python Interactions_search.py -x complex.pdb -n LIG -c A --ligand-sdf ligand.sdf
```

Use only one reference option, with one ligand per execution. An SDF must contain
exactly one molecule. The complete heavy-atom connectivity must match the PDB;
invalid references or incompatible explicit hydrogens cause an error, not a
silent fallback. No atoms are added, removed or repositioned. Equivalent graph
matches prefer compatibility with explicit H, then existing multiple bonds/charges
in the PDB, with first-match tie breaking and a warning. More than 1000 matches
is rejected as ambiguous. Check symmetric groups and protonation assignments.
The reference supplies chemistry, not stereochemical
validation or a new pose.

With a reference, donor/acceptor detection uses RDKit Lipinski HBA/HBD SMARTS on
a copy with hydrogens collapsed, plus protonated aromatic N-H donors,
and rings must be chemically aromatic as well as planar. Without a reference,
the existing SMARTS/Open Babel donor/acceptor procedure and geometric ring
approximation remain in use. Planarity alone is not proof of aromaticity.

Both paths now include five-member rings. To reproduce the previous PDB-only
ring-size filter, omit the reference and add `--legacy-rings` (more than five atoms).
Output filenames are unchanged. Interaction CSV schema version 2 appends a `Reason`
column; existing columns retain their order. Readers requiring an exact column list
must include the new column.

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

### Mode 4 — Site bias (no ligand)

Given a receptor and an arbitrary coordinate (e.g. where a ligand is expected to sit —
a docking box centre, a cavity found by another tool), finds every receptor residue
within a radius (default 10 Å) of that point and exports their acceptor/donor/aromatic
points as a `.bpf` + dummy PDB, using the same receptor-site logic as the regular
pipeline (`active_site_residues()` / `Coordenadas_interes_receptor()` in
`receptor_site.py`) but centred on the given point instead of a real ligand's centre of
mass. No `-l`/ligand input needed.

```bash
python Interactions_search.py -r protein.pdb -c A --site-point 12.3 -4.5 30.1
python Interactions_search.py -r protein.pdb -c A --site-point 12.3 -4.5 30.1 --site-radius 8
python Interactions_search.py -r protein.pdb -c A --site-point 12.3 -4.5 30.1 --site-method ideal
```

Useful for generating GOLD bias points ahead of docking, before a ligand pose exists.
Output goes to `<receptor>_site_<x>_<y>_<z>/`: `<receptor>_site[_ideal].bpf`,
`<receptor>_site[_ideal]_bias.pdb` (dummy atoms, same `DON`/`ACC`/`ARO` convention as
the regular ligand `.bpf` — see [Bias Probe File](#bias-probe-file)), and a copy of the
receptor PDB.

| Argument | Description |
|---|---|
| `--site-point X Y Z` | Coordinate to search around. Requires `-r`; incompatible with `-x`/`-l`. |
| `--site-radius` | Search radius in Å (default: `10.0`). |
| `--site-method` | `atom` (default) or `ideal` — see below. |

Two independent point-placement methods, selected by `--site-method`:

- **`atom`** (default) — one bias point per receptor acceptor/donor atom or aromatic-ring
  centroid, placed **at the atom's own coordinate** (`receptor_site.py`, same logic used
  internally by the regular pipeline's active-site search). Works on any receptor PDB.
- **`ideal`** — a **fan** of points per group: instead of the atom's own position, computes
  where the ligand's *complementary* atom would ideally sit — H-bond geometry (distance +
  angle + dihedral) for acceptors/donors, plus stacked and parallel-displaced positions for
  aromatic rings (`src/interactions_search/ideal_sites.py`, a Python 3 port of the
  standalone `ideal_interaction_sites.py` script at the repo root, not otherwise wired into
  the package). Produces several points per acceptor/donor group (e.g. 5 around a
  carbonyl/carboxylate oxygen at 120°/150°/180°/210°/240°) and 14 per aromatic ring, instead
  of a single point. **Requires the receptor PDB to have explicit hydrogens** (Maestro/Amber
  naming — `HNE`/`HH11`/`HH12`/`HH21`/`HH22` for ARG, `HD21`/`HD22` ASN, `HE21`/`HE22` GLN,
  `HG`/`HG1`/`HH` SER/THR/TYR, `HE1` TRP, `HZ1`/`HZ2`/`HZ3` LYS, `HD1`/`HE2` HIS with explicit
  `HIE`/`HID`/`HIP` protonation) for the receptor-donor groups; acceptor groups based on
  heavy atoms only (carbonyl, carboxylate, amide, imidazole) and aromatic rings don't need H
  and still work without one — points needing a missing H atom are silently skipped rather
  than erroring.

This is receptor-only in both cases — it never touches the ligand side: the regular
per-pair `.bpf` (`options.bias`, see [Bias Probe File](#bias-probe-file)) keeps using the
ligand's real hot-point positions, where an "ideal" position doesn't apply since the actual
atom position is already known.

If the complex PDB has no `HETATM` records at all (ligand saved as `ATOM`, e.g. some CHARMM/AMBER-prepped structures), use `-f` to tell `split_pdb` which residue name(s) to treat as ligand:

```bash
# Ligand TF3 is stored as ATOM records, not HETATM
python Interactions_search.py -x complex.pdb -c A -f TF3

# Multiple forced ligand names; still use -n to pick one for this run
python Interactions_search.py -x complex.pdb -c A -f TF3 7FW -n TF3
```

### Mode 5 — Probe (simulate interactions at given coordinates)

Given a receptor and one or more coordinates, simulates which interactions a ligand
atom/group would make if it sat at each point — no ligand PDB needed. Each point is
tested in one or more **roles** (`--probe-type`, default `all`) against the complementary
receptor groups, with the same thresholds as the regular pipeline:

| Probe role | Receptor partner | Output `Type` | Validation |
|---|---|---|---|
| `acceptor` | donors (YAML `donors`) | `acceptor` | `Dist < Distances_Hidrogen_Bonds`; if the donor is an explicit H, also the D-H···probe angle (at the H) within `Angle_Hidrogen_Bonds_Min/Max` |
| `donor` | acceptors (YAML `acceptors`) | `donor` | distance + probe-Acceptor-Antecedent angle, same as the regular pipeline |
| `aromatic` | TYR/PHE/TRP ring centroids; ARG/LYS/HIS cations | `aromatic`, `pi_cation` | distance only (`Distances_Aromatic`, default 5.5 Å; `Distances_Pi_Cation`, default 5.0 Å) |
| `hydrophobic` | apolar C atoms, collapsed per residue | `hydrophobic` | `Distances_Hidrofobica` |
| `cation` | ASP/GLU carboxylate O; aromatic ring centroids | `salt_bridge`, `pi_cation` | `Distances_Salt_Bridge` / `Distances_Pi_Cation` (defaults 4.0 / 5.0 Å) |
| `anion` | ARG/LYS/HIP N | `salt_bridge` | `Distances_Salt_Bridge` (default 4.0 Å) |

```bash
# one point, all roles
python Interactions_search.py -r protein.pdb -c A --probe 12.3 -4.5 30.1

# several points, only as H-bond donor/acceptor
python Interactions_search.py -r protein.pdb -c A --probe 12.3 -4.5 30.1 --probe 10 -2 28 \
    --probe-type donor acceptor

# points from a file
python Interactions_search.py -r protein.pdb -c A --probe-file points.csv
python Interactions_search.py -r protein.pdb -c A --probe-file protein_site_ideal.bpf
```

`--probe-file` accepts: a `.bpf` (e.g. the Mode 4 output — its `don`/`acc`/`aro` column sets
each point's role), a `.pdb` (one point per `ATOM`/`HETATM`; resname `DON`/`ACC`/`ARO`/`HPH`/`CAT`/`ANI`
sets the role), or any text/CSV with `x y z [role]` rows (header lines are skipped). Points
without their own role use `--probe-type`. `--probe` and `--probe-file` can be combined.

Notes:
- A point can't define a ring plane, so for `aromatic`/`pi_cation` rows `Angle` is the angle
  between the receptor ring normal and the centroid→probe vector (0° = probe over the ring
  face, 90° = in the ring plane), reported for inspection but not used for validation.
- For acceptor probes, a receptor heavy-atom donor row is checked by distance only.
  A receptor explicit-H donor row also checks D-H···probe at H. This differs from
  the ligand acceptor-side angle, which uses the ligand's antecedent atom that a
  bare point doesn't have. An explicit-H row without a usable parent is rejected.
- Donor probes use the same donor–acceptor–antecedent angle as ligand donors.
  A missing antecedent produces `Angle = NaN`, `Interaction = No` in both modes;
  missing/degenerate required geometry is not accepted by distance alone.
- H-bond candidates use `max(Hydrogen_Bond_Search_Distance, Distances_Hidrogen_Bonds)`;
  the old fixed 4 Å prefilter no longer truncates a larger configured final cutoff.
- Every point is also checked for **steric clashes**: receptor heavy atoms within
  `Probe_Clash_Distance` (default 2.5 Å) are
  reported as `Type = clash`, `Interaction = Clash` and flagged in the console summary — a
  probe there would sit inside the protein. Mode 4 `ideal` points are placed at the H···A
  distance (1.9 Å, i.e. where the ligand's *H* would be), so treated as heavy-atom probes they
  clash with their own partner by construction.
- The receptor site is built once around the centroid of all points, with radius
  `centroid_distance` + the farthest point's distance to that centroid.

Output goes to `<receptor>_probe_<label>/` (`<label>` = the point's coordinates for a single
`--probe`, the file stem for `--probe-file`, or `<N>pts`):

| File | Content |
|---|---|
| `Probe_<rec>_all.csv` | Every candidate contact, including `No` and `Clash` rows |
| `Probe_<rec>_true.csv` | Only `Interaction == 'Yes'` |
| `<rec>_probe_points.pdb` | Probe points as dummy atoms (resname by role, `PRB` if tested with several), resid = probe number |
| `vmd_probe_<rec>.tcl` | Receptor + interacting residues + probe spheres + dashed lines per validated interaction (if `vmd_output: 'Yes'`) |
| `<rec>.pdb` | Copy of the receptor PDB (loaded by the `.tcl`) |

CSV columns: `Probe` (point number), `Probe_Type` (role tested, `-` for clashes),
`Probe_X/Y/Z`, `Pos R`, `Res`, `Atom`, `Dist`, `Type`, `Angle`, `Interaction`
(`Yes`/`No`/`Clash`), `Rec_X/Y/Z` (receptor atom, ring centroid, or mean of the collapsed
hydrophobic atoms), followed by `Reason` (decision code).

### Mode 6 — Pockets from MD hotspots

Builds the pocket around each group of hotspots from a cosolvent MD analysis (acceptor /
donor / hydrophobic clusters) and fills it with an AutoDock-style grid annotated with
properties.

```bash
python Interactions_search.py -r 3mss_complex.pdb -c B --hotspots results_global/ --exclude-res STI MS7
```

`--hotspots` points to a directory with `acceptors/`, `donors/` and/or `hydrophobics/`
subfolders, each with `clusters.csv` (`ws_id, x, y, z, R90_A, DG, occ_prob, ...`) and,
optionally, `cluster_points.pdb` (raw probe positions, resid = `ws_id`) and a `grid_*.dx`
ΔG map. The receptor must be in the same frame as the hotspots. `--exclude-res` removes
residues from the receptor — needed when ligands are stored as `ATOM` records, otherwise
they fill the pocket and every grid point there is discarded as a clash.

Per site:

1. **Sites** — hotspot centers are grouped by single linkage (`link_distance`, 8 Å); sites
   are numbered by summed ΔG (most favourable = 1). Sites with fewer than `min_hotspots`
   are listed but not built. Each hotspot contributes its `cluster_points.pdb` points within
   `R90` of its center (just the center if there are none).
2. **Pocket residues** — every residue with a heavy atom within `residue_cutoff` (4 Å) of
   any hotspot point of the site.
3. **Grid** (`grid_spacing`, 0.375 Å) over the box enclosing those residues and points. A
   point is kept if it is inside the convex hull of the pocket residues' heavy atoms, more
   than `grid_clash` (2.6 Å) from every receptor heavy atom, buried (at least
   `grid_buriedness` = 40% of 30 rays hit the receptor within 10 Å, LIGSITE-style), and
   connected (26-neighbourhood) to a hotspot point.
4. **Annotation** of each grid point: ΔG from each `.dx` map at that point, `Best_Type`
   (type with the lowest ΔG if it reaches `grid_dg_threshold`, -1 kcal/mol, else `none`),
   and the receptor environment: number of receptor donors / acceptors (within
   `Distances_Hidrogen_Bonds`), hydrophobic atoms (`Distancia_Hidrofobica`), aromatic ring
   centroids (`Distances_Aromatic`), cations and anions (`Distances_Salt_Bridge`, default
   4.0 Å) — the same partners and
   thresholds as the probe mode (Mode 5), distance only.

Default thresholds were calibrated on 3MSS: the grid of the ATP site contains 100% of the
STI (imatinib) heavy atoms; for the myristate site it contains the buried part of MS7 and
leaves out its solvent-exposed head. The parameters live in the `hotspot_pocket:` section
of `Interacciones_variables.yml`.

Output goes to `<receptor>_hotspot_pockets/`:

| File | Content |
|---|---|
| `sites_summary.csv` | One row per site: hotspot counts by type, ΣΔG, best ΔG, center, residues, grid points, `Volume_A3` (points × spacing³), fraction of grid points per `Best_Type`, `Built` |
| `site_<n>/hotspots.csv` | Hotspots of the site |
| `site_<n>/residues.csv` | Pocket residues: `Pos`, `Residue`, `Min_Dist` to a hotspot point, `N_Atoms` within the cutoff, `Hotspots`, `Types` |
| `site_<n>/grid.csv` | Grid points: `X/Y/Z`, `Buriedness`, `DG_acceptor/donor/hydrophobic`, `Best_Type`, `Best_DG`, `N_Rec_Donors/Acceptors/Hydrophobic/Aromatic/Cations/Anions` |
| `site_<n>/grid.pdb` | Grid points as dummy atoms: resname `ACC`/`DON`/`HPH`/`NON` (Best_Type), occupancy = buriedness, B-factor = Best_DG |
| `site_<n>/hotspots.pdb` | Hotspot centers: resname by type, occupancy = R90, B-factor = ΔG |
| `site_<n>/pocket_mask.dx` | Pocket mask (1 inside, 0 outside); an isosurface at 0.5 shows the pocket shape in VMD/PyMOL/Chimera |
| `site_<n>/vmd_site_<n>.tcl` | Receptor + pocket residues + grid points coloured by type + hotspots + mask isosurface |

Notes:
- The `.dx` maps are sparse at voxel level: most of the pocket volume has ΔG ≈ 0 (no
  cosolvent preference), so most points are `Best_Type = none`; the typed points are the
  cores of the hotspots. The `N_Rec_*` columns describe the chemistry everywhere else.
- Grid points outside a `.dx` box get `NaN` in that map's column.

### Arguments

| Argument | Description |
|---|---|
| `-x / --complex` | Complex PDB. Alternative to `-r`. |
| `-r / --receptor_pdb` | Receptor PDB (already separated). |
| `-l / --ligand_input` | Ligand PDB(s). Accepts one or several (batch). |
| `-c / --chain_receptor` | Protein chain (e.g. `A`). |
| `-n / --lig_name` | HETATM name when multiple groups exist in `--complex`. |
| `-f / --force_ligand` | Residue name(s) to treat as ligand even if stored as `ATOM` instead of `HETATM` in `--complex`. |
| `--config` | Path to the YAML config (default: `Interacciones_variables.yml` at the repo root). |
| `--probe X Y Z` | Probe mode (Mode 5). Repeatable. Requires `-r`; incompatible with `-x`/`-l`/`--site-point`. |
| `--probe-file` | Probe points from a `.bpf`, `.pdb` or `x y z [role]` text/CSV file (Mode 5). |
| `--hotspots DIR` | Hotspot mode (Mode 6). Requires `-r`; incompatible with `-x`/`-l`/`--site-point`/`--probe`. |
| `--exclude-res` | Resnames removed from the receptor in hotspot mode (e.g. ligands stored as `ATOM`). |
| `--probe-type` | Role(s) for points without their own: `all` (default), `acceptor`, `donor`, `aromatic`, `hydrophobic`, `cation`, `anion`. |

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
    └── Aromatic rings:    detected via ring_info, filtered by size >= 5 and by
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
    ├── Salt bridge:   ± group lig ↔ ∓ group rec       (Distances_Salt_Bridge; default 4.0 Å)
    └── π-cation:      lig ring ↔ ARG/LYS/HIS rec      (Distances_Pi_Cation; default 5.0 Å)
        │
        ▼
[7] Angle validation
    ├── H-bond:    D-A···Antecedent angle between 100° and 200°
    ├── Aromatic:  angle between ring planes
    │               angle < Aromatic_Parallel_Max (default 30°) → parallel / sandwich
    │               angle > Aromatic_TShaped_Min (default 60°, up to 90°) → T-shaped
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
  Hydrogen_Bond_Search_Distance: 4.0  # candidate radius, expanded to final H-bond cutoff if needed
  Distances_Salt_Bridge:   4.0   # Å — shared by ligand/probe and hotspot environment counts
  Distances_Pi_Cation:     5.0   # Å — shared by ligand and probe
  Probe_Clash_Distance:    2.5   # Å — probe clashes, independent of hotspot grid_clash
  centroid_distance:       12.0   # Å — active site search radius
  Distances_C_Simple:       1.54  # Å — C-C single bond (reference)
  Distances_C_Doble:        2.56  # Å — C=C double bond (reference)

angulos:
  Angle_Hidrogen_Bonds_Min: 100    # ° — minimum Donor-Acceptor-Antecedent angle
  Angle_Hidrogen_Bonds_Max: 180    # ° — maximum angle (180° is the geometric ceiling)
  Aromatic_Parallel_Max: 30       # ° — accept inter-plane angle strictly below this
  Aromatic_TShaped_Min:  60       # ° — accept inter-plane angle strictly above this

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

hotspot_pocket:        # Mode 6 (--hotspots)
  link_distance:     8.0     # Å — single-linkage distance between hotspot centers
  min_hotspots:      3       # sites with fewer hotspots are not built
  residue_cutoff:    4.0     # Å — residue heavy atom to hotspot point
  grid_spacing:      0.375   # Å
  grid_clash:        2.6     # Å — min distance from grid point to receptor heavy atom
  grid_buriedness:   0.4     # 0-1 — min fraction of buried rays
  grid_dg_threshold: -1.0    # kcal/mol — min ΔG for Best_Type != 'none'

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

## Side-chain chi angles

`chi_angles.py` computes side-chain torsion angles (chi1-chi5, in degrees, standard IUPAC
dihedral convention around the B-C bond of each `A-B-C-D` atom quartet), using the atom
definitions in [`chi_angles.json`](src/interactions_search/data/chi_angles.json) (e.g. LEU
chi1 = `N-CA-CB-CG`, chi2 = `CA-CB-CG-CD1`). Coordinates are taken from `DF_Active_Site`,
which already holds every atom of each active-site residue (not just the H-bond/aromatic
points of interest), so no extra PDB read is needed. A chi is left blank (`None`) when the
residue has no such angle (e.g. ALA has none, VAL only has chi1) or when one of its 4 atoms
is missing from the PDB (common for partially resolved side chains in crystal structures).
Histidine protonation variants (`HID`/`HIE`/`HIP`) reuse the `HIS` atom definitions.

Two functions cover different scopes:

- **`compute_pocket_chi_angles()`** — one row per (pocket, residue) for every hydrophobic
  pocket *candidate* fragment found by `search_hydrophobic_pockets()`, whether or not it
  qualifies as `Is_Pocket == Yes` (that column is included in the output so rows can be
  filtered afterwards). A fragment touched by a single residue, for example, never clears
  the `Coverage_R` bar on its own but its residue's rotamer is still reported.
- **`compute_active_site_chi_angles()`** — one row per residue for **every** residue in
  `DF_Active_Site`, i.e. anything within `centroid_distance` Å of the ligand's centre of mass
  (see `active_site_residues()` in `receptor_site.py`), regardless of whether it forms any validated interaction or
  belongs to a hydrophobic pocket at all. This is the broadest view: rotamers for the whole
  analysed neighbourhood around the ligand.

`chi_angles.json` lives inside the package (`src/interactions_search/data/`), not at the
repo root like `Interacciones_variables.yml`: it's static IUPAC reference data the user
never edits per-project, not tunable config, so `chi_angles.py` loads it as a packaged
resource via `importlib.resources` (declared in `pyproject.toml` under
`[tool.setuptools.package-data]`). This means it ships correctly with a `pip install` from
a wheel, unlike a file that only exists in a cloned repo.

`plot_chi_profile()` (`src/interactions_search/plotting.py`) renders one PNG per chi
(chi1-chi5) from the `ActiveSite_<rec>_<lig>_chi.csv` data: X axis = active-site residues
(labelled `<residue><pos>`), Y axis = that chi in degrees (-180°, 180°], same scale as the
Ramachandran plot. Residues without that particular chi are skipped, and no PNG is written
for a chi that has zero values across the whole active site (e.g. `chi3.png` is only
generated if at least one active-site residue — ARG, GLN, GLU, LYS, MET... — has a chi3).
Generated only if `options.volume_plot: 'Yes'`.

Each point is colored by the residue's physicochemical class (`_AA_CLASS_COLOR` in
`plotting.py`): nonpolar/hydrophobic (blue), polar uncharged (orange), or charged —
acidic/basic together (aqua/green), with a legend. Only 3 classes, not the finer
5-way split (nonpolar / aromatic / polar / acidic / basic) some textbooks use: a
scatter plot compares every pair of colors at once, and past 3 categories no ordering
of the project's validated categorical palette clears the normal-vision separation
floor (see `dataviz`'s `references/palette.md`) — cutting to 3 keeps every pair
legible instead of shipping colors two people in the room can't reliably tell apart.
The X-axis label already gives the exact residue, so color only needs to carry the
coarser class.

---

## Ramachandran (backbone phi/psi) angles

`compute_active_site_phi_psi()` (`src/interactions_search/ramachandran.py`) computes the
backbone torsion angles phi and psi (degrees, IUPAC convention) for every residue in
`DF_Active_Site` — same scope as `compute_active_site_chi_angles()`, i.e. every residue
within `centroid_distance` Å of the ligand, not just ones forming a validated interaction.

Unlike chi (which only needs atoms from the residue itself), phi and psi need the
**neighbouring** residues' backbone atoms — phi = C(i-1)-N(i)-CA(i)-C(i), psi =
N(i)-CA(i)-C(i)-N(i+1) — so the previous/next residue is looked up in the full receptor
chain (`structure[0][chain_receptor]`), not in `DF_Active_Site`, in case that neighbour
falls just outside the active-site radius. A residue's phi and/or psi is left as `None`
when:
- it's the first/last residue of the chain (no previous/next residue at all),
- there's a numbering gap (the neighbour found isn't `pos ± 1`, i.e. a missing residue in
  the crystal structure — no real peptide bond to measure), or
- one of the four required backbone atoms (`N`/`CA`/`C`) is missing from the PDB.

`plot_ramachandran()` (`src/interactions_search/plotting.py`) renders a phi-vs-psi scatter
PNG (rows with `None` dropped), labelling each point with `<residue><pos>` and marking
GLY (▲, no side-chain steric restriction — often falls outside the usual favoured regions)
and PRO (■, ring-constrained phi) separately from the rest (●), since those are the
expected outliers of a standard Ramachandran plot. Generated only if `options.volume_plot:
'Yes'` (same flag that gates the hull PNGs).

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
| `<folder>/Interaction_<rec>_<lig>_threshold.csv` | Filtered by each interaction type's own distance cutoff; angles are not required here |
| `<folder>/Interaction_<rec>_<lig>_true.csv` | Validated by distance and angle |
| `<folder>/Pockets_<rec>_<lig>.csv` | Hydrophobic pocket candidates (see "Hydrophobic Pockets" above), one row per ligand fragment |
| `<folder>/Pockets_<rec>_<lig>_chi.csv` | Side-chain chi angles (chi1-chi5, °) of the residues in every hydrophobic pocket candidate fragment (validated or not, tagged by `Is_Pocket`), one row per (pocket, residue) — see "Side-chain chi angles" below |
| `<folder>/ActiveSite_<rec>_<lig>_chi.csv` | Side-chain chi angles (chi1-chi5, °) of **every** residue in the active site, regardless of interaction/pocket status, one row per residue — see "Side-chain chi angles" below |
| `<folder>/ActiveSite_<rec>_<lig>_chi<N>.png` | Scatter plot of chi`<N>` (`N` = 1-5) across the active-site residues that have it, one PNG per chi that has at least one value (if `volume_plot: 'Yes'`) |
| `<folder>/ActiveSite_<rec>_<lig>_ramachandran.csv` | Backbone phi/psi angles (°) of **every** residue in the active site, one row per residue — see "Ramachandran (backbone phi/psi) angles" below |
| `<folder>/ActiveSite_<rec>_<lig>_ramachandran.png` | Ramachandran scatter plot (phi vs psi) of the active-site residues (if `volume_plot: 'Yes'`) |
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
| `Reason` | Decision code; multiple failed criteria are separated by `;` |

### Decision reasons and run records

Ligand and probe interaction CSVs explain each reported candidate using these codes:

| Code | Meaning |
|---|---|
| `distance_outside_cutoff` | Distance does not satisfy the strict cutoff |
| `angle_outside_range` | Angle does not satisfy the configured range |
| `missing_required_geometry` | Required geometric information is unavailable |
| `distance_and_angle_pass` | Both required criteria pass |
| `distance_pass` | A distance-only contact passes |
| `distance_only_no_probe_orientation` | Distance passes for an unoriented probe |
| `steric_clash` | Probe clashes with the receptor |

Reasons use unrounded geometry. Only candidates found by each detector are reported;
an absent pair is not an explicit rejection row.

Every ligand, probe, site-bias and hotspot analysis writes `config_used.yml` and
`run_metadata.json` in its output folder. Metadata includes effective parameters,
CLI arguments when available, input SHA-256 hashes before and after processing,
software versions, source fingerprint, timestamps, completion status and hashes of
new or modified outputs inside that folder. External cumulative CSVs are excluded.

`run_history/<run_id>/` preserves each configuration and metadata record, plus
`ligand_reference.sdf` when a chemical reference was supplied. The saved SDF retains
template atom order. Root-level records describe the latest attempt. Failed attempts
are marked `failed` and do not claim unchanged outputs from earlier runs. Handled
interruptions are marked `interrupted`; a forcibly killed process can remain `running`.
CLI validation errors before analysis starts do not create a run record.

To repeat an analysis, retain the original inputs and environment, use the saved
YAML with `--config`, restore mode-specific arguments from metadata and, when present,
use the archived SDF with `--ligand-sdf`. History does **not** archive all input or
result files: repeated analyses still overwrite ordinary outputs. See
[SOP section 13](docs/SOP.md#13-motivos-y-registro-automático-de-corridas).

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
├── config_used.yml           ← effective configuration, latest attempt
├── run_metadata.json         ← provenance and status, latest attempt
├── run_history/<run_id>/     ← configuration, metadata and optional reference SDF
├── <receptor>.pdb             ← copy of the receptor PDB
├── <ligand>.pdb               ← copy of the ligand PDB
├── <ligand>_old.pdb           ← pre-cleanup copy (remove_bias)
├── Interaction_*_all.csv      ← all interactions, no filter
├── Interaction_*_threshold.csv← filtered by distance
├── Interaction_*_true.csv     ← validated by distance + angle
├── Pockets_*.csv              ← hydrophobic pocket candidates
├── Pockets_*_chi.csv          ← chi angles of pocket-fragment residues (validated or not)
├── ActiveSite_*_chi.csv       ← chi angles of every active-site residue
├── ActiveSite_*_chi<N>.png    ← scatter of chi<N> across active-site residues (if volume_plot: Yes, one per chi with data)
├── ActiveSite_*_ramachandran.csv ← phi/psi angles of every active-site residue
├── ActiveSite_*_ramachandran.png ← Ramachandran scatter plot (if volume_plot: Yes)
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

The suite also includes offline crystallographic validation against ABL–imatinib
(1IEP) and streptavidin–biotin (1STP): 16 named contact checks, independent geometry,
rigid transforms, displaced ligands and incomplete receptor rings. Protonation
tests cover neutral/charged groups with and without explicit H. The complete
suite currently has **170 passing tests**. See
[`docs/VALIDACION_QUIMICA.md`](docs/VALIDACION_QUIMICA.md) for provenance, tolerances,
corrections and limitations. The distributed YAML now includes SER OG, THR OG1
and TYR OH donor atoms for crystals without explicit H; custom YAML tables remain
authoritative. Incomplete receptor rings are skipped instead of reconstructed.

---

## Notes

- The script should be run from the directory containing the PDB files, or use absolute paths.
- For batch analysis of multiple ligands, the script can be called in a shell loop; with `cumulative_output: 'Yes'`, `Interactions_close.csv` and `CM_all.csv` are appended automatically across runs (set to `'No'` to disable).
- Non-standard residues not listed in `acceptors` / `donors` in the YAML are silently skipped.
- Ligand rings must contain at least 5 atoms and be planar (RMSD to the best-fit plane ≤ `Ring_Planarity_RMSD_Max`). With a chemical reference, aromaticity is also required. PDB-only mode retains the geometric approximation; `--legacy-rings` restores the previous size filter (> 5 atoms).
- `geometry.dihedral_angle()` (used for chi and phi/psi) had a sign-convention bug in earlier versions — its output was the exact negative of the standard IUPAC/Biopython/PyMOL convention. Fixed and verified against `Bio.PDB.vectors.calc_dihedral()` across a full chain. If you have chi/phi/psi CSVs generated before this fix, their angles are sign-flipped relative to the current output.

---

## See also

- [`docs/SOP.md`](docs/SOP.md) — procedimiento operativo: PDB solo, referencias químicas opcionales, compatibilidad, errores y validación.
- [`docs/VALIDACION_QUIMICA.md`](docs/VALIDACION_QUIMICA.md) — validación cristalográfica y química, casos de control y hallazgos de la etapa 3.
- [`docs/BASELINE.md`](docs/BASELINE.md) — full behavioral snapshot: pipeline details, config keys, CSV schema, regression check command
- [`Interacciones_variables.yml`](Interacciones_variables.yml) — live configuration file with all thresholds and per-residue donor/acceptor tables
