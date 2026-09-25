# Progreso: modo sondeo + pockets desde hotspots de dinámica

Estado al 2026-09-24. Todo esto está **sin commitear** en `main` (working tree), sobre el commit
`5f619ce`. La documentación de usuario (en inglés) está en el `README.md`, secciones **Mode 5** y
**Mode 6**. Este archivo es la bitácora para retomar el trabajo.

---

## 1. Qué se hizo

### 1.1 Modo sondeo (`--probe`, "Mode 5")

Simula qué interacciones haría un átomo o grupo del ligando si estuviera en una coordenada dada,
sin necesitar un ligando.

```bash
python Interactions_search.py -r prot.pdb -c A --probe X Y Z [--probe X Y Z ...] [--probe-type donor acceptor ...]
python Interactions_search.py -r prot.pdb -c A --probe-file puntos.bpf|.pdb|.csv
```

- **Roles** (`PROBE_TYPES`): `acceptor`, `donor`, `aromatic`, `hydrophobic`, `cation`, `anion`.
  El default es `all`. Umbrales: los mismos del YAML y de `contacts.py` que usa `analyze_pair()`.
- **Ángulos:**
  - Sonda donor: usa el ángulo sonda–Aceptor–Antecedente, igual que el pipeline.
  - Sonda aceptor frente a un H explícito del receptor: usa el ángulo D‑H···sonda en el H. Es más
    estricto que el pipeline, que necesita el antecedente del ligando, y una sonda no lo tiene.
  - Donor pesado sin H: se valida solo por distancia.
- **Aromático:** un punto no define un plano. `Angle` es el ángulo entre la normal del anillo del
  receptor y el vector centroide→sonda; se reporta pero no se usa para validar.
- **Choques:** átomos pesados a menos de 2.5 Å → `Type=clash`, `Interaction=Clash`
  (constante `_CLASH_DIST`).
- **Tipo desde el archivo:** en un `.bpf` la columna `don/acc/aro` fija el rol. En un `.pdb` lo fija
  el resname `DON/ACC/ARO/HPH/CAT/ANI`. En CSV/texto, la última columna.
- **Salida** en `<rec>_probe_<label>/`: `Probe_<rec>_all.csv`, `Probe_<rec>_true.csv`,
  `<rec>_probe_points.pdb`, `vmd_probe_<rec>.tcl`.

**Validación** con el FAD de `LucAI/ankros_holo_fmn_noDNA.pdb` (cadena A, protonado): sondas en
6 átomos del FAD con su rol real. Coincide con `analyze_pair` salvo:
- O1A–THR231: el H apunta hacia otro lado (D‑H···A = 63°). El sondeo lo rechaza; el pipeline lo
  acepta por usar otro ángulo.
- ALA334: queda a 3.85 Å de C7M y el sondeo la detecta. El pipeline normal la pierde porque su
  centro de masa está a unos 12.2 Å del centro del FAD, afuera de `centroid_distance` = 12.
  → **Limitación del pipeline en ligandos largos.**

Sondear los puntos `--site-method ideal` marca choque en casi todos, y es real por dos motivos:
esos puntos están a 1.9 Å del compañero (posición del **H** del ligando, no del átomo pesado), y
además se generan para residuos enterrados.

### 1.2 Pockets a partir de hotspots de dinámica (`--hotspots`, "Mode 6")

A partir de los clusters de hotspots de una dinámica con cosolvente (EOH), agrupa en sitios, saca
los residuos a 4 Å y arma una grilla estilo AutoDock anotada con propiedades.

```bash
python Interactions_search.py -r Inputs/3mss_complex.pdb -c B \
    --hotspots Inputs/results_global --exclude-res STI MS7
```

**Formato de entrada** (ejemplo en `Inputs/results_global/`): subcarpetas `acceptors/`, `donors/`
y `hydrophobics/`, cada una con:
- `clusters.csv`: `ws_id, x, y, z, R90_A, DG, occ_prob, ...`. `results.csv` es idéntico y no se usa.
- `cluster_points.pdb`: posiciones crudas del EOH; resid = `ws_id`.
- `grid_EOH_<tipo>.dx`: ΔG por vóxel, 0.375 Å, 100³, valores ≤ 0.

**Algoritmo** (`hotspot_pocket.py`):
1. **Sitios:** single-linkage de los centros con `link_distance`. Los sitios se ordenan por ΣΔG.
   Los que tienen menos de `min_hotspots` se listan pero no se construyen.
2. **Puntos por hotspot:** los de `cluster_points.pdb` dentro de `R90` del centro.
3. **Residuos:** los que tienen un átomo pesado a ≤ `residue_cutoff` de algún punto.
4. **Grilla** con `grid_spacing`. Un punto queda si cumple todo esto:
   - está dentro de la envolvente convexa de los átomos pesados de esos residuos,
   - está a más de `grid_clash` de todo átomo pesado,
   - tiene enterramiento ≥ `grid_buriedness` (30 rayos, pasos de 1.5 Å hasta 10 Å, choque a < 2 Å),
   - está conectado (vecindad 26) a algún punto de hotspot a ≤ 2 Å.
5. **Anotación:**
   - `DG_<tipo>` del `.dx`, tomando el vóxel más cercano.
   - `Best_Type`: el tipo de menor ΔG si llega a ≤ `grid_dg_threshold`; si no, `none`.
   - `N_Rec_*`: cantidad de donores, aceptores, átomos hidrofóbicos, anillos, cationes y aniones
     del receptor a distancia de interacción.

**Parámetros** en la sección `hotspot_pocket:` del YAML (pydantic `HotspotPocket` en `config.py`):

| Parámetro | Default | Cómo se eligió |
|---|---|---|
| `link_distance` | 8.0 Å | Con 6–7 Å el sitio de STI se parte en 2; con 8 Å STI y MS7 quedan cada uno en un sitio propio y separados. |
| `min_hotspots` | 3 | |
| `residue_cutoff` | 4.0 Å | Lo pidió el usuario. |
| `grid_spacing` | 0.375 Å | Mismo espaciado que AutoDock y que los `.dx`. |
| `grid_clash` | 2.6 Å | El átomo de STI más cercano a la proteína está a 2.72 Å. |
| `grid_buriedness` | 0.4 | Puntos en solvente: mediana 0.17, p99 0.37. STI: mínimo 0.67. MS7: p10 = 0.27, porque tiene partes expuestas. |
| `grid_dg_threshold` | −1.0 kcal/mol | |

**Salida** en `<rec>_hotspot_pockets/`: `sites_summary.csv`, y en `site_<n>/`: `hotspots.csv`,
`residues.csv`, `grid.csv`, `grid.pdb` (resname = Best_Type, occupancy = enterramiento,
B = Best_DG), `hotspots.pdb`, `pocket_mask.dx` (isosuperficie a 0.5 = forma del pocket),
`vmd_site_<n>.tcl` y una copia del receptor.

**Resultado con 3MSS** (79 hotspots → 7 sitios, 5 construidos, ~15 s):

| Sitio | Hotspots | Residuos | Volumen | Qué es |
|---|---|---|---|---|
| 1 | 32 (12 acc, 13 don, 7 hyd), ΣΔG −86.7 | 45 | 997 Å³ | Bolsillo de ATP / STI: THR315, GLU286, MET318, ASP381 (DFG) |
| 2 | 31 (11 acc, 9 don, 11 hyd), ΣΔG −79.6 | 42 | 1035 Å³ | Bolsillo del miristato / MS7 |
| 3 | 6 | 16 | 221 Å³ | |
| 4 | 4 | 8 | 37 Å³ | |
| 5 | 3 | 10 | 13 Å³ | |

**Validación:** la grilla del sitio 1 contiene el 100% de los átomos pesados de STI (máx. 0.29 Å
al punto más cercano). El sitio 2 contiene toda la parte enterrada de MS7 (C9–C32). Quedan afuera
N1, C4, C34, O35, N36 y C38: la cabeza expuesta al solvente, a 4–8 Å de la proteína, con
enterramiento de 0.13–0.37.

**Limitaciones vistas:**
- El 95–97% de los puntos queda `Best_Type = none`. Los `.dx` son dispersos (mediana ΔG ≈ 0) y
  solo los centros de los hotspots bajan de −1. El resto lo describen las columnas `N_Rec_*`.
- La caja de los `.dx` (origin −24.5, −20.1, −19.7; 37.5 Å de lado) no cubre el 20% de la grilla
  del sitio 2, que queda con `DG_* = NaN`.
- En `3mss_complex.pdb` los ligandos STI y MS7 están guardados como `ATOM`; hay que excluirlos
  con `--exclude-res`.

---

## 2. Archivos tocados

| Archivo | Cambio |
|---|---|
| `src/interactions_search/probe.py` | **nuevo**: `probe_interactions`, `read_probe_file`, `write_probe_pdb`, `PROBE_TYPES` |
| `src/interactions_search/hotspot_pocket.py` | **nuevo**: `load_hotspots`, `group_sites`, `read_dx`/`write_dx`, `analyze_hotspot_pockets` |
| `src/interactions_search/pipeline.py` | `analyze_probe()` + import de `scripting_vmd_probe` |
| `src/interactions_search/vmd.py` | `scripting_vmd_probe()` |
| `src/interactions_search/cli.py` | opciones `--probe`, `--probe-file`, `--probe-type`, `--hotspots`, `--exclude-res`; el armado de `cfg` pasó a `_build_cfg(config_path)` |
| `src/interactions_search/config.py` | modelo `HotspotPocket` y campo `hotspot_pocket` en `InteractionConfig` |
| `Interacciones_variables.yml` | sección `hotspot_pocket:` |
| `README.md` | módulos nuevos en la tabla, Mode 5, Mode 6, argumentos nuevos, bloque YAML |
| `tests/test_probe.py`, `tests/test_hotspot_pocket.py` | **nuevos** (10 tests, en proceso, sin subprocess) |
| `Inputs/` | datos de ejemplo del usuario (sin seguimiento en git; decidir si se commitean o van a `.gitignore`) |

---

## 3. Entorno y cómo correr

- No existe `.venv/`, así que `tests/test_smoke.py` falla con `FileNotFoundError`: tiene
  hardcodeado `.venv/bin/python`. Los demás tests andan.
- El Python base de anaconda **importa otra copia del paquete**
  (`/media/.../Laboratorio/Interacciones/Interactions-search/src`). Para usar este repo:

```bash
PYTHONPATH=src python -m pytest -q tests/test_probe.py tests/test_hotspot_pocket.py tests/config
PYTHONPATH=src python Interactions_search.py ...
```

- ruff: los archivos nuevos pasan. El resto del paquete ya tenía unos 100 avisos (líneas largas,
  imports), que no se tocaron.

---

## 4. Próximos pasos posibles

- [ ] **`Best_Type` por hotspot:** asignarle a cada punto el tipo y ΔG del hotspot más cercano
      dentro de su R90, en vez del valor por vóxel del `.dx`. Bajaría mucho el 95% de `none`.
- [ ] **Caja de docking:** generar el `.gpf` de AutoDock por sitio (centro, `npts` pares,
      spacing 0.375) que envuelva la grilla.
- [ ] Revisar la grilla en VMD: `vmd -e vmd_site_1.tcl` dentro de `site_1/`.
- [ ] Decidir qué hacer con `Inputs/`: commitear como ejemplo o agregar a `.gitignore`.
- [ ] Commit de todo lo anterior (hoy está sin commitear).

### Pendientes detectados al revisar el código (no tocados)

- En el README:
  - La instalación "script-only" no funciona: el shim necesita el paquete y a la lista le faltan
    `pydantic` y `openbabel-wheel`.
  - El ejemplo de Python API usa `load_config`, que no está exportado.
  - Menciona `ideal_interaction_sites.py`, que no existe.
  - El diagrama del pipeline dice 100–200° (el código usa 180°) y el SMARTS `[N;H1]`, que ya se
    sacó.
  - El párrafo de `-f` está dentro del Mode 4.
- `io_pdb.remove_bias` modifica **in situ** el PDB de entrada del ligando.
- `geometry.get_aromatic_coord` tira `IndexError` si falta un átomo del anillo (cristales
  incompletos).
- Puentes salinos: los SMARTS dependen de cargas u órdenes de enlace que RDKit no percibe desde un
  PDB, así que casi nunca se detectan.
- π-catión trata a toda HIS como catión.
- `align.py`: `split('.')[0]` rompe con rutas como `./x.pdb`; habría que usar `Path.stem`.
- `config._DEFAULT_CONFIG_PATH` solo funciona con instalación editable.
- `analyze_pair`: los residuos cuyo centro de masa cae apenas afuera de `centroid_distance` se
  pierden aunque tengan contacto (caso ALA334 del FAD). Posible arreglo: seleccionar por distancia
  mínima átomo–ligando en vez de por centro de masa.
