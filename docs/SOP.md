# SOP — Análisis de interacciones proteína–ligando

Actualizado: 2026-09-25.

## 1. Objetivo y alcance

Analizar contactos proteína–ligando desde estructuras PDB y generar tablas y
visualizaciones reproducibles. Este procedimiento incluye la referencia química
opcional, la detección de anillos de cinco miembros, la opción de conservar el
filtro de anillos anterior y la unificación de criterios de interacción descrita
en la sección 11.

El uso con PDB solo sigue disponible. No se requiere proporcionar SDF ni SMILES.
Los modos sin ligando (site bias, probes y hotspots) se describen en el
[README](../README.md); no admiten las opciones de referencia química del ligando.
Las reglas de distancia compartidas sí se aplican a probes y a los recuentos de
entorno de hotspots según la sección 11.
[BASELINE.md](BASELINE.md) conserva una instantánea histórica y no sustituye este
procedimiento actualizado.

## 2. Preparación

Desde la raíz del repositorio, usar Python 3.11 o posterior con las dependencias
del proyecto instaladas:

```bash
python -m pip install -e '.[dev]'
```

Preparar el receptor y el ligando en el mismo sistema de coordenadas y seleccionar
la cadena del receptor. Alternativamente, proporcionar un PDB del complejo.
Los nombres de archivo y de cadena de los ejemplos deben reemplazarse por los
del sistema analizado.

Usar copias de trabajo de los PDB: la limpieza preexistente `remove_bias()` elimina
del archivo de ligando las líneas que contienen ` CM ` y guarda una copia
`<ligando>_old.pdb` en la carpeta de resultados. La aplicación de la referencia
química ocurre en memoria y no reescribe enlaces ni cargas del PDB de entrada.

Conservar los archivos de entrada originales y el entorno utilizado. Cada análisis
genera automáticamente `config_used.yml`, `run_metadata.json` y un historial por
corrida con la referencia química, si se utiliza (ver sección 13). El historial
no reemplaza una copia de los archivos de entrada ni de los resultados.
Puede seleccionarse otra configuración mediante `--config ruta/config.yml`.

## 3. Elegir el modo químico

| Entrada | Donores y aceptores | Anillos | Coordenadas |
|---|---|---|---|
| Solo PDB | Procedimiento existente de SMARTS y filtro de hidrógenos con Open Babel | Planaridad y tamaño mínimo de 5 átomos | PDB |
| PDB + SMILES o SDF | Definiciones de características químicas de RDKit sobre la molécula corregida | Aromaticidad asignada, planaridad y tamaño mínimo de 5 átomos | PDB |
| Solo PDB + `--legacy-rings` | Procedimiento existente | Planaridad y más de 5 átomos, como antes | PDB |

La planaridad usa `Ring_Planarity_RMSD_Max` del YAML. Sin referencia, se conserva
una aproximación geométrica: un anillo plano no queda demostrado químicamente
aromático por ese criterio. Con referencia, todos los átomos del anillo deben
estar marcados como aromáticos por RDKit.

## 4. Ejecutar con PDB solo

```bash
# Receptor y ligando separados
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A

# Varios ligandos, sin referencia química
python Interactions_search.py -r protein.pdb -l lig1.pdb lig2.pdb -c A

# Complejo: seleccionar el ligando por nombre de residuo
python Interactions_search.py -x complex.pdb -c A -n LIG

# Ligando guardado como registros ATOM
python Interactions_search.py -x complex.pdb -c A -f LIG -n LIG
```

Se mantienen las entradas y las salidas del flujo previo. La diferencia por
defecto es que ahora también se consideran anillos planos de cinco miembros,
por lo que pueden aparecer contactos adicionales.

## 5. Ejecutar con referencia química opcional

La referencia debe representar el mismo ligando, incluida la protonación deseada.
El siguiente SMILES es únicamente un ejemplo para acetamida; reemplazarlo por el
SMILES real del ligando:

```bash
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A --ligand-smiles 'CC(=O)N'
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A --ligand-sdf ligand.sdf
python Interactions_search.py -x complex.pdb -n LIG -c A --ligand-sdf ligand.sdf
```

Condiciones de uso:

- Usar solo una opción: `--ligand-smiles` o `--ligand-sdf`.
- El SDF debe contener exactamente una molécula válida.
- Analizar un solo ligando por ejecución cuando se proporciona referencia. Para
  ligandos diferentes, ejecutar cada uno con su referencia correspondiente.
- Debe coincidir el grafo completo de átomos pesados: elementos y conectividad.
  La referencia corrige órdenes de enlace; no repara una conectividad PDB errónea.
- Los hidrógenos explícitos del PDB deben ser compatibles con la referencia.

El programa transfiere órdenes de enlace, cargas formales, aromaticidad y estado
de hidrógenos al objeto molecular usado por el análisis. Conserva los átomos,
su orden, nombres, seriales y coordenadas del PDB; no agrega hidrógenos con nuevas
coordenadas ni usa la conformación del SDF. La referencia no valida la
estereoquímica de la pose.

Los donores y aceptores se identifican con SMARTS HBA/HBD de `rdkit.Chem.Lipinski`
sobre una copia con H colapsados, añadiendo los N aromáticos protonados con H.
La molécula original conserva sus hidrógenos e índices. Esta corrección de la
etapa 3 reemplaza `BaseFeatures.fdef`, que incluía donores potencialmente
ionizables sin H en el estado suministrado. El filtro de hidrógenos inferidos
con Open Babel no reemplaza la química de referencia. La molécula
corregida también se usa en los análisis posteriores que consumen su química.
Las reglas de contactos del receptor y sus umbrales no cambiaron en esta mejora.

Si existen varias correspondencias del grafo, se prioriza compatibilidad con H
explícitos y después con enlaces múltiples/cargas del PDB. En empates se usa la
primera y se emite una advertencia de ambigüedad. Más de 1000 correspondencias
produce un error. Revisar la asignación en grupos simétricos y grupos protonables;
la coincidencia del grafo no garantiza una asignación única de carga/protones
a los nombres de átomos del PDB.

## 6. Recuperar el filtro de anillos anterior

```bash
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A --legacy-rings
```

Para conservar el procedimiento químico previo, ejecutar sin referencia y con
`--legacy-rings`. Esta opción vuelve a exigir más de cinco átomos por anillo.
Puede combinarse con una referencia, pero en ese caso solo restaura el filtro
de tamaño: la percepción química continúa usando la referencia.

Para comparar resultados, usar directorios de trabajo separados con copias de
los mismos PDB y el mismo YAML. Repetir una corrida con los mismos nombres en el
mismo directorio puede sobrescribir archivos por par; los CSV acumulativos,
si están habilitados, agregan filas.

## 7. Revisar las salidas

Cada par genera una carpeta `<receptor>_<ligando>/`. Los nombres y esquemas de
salida se conservan:

- `Interaction_*_all.csv`: contactos candidatos dentro del alcance de búsqueda.
- `Interaction_*_threshold.csv`: contactos filtrados por la distancia de su propio
  tipo de interacción; puede incluir contactos que no pasan el criterio angular.
- `Interaction_*_true.csv`: contactos que pasaron los criterios implementados.
- `summary.csv`: resumen de interacciones por tipo.
- Imágenes del ligando y scripts VMD cuando las opciones correspondientes están
  habilitadas; los anillos de cinco miembros también se incluyen en la imagen.
- Salidas de pockets, ángulos y bias según el flujo y la configuración existentes.

Comprobar que los nombres de átomos y residuos correspondan al sistema y revisar
las diferencias de donores, aceptores, anillos y contactos cargados al incorporar
una referencia. `Interaction = Yes` indica que un contacto pasó las reglas del
programa; el análisis no calcula una energía de unión.

## 8. Errores y acciones

| Mensaje o situación | Acción |
|---|---|
| SMILES inválido o vacío | Revisar la cadena SMILES y pasarla entre comillas. |
| SDF inexistente, inválido o con varias moléculas | Proporcionar un archivo válido con exactamente una molécula. |
| Grafo de átomos pesados incompatible | Confirmar que la referencia sea del mismo ligando y revisar la conectividad leída del PDB. |
| Hidrógenos explícitos incompatibles | Revisar protonación y estructura de ambas entradas antes de repetir. |
| Varias correspondencias | Revisar qué átomos reciben las cargas y protonación en grupos simétricos. |
| Referencia con batch o con un modo sin ligando | Ejecutar por ligando o retirar las opciones que no corresponden al modo. |

Una referencia incompatible produce un error; no se vuelve silenciosamente a la
inferencia desde PDB. Pueden quedar archivos parciales de una ejecución fallida;
su presencia no implica que el análisis haya terminado correctamente.

## 9. Validación registrada

En la implementación de la referencia química se ejecutaron **43 pruebas, todas aprobadas**.
Se verificaron:

- Conservación de coordenadas, nombres y seriales con átomos reordenados.
- Referencias con y sin hidrógenos explícitos; transferencia de enlaces y cargas.
- Identificación de donores y aceptores en amida, amonio y un anillo con `[nH]`.
- Inclusión de anillos de cinco miembros y recuperación del filtro anterior.
- Rechazo de un anillo plano no aromático cuando hay referencia química.
- Rechazo de referencias incompatibles y combinaciones de CLI inválidas.
- Ejecución completa con SMILES y SDF, y pruebas existentes de PDB solo,
  configuración, probes y hotspots.

Las pruebas de integración usan copias temporales de los PDB y el intérprete
activo, sin depender de una ruta fija `.venv/bin/python` ni borrar salidas del
usuario. Las advertencias por correspondencias múltiples en moléculas simétricas
son esperadas en algunos casos de prueba.

Para repetir la suite desde la raíz del repositorio:

```bash
PYTHONPATH="$PWD/src" python -m pytest tests -q
```

También pasaron `git diff --check` y la revisión de Ruff sobre el módulo y las
pruebas nuevos. Estas comprobaciones validan los casos cubiertos; no constituyen
una validación química exhaustiva para todos los ligandos.

## 10. Archivos de implementación

- [ligand_chemistry.py](../src/interactions_search/ligand_chemistry.py): lectura y
  aplicación de la referencia.
- [ligand_hotpoints.py](../src/interactions_search/ligand_hotpoints.py): percepción
  química, selección y visualización de anillos.
- [pipeline.py](../src/interactions_search/pipeline.py) y
  [cli.py](../src/interactions_search/cli.py): integración y parámetros.
- [test_ligand_chemistry.py](../tests/test_ligand_chemistry.py) y
  [test_smoke.py](../tests/test_smoke.py): validación química y compatibilidad.

La refactorización general adicional del núcleo y el análisis de sensibilidad
de hotspots siguen pendientes. La etapa 5 separa la validación geométrica y la
clasificación de contactos del flujo orquestador (sección 14). Los motivos de
decisión y el registro automático están implementados en la sección 13. La unificación de
criterios entre modos está implementada y se describe a continuación.

## 11. Criterios compartidos entre ligando y probe

### 11.1. Configuración

Los detectores comparten búsqueda de vecinos, selección de átomos del receptor
para contactos hidrofóbicos/cargados y reglas de validación en
[interaction_rules.py](../src/interactions_search/interaction_rules.py).
La CLI carga una sola vez la configuración de interacción para construir el
diccionario que consumen los modos de ligando y probe.

Los archivos YAML anteriores siguen funcionando: las claves nuevas son opcionales
y reciben los valores predeterminados indicados abajo. Para personalizarlas,
editar estas entradas en una copia del YAML existente, conservando las tablas de
donores, aceptores y antecedentes:

```yaml
distancias:
  Distances_Hidrogen_Bonds: 3.2
  Hydrogen_Bond_Search_Distance: 4.0
  Distances_Aromatic: 5.5
  Distances_Hidrofobica: 4.0
  Distances_Salt_Bridge: 4.0
  Distances_Pi_Cation: 5.0
  Probe_Clash_Distance: 2.5

angulos:
  Angle_Hidrogen_Bonds_Min: 100
  Angle_Hidrogen_Bonds_Max: 180
  Aromatic_Parallel_Max: 30
  Aromatic_TShaped_Min: 60
```

| Parámetro nuevo | Uso |
|---|---|
| `Hydrogen_Bond_Search_Distance` | Radio de candidatos H-bond, en Å, para ligando y probe. El radio efectivo es el máximo entre este valor y `Distances_Hidrogen_Bonds`. |
| `Distances_Salt_Bridge` | Umbral de puentes salinos, en Å, para ligando y probe; también recuentos de cationes/aniones en el entorno de hotspots. |
| `Distances_Pi_Cation` | Umbral π-catión, en Å, para ligando y probe. |
| `Probe_Clash_Distance` | Umbral de choque de una sonda con átomos pesados, en Å. No sustituye `hotspot_pocket.grid_clash`. |
| `Aromatic_Parallel_Max` | Ángulo entre planos aromáticos estrictamente menor que este valor, en grados. |
| `Aromatic_TShaped_Min` | Ángulo entre planos aromáticos estrictamente mayor que este valor, hasta 90°. |

Las distancias nuevas deben ser positivas y finitas. El mínimo angular de H-bond
debe ser menor que el máximo; el límite paralelo aromático debe ser menor que el
límite de T-shaped, ambos entre 0 y 90°.

### 11.2. Límites y búsqueda de candidatos

Todos los umbrales de contacto usan `distancia < cutoff`: un contacto exactamente
en el límite queda excluido. Los ángulos H-bond usan
`mínimo < ángulo <= máximo`; los límites aromáticos de 30° y 60° quedan excluidos
con la configuración predeterminada.

El filtro fijo de 4 Å ya no limita los H-bonds. Por ejemplo, si se establece
`Distances_Hidrogen_Bonds: 5.0` manteniendo el radio de candidatos en 4.0, la
búsqueda se amplía a 5.0. Un contacto a 4.5 Å llega a la validación angular. Este
ejemplo prueba la aplicación del parámetro, no recomienda 5 Å como criterio
químico de uso general.

La selección espacial del sitio sigue controlada por `centroid_distance`; aumentar
un cutoff de contacto no incorpora automáticamente residuos fuera de ese sitio.
Para comparar modos, usar una selección espacial que contenga los mismos residuos
relevantes. El modo probe centra su selección en las coordenadas de las sondas.

La decisión se toma antes del redondeo de salida. Un contacto justo por debajo de
un límite puede aparecer redondeado al valor del límite en el CSV y seguir siendo
válido. Para contactos hidrofóbicos agrupados se conserva el promedio de las
distancias aceptadas, como antes.

### 11.3. Diferencias geométricas explícitas

| Caso | Geometría y decisión |
|---|---|
| Ligando donor o sonda donor | Mismo ángulo donor–aceptor–antecedente del receptor, usando `acceptors_antecedent` o `C` como nombre por defecto. Si falta el antecedente, se informa `Angle = NaN` e `Interaction = No`. |
| Ligando aceptor | Se conserva el ángulo donor del receptor–aceptor del ligando–antecedente del ligando. |
| Sonda aceptora frente a un donor pesado | Solo distancia: la sonda no define un antecedente propio. |
| Sonda aceptora frente a una fila de donor H explícito | Distancia y ángulo D–H–sonda; se busca el átomo pesado más cercano del residuo como padre. Si no se puede calcular la geometría, se rechaza. |
| Dos anillos reales | Distancia entre centroides y ángulo entre planos, con límites del YAML. |
| Sonda aromática | Solo distancia; un punto no define un plano. El ángulo respecto a la cara del anillo receptor, cuando existe, es informativo. |
| Hidrofóbicos, puentes salinos y π-catión | Validación por distancia sobre los grupos seleccionados; sin filtro angular nuevo. |

Una geometría requerida ausente o degenerada no se transforma en una validación
por distancia. La aceptación solo por distancia de una sonda sin orientación es
una decisión explícita del modo, no una recuperación silenciosa de errores.
La etapa 2 conservó el esquema CSV. La etapa 4 agrega `Reason` (sección 13).

### 11.4. Procedimiento de ejecución y revisión

Usar la misma copia del YAML para el ligando y para sondas ubicadas en puntos
comparables de la pose:

```bash
python Interactions_search.py -r protein.pdb -l ligand.pdb -c A --config config.yml
python Interactions_search.py -r protein.pdb -c A --probe 0 0 0 --probe-type donor --config config.yml
```

Reemplazar `(0, 0, 0)` por la coordenada del donor que se desea comparar. Las
referencias opcionales de la sección 5 siguen disponibles en el comando con
ligando. No se necesita cambiar el comando habitual si se usan los valores
predeterminados.

Al revisar los resultados:

1. `_all` conserva los candidatos dentro del radio de búsqueda, incluidos
   H-bonds que no pasan distancia o ángulo final.
2. `_threshold` del análisis con ligando aplica el cutoff correspondiente a cada
   tipo. Ya no usa `Distances_Aromatic` como filtro global. Puede diferir de
   archivos anteriores incluso con el mismo YAML, por esta corrección.
3. `_true` contiene solo `Interaction = Yes`.
4. Las sondas mantienen sus archivos `_all` y `_true`. Un choque se informa como
   una fila separada `Clash`; no elimina automáticamente las otras filas del punto.
5. Los recuentos de entorno de hotspots usan los cutoffs compartidos y límites
   estrictos. Continúan siendo recuentos por distancia, sin validación angular ni
   cambios en los criterios geométricos de construcción de la cavidad.

### 11.5. Validación de esta etapa

Se ejecutó la suite completa: **98 pruebas aprobadas** (las 43 existentes y 55
nuevas). Las nuevas pruebas están en
[test_interaction_rules.py](../tests/test_interaction_rules.py) y cubren:

- Ejecuciones completas de ligando y probe con el mismo H-bond, incluyendo
  distancias mayores que el antiguo prefiltro, límites exactos y antecedente ausente.
- Equivalencia de detección hidrofóbica, salina y π-catión en los casos comparables,
  incluyendo el redondeo de distancias próximas al límite.
- Aplicación de los cutoffs salino y π-catión del YAML en el pipeline y conservación
  de esos contactos en `_threshold` cuando superan el umbral aromático.
- Aplicación de ángulos aromáticos personalizados en un par real de anillos.
- Cutoff de choque configurable, recuentos de hotspots, rangos inválidos y
  compatibilidad con YAML sin las nuevas claves.

Comando: `PYTHONPATH="$PWD/src" python -m pytest tests -q`.
Pasaron también Ruff para el módulo y las pruebas nuevos, la comprobación de
nombres no definidos en `src/interactions_search` y `git diff --check`.
Las advertencias de correspondencias químicas múltiples en casos simétricos
siguen siendo esperadas.

## 12. Validación con estructuras reales y grupos químicos

La etapa 3 amplía la suite de 98 a **151 pruebas aprobadas**. Consultar el
[informe de validación química](VALIDACION_QUIMICA.md) y la
[procedencia de las fixtures](../tests/fixtures/real/README.md).

Se usan dos recortes cristalográficos originales, con coordenadas sin optimizar:
ABL–imatinib (1IEP, cadena A) y estreptavidina–biotina (1STP, cadena A). Las pruebas
son offline y comprueban 16 contactos seleccionados, con geometría independiente,
además de controles negativos y transformaciones rígidas. No deben interpretarse
como una estimación global de precisión/recall o de energía de unión.

La validación detectó y corrigió donores incorrectos en aminas terciarias neutras,
inconsistencias de ácidos con H explícitos, N aromáticos protonados omitidos,
elección de correspondencias frente a evidencia existente y errores al leer
anillos incompletos del receptor. Los detalles y las tolerancias están en el
informe; los contactos ambiguos del carboxilo de BTN no se fijan como referencia
nominal única.

El YAML distribuido ahora incluye los átomos donores pesados SER OG, THR OG1 y
TYR OH para cristales sin H. Para reproducir esos casos con una configuración
personalizada, revisar estas entradas sin eliminar otras tablas:

```yaml
donors:
  SER: [N, HG, OG]
  THR: [N, HG1, OG1]
  TYR: [N, HH, OH]
```

Esto puede agregar contactos respecto a resultados anteriores. En receptores
protonados pueden coexistir filas del H y del átomo pesado. Los anillos del
receptor con átomos requeridos ausentes se omiten, conservando el análisis de los
demás contactos. No se reconstruyen átomos faltantes.

Ejecutar la suite antes de comparar resultados o cambiar criterios:

```bash
PYTHONPATH="$PWD/src" python -m pytest tests -q
```

## 14. Separación de validación y clasificación

En `pipeline.py`, el análisis de ligando delega ahora la validación angular a
`_validate_interaction_angles()` y la asignación de `Interaction` y `Reason` a
`_classify_interactions()`. La función principal conserva la secuencia de
preparación, búsqueda, exportación y resumen. Las reglas, el orden de filas y
las columnas permanecen iguales; los helpers son privados del módulo.

## 15. Exportar pockets de varios sitios a PDB globales

Después de generar los pockets con `--hotspots`, combinar las grillas de los
sitios construidos:

```bash
python scripts/export_hotspot_pockets.py 3mss_complex_hotspot_pockets \
  --protein-pdb Inputs/3mss_complex.pdb --chain B
```

El script escribe dos archivos en esa carpeta:

- `pockets_global_spheres.pdb`: todos los centros de voxel de cada pocket.
- `pockets_global_surface.pdb`: puntos remuestreados cada ~0.5 Å sobre la
  isosuperficie interpolada de cada grilla; ofrece una superficie más suave y
  compacta que usar los voxels del borde.
- `pockets_global_surface_with_protein.pdb`: por cada sitio, sus puntos de
  superficie y solo los átomos de los residuos proteicos que lo forman. Cada sitio
  usa una cadena propia (A = sitio 1, B = sitio 2, etc.); los resids de proteína
  conservan la numeración de origen.
- Escenas VMD `view_pockets_spheres.tcl`, `view_pockets_surface.tcl` y,
  cuando se pasa `--protein-pdb`, `view_pockets_with_protein.tcl`; estas aplican
  representaciones visuales distintas porque el PDB por sí solo no las define.

En los tres PDB, cada sitio ocupa una cadena: A = sitio 1, B = sitio 2, etc.
Seleccionar `chain A` muestra el sitio 1; `chain B`, el sitio 2. En los PDB solo
de pockets, todos los puntos de esa cadena están en `resid 1`; en el combinado,
usar `chain A and resname PCK` para aislar la grilla y dejar el fragmento proteico
aparte. El nombre de átomo conserva el tipo más favorable: empieza con A
(aceptor), D (donor), H (hidrofóbico) o N (sin tipo asignado). `occupancy`
contiene enterramiento y el factor B contiene Best_DG. La cadena basta para
aislar el sitio en VMD, Maestro y ChimeraX.

El PDB combinado contiene solo los residuos listados en `residues.csv` para los
sitios construidos, copiados a la cadena de cada sitio; residuos compartidos entre
sitios aparecen en ambas cadenas. Para el ejemplo 3MMS son 116 posiciones únicas
(1.916 átomos de proteína en total antes de duplicar las compartidas). Si no se
pasa `--protein-pdb`, se generan solo los dos PDB de pockets.

El archivo `surface` guarda puntos de la isosuperficie interpolada, no las caras
de una malla triangular. Se genera con `scikit-image`; la aplicación calcula la
superficie al visualizar esos puntos. El archivo `spheres` conserva el volumen
completo como puntos y permite mostrar los voxels como esferas. Las coordenadas
de ambos siguen la grilla original (la superficie se interpola entre voxels).
Instalar el extra opcional si aún no está disponible: `python -m pip install -e '.[surfaces]'`.
Para ver la diferencia en VMD, ejecutar desde la carpeta de resultados
`vmd -e view_pockets_spheres.tcl` o `vmd -e view_pockets_surface.tcl`. En Maestro
y ChimeraX hay que asignar manualmente el estilo de esfera o superficie al modelo;
el tipo de representación no forma parte del formato PDB.

Registro de esta etapa: **151 passed**, 40 advertencias esperadas de
correspondencias múltiples, sin fallos esperados (`xfail`). Las referencias de
las etapas anteriores (43 y 98 pruebas) se conservan como registros históricos.

## 13. Motivos y registro automático de corridas

### 13.1. Interpretar los contactos

Los CSV de interacciones de ligando y probe agregan `Reason` al final de las
columnas anteriores (esquema 2). Los nombres de archivo se mantienen. Los lectores
que exijan una lista exacta de columnas deben incorporar esta columna.

| Código | Interpretación |
|---|---|
| `distance_outside_cutoff` | La distancia no cumple el límite estricto |
| `angle_outside_range` | El ángulo no cumple el intervalo configurado |
| `missing_required_geometry` | Falta información geométrica necesaria |
| `distance_and_angle_pass` | Cumple distancia y ángulo |
| `distance_pass` | Cumple el criterio de un contacto evaluado por distancia |
| `distance_only_no_probe_orientation` | Cumple distancia; el probe no tiene orientación |
| `steric_clash` | El probe presenta un choque estérico |

Si falla más de un criterio se registran los códigos separados por `;`. La
evaluación usa valores sin redondear, aunque el CSV muestre distancias o ángulos
redondeados. Solo se explican los candidatos encontrados por los detectores:
la ausencia de una pareja en `_all.csv` no equivale a una fila rechazada.

### 13.2. Archivos generados automáticamente

Los modos ligando, probe, site-bias y hotspots escriben en su carpeta de resultados:

- `config_used.yml`: configuración efectiva, reutilizable con `--config`.
- `run_metadata.json`: identificador, modo, parámetros efectivos, argumentos CLI
  cuando estén disponibles, fechas UTC, duración, versiones de dependencias,
  huella del código y revisión Git cuando esté disponible.
- `run_history/<run_id>/`: copia de la configuración y los metadatos de cada intento.
  Si se utilizó una referencia química, incluye `ligand_reference.sdf`, conservando
  el orden de átomos de la plantilla para mantener el desempate entre mapeos.

Los metadatos contienen rutas y SHA-256 de entradas antes y después del análisis
(incluida la limpieza preexistente del ligando), y de archivos nuevos o modificados
dentro de la carpeta de resultados. Los CSV acumulativos externos quedan fuera
de ese inventario. Se registra también la distancia efectiva de búsqueda de H-bonds.

El estado pasa de `running` a `completed`, `failed` o `interrupted`. Un fallo
incluye el tipo y mensaje del error; no atribuye a ese intento archivos anteriores
que no cambiaron. Un ligando que RDKit no puede leer produce un error explícito.
Los errores de validación CLI anteriores al inicio del análisis no generan registro;
una terminación forzada del proceso puede dejar el estado `running`.

### 13.3. Repetir una corrida

1. Conservar los archivos de entrada originales y el entorno de software.
2. Elegir el registro deseado dentro de `run_history/` y comprobar sus hashes.
3. Usar su `config_used.yml` con `--config` y recuperar de los metadatos el modo,
   cadena, coordenadas, radios y demás argumentos que correspondan.
4. Cuando exista una referencia archivada, usar ese SDF con `--ligand-sdf`;
   conservar también el valor registrado de `legacy_rings`.
5. Ejecutar desde otra carpeta para conservar los resultados anteriores.

El historial conserva configuración, metadatos y referencia química, **no** todas
las entradas ni todos los resultados. Los registros de la raíz describen el último
intento y los archivos de salida habituales siguen pudiendo sobrescribirse.
La configuración por sí sola no sustituye los argumentos del modo ni el entorno.

### 13.4. Verificación

Resultado de la etapa 4: **170 pruebas aprobadas**, con 40 advertencias esperadas
de correspondencias químicas múltiples.

La suite incorpora pruebas de motivos, geometría sin redondear, configuración
recargable, hashes, historial, fallos con archivos previos, limpieza de entradas,
referencia SDF y registros de los cuatro modos. Ejecutar:

```bash
PYTHONPATH="$PWD/src" python -m pytest tests -q
```
