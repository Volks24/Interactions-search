# Validación química y cristalográfica — etapa 3

Fecha: 2026-09-25. Resultado: **151 pruebas aprobadas**, sin fallos ni pruebas
marcadas como fallo esperado. Esta etapa incorpora 53 pruebas a las 98 previas.

## Datos y criterio de referencia

Se usaron dos estructuras de difracción de rayos X de familias diferentes:

- [1IEP, ABL–imatinib](https://www.rcsb.org/structure/1IEP), cadena A, STI 201;
  resolución reportada de 2.10 Å. Referencia primaria: Nagar et al., *Cancer
  Research* 62 (2002), 4236–4243, vinculada en la entrada RCSB.
- [1STP, estreptavidina–biotina](https://www.rcsb.org/structure/1STP), cadena A,
  BTN 300; resolución reportada de 2.60 Å. Referencia primaria: Weber et al.,
  *Science* 243 (1989), 85–88,
  [DOI 10.1126/science.2911722](https://doi.org/10.1126/science.2911722).

Los PDB originales y las definiciones químicas del CCD se descargaron de RCSB.
Se extrajeron sitios pequeños y se guardaron sus fuentes, selección y SHA-256 en
[manifest.json](../tests/fixtures/real/manifest.json). La suite usa esas copias
locales; no necesita red, `Inputs/` ni un servidor de análisis externo.

Los contactos de control fueron seleccionados por identidad de átomo y grupo
químico. Un cálculo independiente lee las columnas del PDB y calcula distancias
y ángulos sin importar las funciones geométricas ni de detección del programa.
Las expectativas no son una copia de sus CSV. **Yes/No indica si ese par satisface
las reglas geométricas configuradas; no implica que la cristalografía haya
medido un puente de hidrógeno como una observación binaria.**

## Contactos controlados

Se comprueban 16 contactos concretos: 11 candidatos H-bond, 3 aromáticos y 2
hidrofóbicos. Incluyen 13 positivos y 3 negativos según los criterios por defecto.

| Estructura | Átomo del ligando | Residuo y átomo del receptor | Distancia (Å) | Ángulo (°) | Resultado |
|---|---|---|---:|---:|---|
| 1IEP | N3 | MET318 N | 2.899 | 119.971 | Yes |
| 1IEP | O29 | ASP381 N | 2.901 | 112.100 | Yes |
| 1IEP | N13 | THR315 OG1 | 2.883 | 126.950 | Yes |
| 1IEP | N21 | GLU286 OE2 | 2.996 | 156.650 | Yes |
| 1IEP | N21 | ASP381 O | 3.668 | 132.553 | No: distancia |
| 1STP | O3 | ASN23 ND2 | 2.992 | 114.930 | Yes |
| 1STP | O3 | SER27 OG | 2.771 | 127.997 | Yes |
| 1STP | O3 | TYR43 OH | 2.585 | 126.426 | Yes |
| 1STP | N1 | ASP128 OD2 | 2.784 | 123.622 | Yes |
| 1STP | N2 | SER45 OG | 3.044 | 134.280 | Yes |
| 1STP | N1 | ASP128 OD1 | 3.683 | 77.514 | No: distancia y ángulo |

Las identidades y cifras de H-bond están fijadas en
[contacts.json](../tests/fixtures/real/contacts.json). Se comprueba la aparición
en `_all`, la clasificación y la inclusión/exclusión correspondiente en
`_threshold` y `_true`.

Los controles aromáticos de 1IEP comparan el primer anillo de STI con PHE317 y
el segundo con TYR253 y PHE382. Los dos primeros son positivos; PHE382 se rechaza
por la orientación intermedia. Un ajuste independiente del plano con todos los
átomos verifica la clasificación. Los controles hidrofóbicos son STI C6–TYR253
CE1 y BTN C9–LEU110 CD2.

## Pruebas adicionales

Las 27 pruebas de [test_real_structures.py](../tests/test_real_structures.py)
incluyen los contactos anteriores, integridad de fixtures, identidad de grupos
químicos/anillos, ejecución con PDB solo y estos controles:

- Rotación rígida de 90° y traslación del complejo: se conservan las identidades
  y clasificaciones en los tres CSV.
- Traslación del ligando en 100 Å por eje respecto del receptor: no quedan contactos.
- Eliminación de CE1 de PHE317: ese anillo deja de generar un centro aromático,
  mientras continúan los contactos válidos del resto del sitio.
- STI neutro: las aminas terciarias N48/N51 son aceptores sin H donor; el N de
  amida no es aceptor. Se verifican por nombre los cuatro anillos aromáticos.
- BTN: los N ureido donan pero no aceptan; sus dos anillos de cinco miembros no
  se clasifican como aromáticos con referencia química.

Se añadieron también 26 pruebas a
[test_ligand_chemistry.py](../tests/test_ligand_chemistry.py): 12 grupos/estados de
protonación, con y sin H explícitos, y dos controles de correspondencia de
oxígenos equivalentes. Cubren aminas terciarias neutras y protonadas, amonio
cuaternario, ácido/carboxilato, sulfonamida, nitrilo, piridina/piridinio, pirrol,
amida terciaria y sulfona. Las expectativas de cantidad de donores/aceptores son
explícitas; se exige que todo donor tenga H en el estado químico asignado.

## Errores detectados y corregidos

1. **Donores potenciales frente a donores del estado suministrado.**
   `BaseFeatures.fdef` contempla aminas terciarias como posibles donores. En el
   modo con referencia se reemplazó por los SMARTS específicos HBA/HBD de
   `rdkit.Chem.Lipinski`, añadiendo el caso de N aromático protonado con H.
   Una amina terciaria neutra ya no se etiqueta como donor.
2. **Diferencias por representación explícita de H.** Los SMARTS se evalúan
   sobre una copia con H colapsados. Esto evita que el enlace O–H de un ácido
   satisfaga una rama de aceptor de alcohol al hacer el matching con H explícitos.
   Los índices se traducen de vuelta a la molécula original, que conserva sus H,
   coordenadas y nombres.
3. **Correspondencias químicas equivalentes.** Se priorizan las correspondencias
   compatibles con H explícitos y luego con enlaces múltiples/cargas existentes
   del PDB. En empates se conserva la primera y se informa la ambigüedad.
   Más de 1000 correspondencias produce un error en vez de elegir un subconjunto
   arbitrario. No se resuelve de forma única un grupo cuando la entrada carece
   de evidencia suficiente.
4. **Donores del receptor ausentes en cristales sin H.** El YAML distribuido
   incorpora OG de SER, OG1 de THR y OH de TYR junto a sus nombres de H. Esto
   recupera, entre otros, los contactos de BTN O3 con SER27 y TYR43.
   Los YAML personalizados conservan sus tablas propias y deben actualizarse
   si se desea esa misma cobertura. En estructuras protonadas pueden coexistir
   filas del átomo pesado y del H, como ya ocurría con otros residuos.
5. **Anillos incompletos del receptor.** Si falta un átomo requerido para el
   centro aromático, se omite ese centro; ya no se produce el `IndexError`
   asociado a ese caso. No se reconstruyen coordenadas faltantes.

El procedimiento químico del ligando en modo PDB solo se conserva. Los cambios
de tablas del receptor y manejo de anillos incompletos sí afectan ambos modos.

## Tolerancias y límites

- H-bonds: 0.00001 Å en distancia y 0.001° en ángulo frente al cálculo
  independiente sobre las mismas coordenadas.
- Hidrofóbicos: 0.0005 Å, por redondeo de salida a tres decimales.
- Aromáticos: 0.01 Å y 2° frente al centro/plano de todos los átomos. El código
  histórico usa cuatro átomos para el centro del receptor, redondeado a 0.001 Å,
  y tres para el plano. Las etiquetas se comparan exactamente.
- Transformación rígida: tolerancia numérica absoluta 0.001, por el redondeo del
  centro del receptor; identidades y clasificaciones deben permanecer iguales.

Son dos sitios cristalográficos recortados, no un benchmark exhaustivo. No se
estima precisión/recall global, energía de unión, dependencia de pH ni contactos
con aguas, metales o cadenas generadas por simetría. La referencia CCD es neutra:
no se exige que STI N51 forme un H-bond como si estuviera protonado. En BTN,
O11/O12 pueden intercambiarse al asignar el ácido neutro si faltan órdenes de
enlace e H; no se fijan contactos nominales de ese carboxilo como verdad única.

Estas pruebas tampoco certifican todas las heurísticas de puentes salinos o
π-catión; la cobertura por distancias de esos detectores está descrita en la
etapa 2 del [SOP](SOP.md).

## Reproducir

Desde la raíz, con las dependencias del proyecto instaladas:

```bash
PYTHONPATH="$PWD/src" python -m pytest tests -q
PYTHONPATH="$PWD/src" python -m pytest tests/test_real_structures.py tests/test_ligand_chemistry.py -q
```

La ejecución completa registrada terminó con **151 passed** y 40 advertencias
esperadas de correspondencias múltiples. No se ocultaron fallos como `xfail`.
Pasaron Ruff para el script de extracción, el módulo de referencia y las pruebas
de esta etapa, y `git diff --check`.

La [guía de las fixtures](../tests/fixtures/real/README.md) explica su reconstrucción
desde las fuentes originales; ese paso es separado de la ejecución de pruebas.
