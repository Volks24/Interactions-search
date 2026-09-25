# Fixtures cristalográficas de validación

Estas copias locales permiten ejecutar las pruebas sin conexión. No proceden de
las salidas del detector ni del PDB preparado del directorio `Inputs/`.

| Entrada RCSB | Selección | Átomos receptor / ligando |
|---|---|---|
| [1IEP](https://www.rcsb.org/structure/1IEP) | Cadena A, STI 201; ABL–imatinib | 504 / 37 |
| [1STP](https://www.rcsb.org/structure/1STP) | Cadena A, BTN 300; estreptavidina–biotina | 378 / 16 |

Se conservaron los residuos completos con al menos un átomo a ≤8 Å del ligando,
los nombres/seriales/coordenadas originales y los registros CONECT del ligando.
Solo se seleccionaron altlocs vacíos o A. No se añadieron H ni se optimizó la
geometría; se excluyeron aguas, otras cadenas y contactos de simetría cristalina.
Estas fixtures son recortes de sitios, no ensamblajes biológicos completos.

`STI.smi` y `BTN.smi` contienen el descriptor `SMILES_CANONICAL` de CACTVS del
Chemical Component Dictionary, descargado el 2026-09-25. Ambos componentes de
referencia son neutros. Ver [STI](https://www.rcsb.org/ligand/STI) y
[BTN](https://www.rcsb.org/ligand/BTN). No se presupone que esa protonación sea la
predominante a un pH experimental concreto.

[manifest.json](manifest.json) registra URLs, SHA-256 de fuentes y fixtures,
identificadores y selección. [contacts.json](contacts.json) contiene pares de
átomos seleccionados manualmente, con distancias/ángulos calculados desde las
coordenadas originales usando operaciones matemáticas independientes del código
de producción. Las etiquetas Yes/No expresan los criterios geométricos
predeterminados del programa; no son una anotación experimental de puentes de H.

## Regeneración controlada

Guardar como `1IEP.pdb`, `1STP.pdb`, `STI.cif` y `BTN.cif` las fuentes de las URLs
del manifiesto en un directorio local. Comprobar sus hashes antes de sustituir
fixtures: una nueva revisión de RCSB exige revisar las expectativas.

Desde la raíz del proyecto:

```bash
python scripts/build_validation_fixtures.py --source-dir /ruta/fuentes --output-dir /tmp/fixtures_reconstruidas
PYTHONPATH="$PWD/src" python -m pytest tests/test_real_structures.py -q
```

El script reconstruye PDB, SMILES y manifiesto, pero deliberadamente **no genera
`contacts.json` ni lee resultados del detector**. Comparar los archivos
reconstruidos antes de actualizar las copias versionadas.

Los contactos del carboxilo de BTN no fijan qué nombre O11/O12 debe recibir el
protón: sin H ni órdenes de enlace explícitos esa correspondencia puede ser
ambigua. Los controles nominales de BTN se centran en su grupo ureido y sus
contactos hidrofóbicos.

Ver [informe de validación](../../../docs/VALIDACION_QUIMICA.md) para cobertura,
tolerancias, correcciones y límites.
