# Verificación independiente - campaña zero-shot 409 aa v2 (FINAL, 26oct 09:16 CET)

Versión FINAL de la campaña zero-shot sobre el target oficial 409 aa (aprobación del usuario
25oct 23:21, verificada en su canal). SUSTITUYE a v1 (`zeroshot_409_v1/`, conservada como
historial): el delta es exactamente el filtrado de gA_3, más la regeneración de
results.json / results_seed2.json / contacts409.json / logs y MANIFEST. Todo lo demás es
byte-idéntico a v1 (ya auditado; ver `zeroshot_409_v1/verificacion_zeroshot_409_v1.md`,
cuyas etiquetas y matices se mantienen íntegros aquí).

TODO ES PREDICCIÓN COMPUTACIONAL. NO validación experimental, NO evidencia de agonismo,
NO submission.

## 1. Cadena de custodia v2 - VERIFICADA

- tar sha256 `5d706427cc621bc8433073953eb16f7cd02e1d3ffe10b011ed05461f93198c3a` (coincide).
- MANIFEST.sha256: 40/40 hashes recomputados OK.
- Delta vs v1 verificado por diff: solo gA_3 (YAML + 2 cifs) y los JSON/log/MANIFEST regenerados.
- gA_3: receptor del YAML == target oficial 409 (sha1 689647a8...); binder == fasta == PDB == results.
- Nota menor: la ruta del cif de gA_3 seed2 omite el nivel `boltz_results_f409_gA_3/`
  (inconsistencia de empaquetado, sin efecto en el contenido).
- Saneado: limpio.

## 2. gA_3 (la novedad) - auditado

- Métricas declaradas (mismas limitaciones que v1: sin confidence_*.json crudos; corroborado
  por logs y pLDDT de los cifs, diff <=0.002): ipTM 0.9310 (seed1) / 0.9234 (seed2),
  confidence 0.792/0.788, pTM 0.913/0.918, complex pLDDT 0.758/0.754.
  ipTM seed1 de gA_3 es el MAYOR de los tres candidatos (.931 > gB_2 .915 > gB_1 .834).
- Contactos: REPRODUCIDOS EXACTOS (13/13 conjuntos del contacts409.json con mi parser).
  gA_3: 56/61 contactos totales y overlap epítopo 30 (seed1) / 35 (seed2), en región
  ortostérica (ECD 24-148 + bolsillo TM 197-233, 296-301, 376-388).
- Estabilidad de pose entre semillas (mi Kabsch, frame core TM 166 CA): binder RMSD 11.6 A;
  Jaccard de contactos 0.48; overlap de epítopo compartido entre semillas 23 de 30/35.
  MISMO MATIZ que gB_1/gB_2: métricas estables, pose concreta no.
- Secuencia (ya verificada en v1, idéntica): MEDGSGKKTTYSPYVDYMINELDMWPVSGFAA, 32 aa,
  SIN tramos de baja complejidad (9% G, max run 1), sin tags, 0 k-mers >=4 exactos vs
  GLP-1/exendin-4/liraglutide/semaglutide/taspoglutide/top02. Es la secuencia de mejor
  calidad de los tres candidatos (gB_1 y gB_2 llevan colas G4).

## 3. Cuadro final de candidatos (predicción, 1 muestra/semilla, sin MSA, sin templates, sin membrana/Gs)

| candidato | ipTM s1/s2 | epítopo s1/s2 | RMSD pose entre semillas | complejidad |
|---|---|---|---|---|
| gA_3 | 0.931/0.923 | 30/35 | 11.6 A | limpia |
| gB_2 | 0.915/0.924 | 48/37 | 7.5 A | cola G4 (19% G) |
| gB_1 | 0.834/0.907 | 40/38 | 12.4 A | cola G4 (19% G) |
| nativo (control) | 0.934/0.944 | 22/44 | 13.9 A | - |
| aleatorio36 (control) | 0.292/0.245 | 0/0 | 3.5 A | - |

## 4. Lectura honesta final

- gA_3 pasa los 3 criterios del investigador en dos semillas, igual que gB_2 y gB_1. Ningún
  candidato supera al control nativo en ipTM, y el nativo recupera en seed2 más epítopo que
  cualquiera (44). El filtro NO discrimina los diseños del péptido nativo.
- Ningún candidato tiene pose geométricamente reproducible entre semillas (7-12 A).
- Esto NO es un diseño validado ni evidencia de unión o agonismo: métricas internas del
  modelo, corroboradas solo en su cadena de custodia y consistencia.
