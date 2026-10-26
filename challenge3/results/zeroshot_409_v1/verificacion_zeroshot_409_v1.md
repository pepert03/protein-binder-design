# Verificación independiente - campaña zero-shot 409 aa v1 (26oct 08:45 CET)

Campaña de diseño de novo sobre el target OFICIAL del challenge 3 (UniProt P43220 24-432,
409 aa), ejecutada por el investigador en Kaggle (aprobación del usuario 25oct 23:21),
auditada de forma independiente antes de publicar.

TODO ES PREDICCIÓN COMPUTACIONAL. NO validación experimental, NO evidencia de agonismo,
NO submission al challenge.

## 1. Cadena de custodia - VERIFICADA

- tar sha256 `c45817c2f7b814d47a40b534c8dec1fb8665168945444d61b06abf439d0e3179` (coincide).
- MANIFEST.sha256: 37/37 hashes recomputados OK.
- Nota: el paquete contiene 39 archivos (37 en MANIFEST + README + MANIFEST); se declararon 41.
  No falta nada según el MANIFEST.
- Receptor de los 6 YAML de filtro: 409 aa, sha1 de secuencia == target oficial
  (`689647a8c8b799bdf43c9b753e1e83ff5446b1f7`) en los 6 casos.
- Binders: YAML == `designs/sequences.fasta` == `designs/*.pdb` (cadena B) == results.json, en los 5 diseños.
- Controles nativo y aleatorio36: idénticos a los canónicos de la tanda 409 (commit 03955d1).
- Saneado de datos privados: limpio (sin wamids, IDs internos, rutas ni secretos).

## 2. Generación de los diseños (logs y script auditados)

- ColabDesign AF2 `protocol="binder"`, `binder_len=32`, SIN secuencia inicial: generación
  de novo real. Los hotspots usados son contactos del TARGET (lado receptor), no un binder
  de partida: la regla zero-shot/de novo prohíbe partir de un binder existente, no usar
  información del target.
- Diseño sobre el crop v3 (316 res, 6 fragmentos); filtro Boltz sobre el 409 completo.
- Trayectorias: gA 3/3 completas; gB 2 completas (la 3ª quedó cortada por presupuesto GPU 36000 s).
- gA_3 fue DISEÑADO pero NUNCA pasó por el filtro Boltz (no hay YAML ni predicción suya);
  queda publicado como diseño sin filtrar (sus métricas AF2: i_ptm 0.495, plddt 0.443).

## 3. Métricas Boltz 2.2.1 (setup idéntico a tandas anteriores: recycling 3, sampling 200, 1 muestra, sin MSA, sin templates)

Declaradas por el investigador; los confidence_*.json CRUDOS no se incluyeron, así que
ipTM/pTM/confidence no son recomputables de forma independiente. Corroboración parcial:
los valores de los logs reproducen results*.json al decimal, y el complex_plddt recomputado
de los .cif (media all-atom del B-factor) difiere <=0.004 de lo declarado en los 6 casos.

| caso | ipTM seed1 | ipTM seed2 | complex pLDDT s1/s2 |
|---|---|---|---|
| gB_2 | 0.9145 | 0.9237 | 0.768 / 0.800 |
| gB_1 | 0.8336 | 0.9068 | 0.767 / 0.775 |
| gA_2 | 0.8162 | 0.8784 | 0.753 / 0.736 |
| gA_1 | 0.3560 | (sin seed2) | 0.769 / - |
| nativo | 0.9336 | 0.9441 | 0.844 / 0.832 |
| aleatorio36 | 0.2918 | 0.2446 | 0.786 / 0.793 |

## 4. Contactos de interfaz - REPRODUCIDOS EXACTOS

Recomputados con parser propio e independiente (misma definición que el investigador:
átomos pesados del receptor a <5.0 A del binder, receptor rn+23 = UniProt). Los 11
conjuntos de contactos y los 11 solapamientos con el epítopo coinciden EXACTAMENTE con
`filter_outputs/contacts409.json`.

Epítopo (`docs/epitope_contacts.tsv`, derivado de 6X18): 57 filas = 56 residuos del
receptor + 1 fila de agua (513, nunca puede contar). Incluye residuos definidos por
contacto con aguas cristalográficas (HOH101-105), no solo con el péptido GLP-1.

Overlap con epítopo (seed1/seed2): gB_2 48/37, gB_1 40/38, gA_2 0/0, gA_1 0/-,
nativo 22/44, aleatorio36 0/0.

## 5. Estabilidad entre semillas (mi análisis, frame core TM 166 CA FT UniProt)

| caso | RMSD CA binder seed1 vs seed2 | Jaccard contactos |
|---|---|---|
| gB_2 | 7.5 A | 0.55 |
| gB_1 | 12.4 A | 0.56 |
| gA_2 | 75.6 A (sitio distinto) | 0.00 |
| nativo | 13.9 A | 0.41 |
| aleatorio36 | 3.5 A | 0.83 |

La estabilidad reclamada se sostiene en MÉTRICAS (ipTM y engagement del epítopo se mantienen
altos en seed2 para gB_2/gB_1) pero NO en geometría: las poses concretas se desplazan
7-12 A entre semillas. El control negativo es geométricamente más estable que los candidatos.

## 6. Comprobaciones de secuencia (regla de novo / zero-shot)

- Sin tags ni secuencias de señal (His6, FLAG, etc.) en ningún diseño.
- Identidad vs péptidos conocidos: CERO k-mers exactos de >=4 residuos contra GLP-1(7-36),
  GLP-1(7-37), exendin-4, liraglutide, semaglutide y taspoglutide; mejor ventana sin gaps
  15-38% de identidad. También distintos de la candidata anterior top02 (0 k-mers).
  Consistente con de novo: NO derivan de binders conocidos.
- Baja complejidad (a anotar, no descalifica): colas con tramos de Gly - gA_1 run G6 (19% G),
  gA_2 run G5 (31% G), gB_1 run G4 (19% G), gB_2 run G4 (19% G). gA_3 sin tramos (9% G).
  Las identidades 50% entre gA_1/gB_2 y gA_2/gB_2 son artefacto de alinear esas colas Gly.

## 7. Lectura honesta

- Veredicto del investigador REPRODUCIDO en sus propios términos: gB_2 y gB_1 pasan sus 3
  criterios en dos semillas (ipTM alto, contactos de epítopo altos, métricas estables);
  gA_2 y gA_1 rechazados (gA_2 por pose inestable/sitio errante, gA_1 por ipTM 0.36).
- Matices que el filtro NO resuelve: (a) la estabilidad es de métricas, no de pose
  (7-12 A entre semillas); (b) el control nativo en seed2 recupera MÁS contactos de
  epítopo que cualquier candidato (44 vs 38/37) con ipTM mayor (0.944): el filtro no
  discrimina los diseños del péptido nativo; (c) colas Gly de baja complejidad en gB_1/gB_2.
- Esto NO es un diseño validado ni evidencia de unión o agonismo: 1 muestra de difusión,
  sin MSA, sin templates, sin membrana/Gs, y métricas de confianza internas del modelo.

## 8. Incidencias de esta verificación (transparencia)

- `env/environment.json` tiene los campos de entorno cdenv/boltzenv VACÍOS (solo registra GPU).
- Faltan los confidence_*.json crudos (ver §3).
- Script de auditoría incluido: `audita_zs1.py`.
