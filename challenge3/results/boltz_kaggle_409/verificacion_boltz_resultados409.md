# Verificación independiente - tanda Boltz con constructo OFICIAL 409 aa (UP 24-432)

Ejecutada por el investigador en Kaggle el 25oct con el target OFICIAL del challenge 3
(UniProt P43220 24-432, 409 aa; ver §8 de `boltz_kaggle/verificacion_boltz_resultados.md`
para la discrepancia de constructo de la tanda anterior). Auditada de forma independiente
antes de publicar. Esta tanda usa el constructo oficial; la tanda de 440 aa
(`boltz_kaggle/`) queda como referencia histórica no-oficial.

NO es validación experimental. NO es una submission al challenge.

## 1. Cadena de custodia - VERIFICADA

- tar recibido: sha1 `4d4302a88c6a41de75280925113a6fe09c50cef6` (coincide con lo declarado).
- 33/33 hashes sha1 de salida recomputados y OK (`out_hashes409.json`).
- Identidad de inputs: las secuencias de los 3 YAML ejecutados son IDÉNTICAS a los DRAFT
  canónicos del repo (`boltz_prep/DRAFT_glp1r_oficial409_*.yaml`, commit 03955d1); solo
  difieren comentarios de cabecera. Receptor 409 aa en los 3 casos; binders 36/30/36 aa.
- Setup idéntico a la tanda 440: Boltz 2.2.1, `--recycling_steps 3 --sampling_steps 200
  --diffusion_samples 1 --devices 1`, SIN MSA (modo single-sequence explícito en el log),
  SIN templates aportados. Log sin errores, 3/3 "Number of failed examples: 0"
  (~2:20 min por predicción, GPU T4 Kaggle).

## 2. Métricas de confianza (recomputadas de confidence_*.json; entre paréntesis tanda 440)

| caso | ipTM | confidence | pTM | complex pLDDT | pLDDT rec / binder (CA) |
|---|---|---|---|---|---|
| top02 | 0.8862 (0.9067) | 0.7774 (0.7378) | 0.9103 (0.8580) | 0.7502 (0.6956) | 0.762 / 0.614 |
| nativo | 0.9337 (0.3662) | 0.8440 (0.6513) | 0.9010 (0.7635) | 0.8216 (0.7226) | 0.831 / 0.695 |
| aleatorio36 | 0.3247 (0.1663) | 0.6974 (0.6268) | 0.8438 (0.7466) | 0.7906 (0.7420) | 0.809 / 0.581 |

El control nativo pasa de ipTM 0.366 a 0.934 con el constructo oficial - el cambio más
grande de la tanda. El ranking top02 > aleatorio36 se mantiene, pero nativo > top02:
el modelo NO discrimina la candidata del péptido nativo.

## 3. Análisis de pose (código propio, `analiza409.py`; referencias EM 6X18 y complejo_v6)

Definiciones idénticas a la tanda 440: core TM = 7 segmentos TRANSMEM FT UniProt P43220
(140-164, 176-201, 228-251, 266-290, 306-328, 349-370, 384-404 = 166 CA); 6X18 cadena R
observada 384 CA (UP 29-423, aguas excluidas); Kabsch sobre CA comunes.

- Receptor vs 6X18: global 8.89 / 6.70 / 9.05 Å (top02/nativo/aleatorio36); core TM
  2.52 / 2.64 / 2.63 Å. El core TM pliega bien; el desvío global viene de ECD y loops.
- NATIVO vs su pose EM real (correspondencia 1:1 justificada, misma secuencia): 14.9 Å
  (frame global 384 CA) / 14.3 Å (frame core TM). MEJORA sustancial frente a 34.3 Å de
  la tanda 440, pero SIGUE SIN reproducir la pose experimental. Sus contactos (47 res)
  cubren ECD (29-145) + boca del bolsillo TM/ECL2 (199-218) + 378-385: la ARQUITECTURA de
  regiones es la correcta (cola del péptido en ECD, extremo hacia el bolsillo), la
  geometría concreta no.
- TOP02 vs pose de diseño (complejo_v6): 21.3 Å CA-RMSD - en la tanda 440 era 5.2 Å.
  Contactos (32 res) en región ortostérica (ECD 27-148, bolsillo 197-234, 298-310,
  371-388). Distancia nearest-CA a la pose EM del nativo: media 8.8 Å, máx 20.8 Å
  (descriptivo, sin mapa 1:1 entre binders de distinta longitud).
- ALEATORIO36: 22 contactos en la cara CITOPLASMÁTICA (254-259, 334-346, 408-430; más un
  residuo aislado en 171) - sitio no ortostérico, igual que en la tanda 440.

## 4. Estabilidad entre constructos (409 vs 440, receptor alineado por core TM)

- Receptor predicho 409 vs 440: 5.03 / 4.67 / 3.22 Å (top02/nativo/aleatorio36).
- Binder/péptido predicho 409 vs 440: top02 15.4 Å, nativo 30.4 Å, aleatorio36 3.6 Å.
- La pose del binder top02 y la del péptido nativo CAMBIAN por completo entre constructos;
  el sitio del control negativo es el único consistente. Con 1 muestra de difusión, sin
  MSA y sin templates, la geometría de unión NO es reproducible.

## 5. Lectura honesta

- Con el constructo oficial, el control nativo obtiene ipTM alto (0.934) y contactos en
  las regiones correctas, pero sigue a 14 Å de su pose EM. ipTM alto NO implica pose
  experimental: es la confianza del modelo en su propia pose.
- top02 mantiene ipTM alto (0.886) y región ortostérica, pero su pose no coincide con la
  del diseño (21 Å) y no es estable entre constructos (15 Å): NO hay evidencia de una
  pose de unión específica y reproducible.
- Esto NO es un diseño validado. Es una comprobación de coherencia con métricas de
  confianza favorables, acotada por: 1 muestra, sin MSA, sin templates, sin contexto de
  membrana/Gs, y poses inestables entre constructos.

## 6. Incidencias de esta verificación (transparencia)

- Saneado de datos privados sobre todo el paquete: sin wamids, IDs internos, rutas de
  usuario ni secretos.
- Script de análisis incluido (`analiza409.py`); resultados crudos en
  `analiza409_result.json`. Log completo en `boltz_run409.log`.
- El set core TM (166 CA) se reconstruyó desde FT TRANSMEM de UniProt P43220
  (https://rest.uniprot.org/uniprotkb/P43220.txt) y coincide con el usado en la tanda 440.
