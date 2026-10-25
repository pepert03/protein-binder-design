# Auditoria preparacion Boltz-2 (NO ejecutada) — verificador, 25 oct (entrega final)

## Mapa de numeracion: VERIFICADO y corregido el encabezado
- Archivo complejo_v6.pdb cadena A resseq r -> auth 6X18 / UniProt P43220 r+28: 316/316 residuos
  EXACTOS (archivo 1 = UP 29, archivo 376 = UP 404). El texto "offset +0" de la prep es incorrecto
  y contradice su propia tabla; el mapa correcto es +28. 6X18 auth = UniProt en esta region.
- 6X18_RP cadena R: 404 res, auth 29..520 (huecos propios; >463 = tags de fusion).

## Huecos del crop: son artefacto del crop, no de la estructura (con cifra corregida)
- Huecos del crop en numeracion auth: 52-60, 130-135, 165-175, 252-265, 329-348 (60 res).
- Presentes en 6X18 R: 49/60 (52-60, 165-175, 252-265, 329-338, 344-348). Ausentes de verdad:
  130-135 (ECD) y 339-343 (ICL3) = 11 res. La prep decia "54 de 60": su propio desglose suma 49;
  verificado 49. Cola 405-423 observada en 6X18 (19/19) y cortada por el crop.
- Consecuencia: los 5 cortes del crop eliminan residuos que la estructura SI observa.

## Concatenados: DESCARTADOS (verificado el motivo)
- 5 uniones peptidicas artificiales con CA-CA reales 16.4/13.7/15.7/12.4/8.0 A (mediana contiguos
  3.82 A) = declarado. Un predictor con ese input modela un constructo ficticio. Descarte correcto.

## DRAFT continuo: secuencia VERIFICADA contra fuente primaria
- Cadena A de los 3 DRAFT == UniProt P43220 (Reviewed, SV=2, 463 aa) posiciones 24-463, 440 aa,
  EXACTO (verificado contra rest.uniprot.org, no contra el archivo local citado).
- Binders == top02/nativo/aleatorio36 auditados. msa:empty en ambas cadenas.
- sha1 medidos: DRAFT top02 fb287f9d, nativo 355ec1b1, aleatorio36 d2eb53f2; DESCARTADO top02
  d3962639, nativo 94ed319d, aleatorio36 9869f1fc; check_resources.py c121ba10; seqs.json 445a5ebc.
- DEFECTO: seqs.json NO se actualizo — sigue llevando la secuencia CONCATENADA descartada como
  "recA". Corregir o regenerar antes de cualquier uso como indice.

## Estado y limites
- NADA ejecutado; sin GPU, sin sesion nueva. Viabilidad RAM/VRAM/pip boltz2 DESCONOCIDA (cap cgroup
  no medido; 2 TB del nodo NO es el limite del contenedor). check_resources.py seguro (solo mide).
- DRAFT NO ES EJECUTABLE-LISTO: constructo pendiente de recomendacion cientifica justificada y
  aprobacion del usuario. Incluye regiones SIN densidad 6X18 (24-28, 130-135, ICL1-3, 339-343,
  cola 405-463), sin contexto de membrana ni Gs: la salida seria pose del modelo, no el EM 6X18.
- Incidente ESMFold registrado en results/esmfold/esmfold_verificacion.md: reportado, no
  reproducible, no invalida los hashes auditados.

## Constructo recomendado 24-463: base factual VERIFICADA (RCSB)
- DBREF de la cabecera oficial 6X18 (files.rcsb.org/header/6X18.pdb): "DBREF 6X18 R 24 463 UNP
  P43220 GLP1R_HUMAN 24 463" — el constructo experimental de la estructura de referencia ES
  UniProt 24-463. SEQADV cadena R: 0 conflictos. VERIFICADO.
- Detalle menor NO verificado aqui: el tag de expresion N-terminal (residuos -8..0, REMARK 465
  segun la prep) no se comprobo residuo a residuo; no afecta al constructo recomendado.
- Mapa formal: archivo = auth-28; auth = UniProt (+0); archivo->UniProt +28. La prep corrige su
  frase "offset +0" (confundia los dos mapeos); conforme con mi verificacion 316/316.
- GLP-1 nativo: DBREF 6X18 P 7-36 -> UNP P01275 98-127, coherente con la cadena P usada.
- Correcciones aplicadas por el verificador al publicar (prosa acordada): recomendacion por
  FIDELIDAD de constructo (no por invalidez de truncar); CPU NO declarado inviable por la tabla
  orientativa Boltz-1 (cap RAM/boltz 2.x sin verificar); saneo de IDs internos de la prosa.
