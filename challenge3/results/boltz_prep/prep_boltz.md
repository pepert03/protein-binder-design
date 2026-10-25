# Preparación Boltz-2 complejo GLP1R (SIN EJECUCIÓN, sin GPU, sin sesión nueva)

Estado: PREPARADO. Nada instalado, descargado ni ejecutado. Viabilidad NO prometida.

## 1. Entrada exacta — CORREGIDA con mapa unívoco (25oct 16:08)
- complejo_v6.pdb (sha1 4fd26c70) cadena A: 316 residuos OBSERVADOS en rango de numeración 1..376
  (no "crop 1..316": son 316 residuos con huecos). MAPA FORMAL (verificado por match exacto de secuencia,
  posición única por fragmento): archivo = auth_6X18 − 28; auth_6X18 = UniProt P43220 (offset auth↔UniProt +0,
  coherente con prepare_target.py); por tanto archivo→UniProt = +28 (archivo 1 = UP 29; archivo 376 = UP 404).
  Corrijo mi frase anterior "offset +0" que confundía los dos mapeos.
- Mapa fragmento -> UniProt (match exacto ÚNICO cada uno):
  archivo 1-23 -> UP 29-51 | 33-101 -> 61-129 | 108-136 -> 136-164 | 148-223 -> 176-251 |
  238-300 -> 266-328 | 321-376 -> 349-404.
- Hallazgo clave: los huecos NO son densidad ausente en 6X18. 6X18 cadena R observa 54 de los 60
  residuos del hueco (52-60, 165-175, 252-265, 329-338+344-348) y toda la cola 405-423 que el crop cortó;
  solo 130-135 (6 aa, ECD) falta de verdad en 6X18. Las discontinuidades son ARTEFACTO DE NUESTRO CROP
  (prepare_target.py), no de la estructura. Topología UniProt de lo cortado: 165-175, 252-265, 329-348 =
  loops CITOPLASMICOS (ICL1-3); 52-60 = ECD; 405-463 = cola citoplasmática.
- Fuente primaria local: /tmp/pose/P43220.txt (sha1 739ecaa8ea56c62a566da3737a25f6de6a532b9d; Reviewed,
  seq version 2, 463 aa). FT SIGNAL 1..23 / CHAIN 24..463; TRANSMEM y TOPO_DOM citados en §6c.

## 2. Fallo anterior exacto (registro interno, 6 oct)
- 6oct 11:34:55 CEST: la sesión de estudiante tenía límite de RAM de host de 5,86 GiB; el piloto llegó
  a descarga de pesos, predicción no completada; se añadió guarda de parada por RAM.
- 6oct 12:19:35 CEST: README de evidencia: "No successful structure or confidence output was produced.
  Pilot: n02_rank02_native vs full human EGFR ECD (621 aa), msa:empty for both chains, Boltz 2.2.1,
  default 3 recycling and 200 sampling steps, 1 sample".
- Tras el fallo los plegados se movieron fuera de esa sesión.
- Interpretación acotada: el fallo fue de RAM de host en AQUELLA sesión (5,86 GiB) con un complejo de ~690 res.
  No es una medida del perfil GPU MIG actual ni del cap de su sesión Terminal (desconocido, ver §4).

## 3. Requisitos actuales del modelo (URLs exactos, consultados 25oct)
- https://raw.githubusercontent.com/jwohlwend/boltz/main/README.md — "All the code and weights are provided
  under MIT license"; instalación "pip install boltz[cuda] -U"; para CPU-only quitar [cuda]; "the CPU version
  is significantly slower than the GPU version".
- https://huggingface.co/api/models/boltz-community/boltz-2?blobs=true — tamaños: boltz2_conf.ckpt
  2.286.561.469 B (2,29 GB, estructura+confianza), boltz2_aff.ckpt 2.062.139.170 B (afinidad, no necesario),
  mols.tar 1.855.662.080 B (CCD). Descarga mínima ≈ 4,2 GB + torch/dependencias.
- https://docs.nvidia.com/nim/bionemo/boltz2/latest/support-matrix.html — NIM pide ≥48 GB VRAM, 64 GB RAM,
  12 cores, 80 GB NVMe. ALCANCE: es el contenedor NIM de producción, NO el paquete pip; no usable para
  dimensionar nuestra sesión.
- https://www.rbvi.ucsf.edu/chimerax/data/boltz-apr2025/boltz_help.html — dimensionado empírico (Boltz vía
  ChimeraX, abr2025, linaje Boltz-1): GPU Linux 24 GB → ~1300 res; 12 GB → ~1000; 8 GB → ~700; CPU 16 GB RAM
  → ~350 res; 32 GB → ~1000; tiempo ~cuadrático. ALCANCE: NO es boltz 2.x pip en nuestro entorno; solo
  orienta el orden de magnitud. 352 res queda en el borde del sobre CPU-16GB.

## 4. Recursos reales observados (sin inferir)
- Perfil Estudiante MIG (medido 24-25oct en sesiones previas): MIG 1g.18gb = 16384 MiB VRAM, 4 cores,
  ≤8 h, una sesión GPU a la vez. nvidia-smi dentro de MIG reporta la tarjeta completa; el slice es 16 GiB.
- RAM de host: nodo muestra 2 TB, pero NO es el límite del contenedor. Cap cgroup NO medido (desconocido).
- Dato nuevo honesto (faltaba en mi resultado ESMFold): en la sesión ESMFold del 25oct, intento 1 (fp32 con .half() posterior)
  murió por RAM con pico ~18 GB observado en vivo durante la conversión de precisión (2× pesos fp32 ~5,6 GB +
  copias); intento 2 (fp16) produjo NaN en compute_tm (transformers 5.17.0, loss.py:103); intento 3 (fp32
  + torch_dtype=float32 + low_cpu_mem_usage=True, parche en esmfold_sctm.py sha c5c934ab) SÍ completó los 5
  plegados. Los logs de los intentos 1-2 no se descargaron (sesión borrada); sobrevive run3.log + env.txt
  (env.txt registra: "torch_dtype=float32, low_cpu_mem_usage=True, fp16 descartado (NaN compute_tm)").
  El pico ~18 GB NO acota el cap (pudo morir por debajo o por encima); solo prueba que un pico de esa
  magnitud no es tolerado en esa sesión.
- Sandbox local (medido hoy con check_resources.py): 1,94 GiB efectivos, 2 cores → NO apto para Boltz-2 352 res.

## 5. Bloqueo RAM (ABIERTO, técnico)
El cap efectivo de RAM de la sesión MIG Terminal es desconocido y es el factor que ya tumbó el piloto 6oct.
Resolverlo exige ejecutar check_resources.py DENTRO de una sesión (abrir sesión o usar GPU requiere
aprobación explícita previa). Hasta medirlo, no se puede declarar viable ni CPU ni GPU.

## 6. Artefactos preparados (sha1)
- check_resources.py c121ba10ca30dff0ee64bc472b1a35928c3f66be — mide cgroup/MemAvailable/cores/disco/GPU en la sesión destino; no instala ni descarga; no decide.
- glp1r_top02.yaml 853466623d5a4f9e538153e93b05990b4e3cbdde — complejo 316+36, msa:empty.
- glp1r_nativo.yaml c4caa2feededf93ce16ce5d299e50021f0f6b9c0 — control techo.
- glp1r_aleatorio36.yaml 3d7355f0674fd5b1c148c320eb2689464e78894b — control suelo.
- seqs.json 3caedeb3303a763fbdfc070da807c31b833b5552 — secuencias exactas y su procedencia. CORREGIDO 25oct por el verificador:
  la entrada recA (concatenado descartado) se sustituye por glp1r_24-463 (secuencia exacta de las
  DRAFT, 440 aa); binders sin cambios. El historial git conserva la versión anterior.

## 6b. YAML concatenados = DESCARTADOS (25oct)
Concatenar los fragmentos del crop crea 5 uniones peptídicas artificiales; aun con caveat NO es ruta
defendible para interfaz nativa. Los 3 YAML concatenados quedan renombrados DESCARTADO_*.yaml, marcados
NO EJECUTAR, solo como registro. No se propondrá ejecutar ese input.

## 6c. Constructo RECOMENDADO con base (fuentes primarias locales, 25oct 16:09)
RECOMENDACIÓN: GLP1R 24-463 completo (cadena madura), DRAFT_glp1r_full_*.yaml ya preparados. Razones:
1. Es EXACTAMENTE el constructo experimental de la estructura de referencia que define la interfaz
   nativa: DBREF 6X18 R 24-463 → UNP P43220 24-463 (cabecera de 6X18.pdb, sha1 8985195e); sin mutaciones
   SEQADV en la cadena R (solo tag de expresión N-terminal, REMARK 465, residuos −8..0).
2. Señal de secreción 1-23 excluida por maduración, no por intuición: FT SIGNAL 1..23 / FT CHAIN 24..463
   en P43220 (fuente primaria).
3. Cola citoplasmática 405-463 INCLUIDA: forma parte del constructo experimental real (6X18 la observa
   hasta 423) y quitarla sería una desviación artificial por el extremo C, el mismo tipo de edición que
   invalidó el crop concatenado. Objetivo interfaz: GLP-1 nativo contacta ECD + bolsillo TM; toda esa
   superficie está dentro de 24-463.
4. Cero uniones artificiales: una sola cadena continua, match exacto verificado contra el SQ de P43220.
Alternativa 24-404 (sin cola citoplasmática) NO elegida: truncar no es inválido per se; se
recomienda 24-463 por fidelidad al constructo experimental, no porque la truncada equivalga al
concatenado descartado (problemas distintos).
Coste declarado (fuente UCSF, linaje Boltz-1, SOLO orden de magnitud, no boltz 2.x pip): la referencia
orientativa sitúa ~350 res en el sobre CPU-16GB y ~700 en GPU 8-16GB; 476 res excede la primera en ESA
tabla, pero NO se declara CPU inviable: boltz 2.x pip no está dimensionado y el cap de RAM de host es
DESCONOCIDO. Nada ejecutado; no listo hasta medir cap en sesión y decidir.

## 7. Qué contestaría y qué no (acotado, con constructo continuo)
ipTM/confianza de complejo = hipótesis de vía sobre si un modelo independiente coloca cada secuencia en la
interfaz diseñada. No discrimina necesariamente sin controles de pose/receptor adicionales; confianza de unión
no es agonismo; y 352 res en CPU con RAM desconocida puede simplemente no correr.
