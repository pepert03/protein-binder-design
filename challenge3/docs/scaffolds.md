# Scaffolds geometricos (Challenge 3, GLP1R) - estado verificado

Fecha: 2026-10-25. TODO aqui es trabajo geometrico local SIN GPU. Un scaffold es una
HIPOTESIS DE POSE, no un diseno: sin secuencia, sin sidechains, sin plegado, sin funcion.
Nada de esto predice union ni activacion. Marco de reglas: de novo/zero-shot (ver
viabilidad_rediseno.md §0); las poses se buscaron SOLO con coordenadas del receptor.

## scaffold_v3_FALLA_QUIMICA.pdb - NO USAR: falla la verificacion quimica

Reportado por el investigador como "geometria REAL de backbone (0 enlaces fuera de
tolerancia)". La verificacion independiente (escritor, 25 oct) dice:

- bb_check.py confirma 0 LONGITUDES de enlace fuera de tolerancia (N-CA 1.46, CA-C 1.52,
  C-O 1.23, C-N peptido 1.33, tol 0.1 A). PERO bb_check solo mide longitudes y una
  holgura de CB aproximada ("simplificacion" del propio script): no mide angulos,
  diedros ni quiralidad.
- Revision de angulos/diedros (escritor, sobre el archivo): CA-CA consecutivos = 1.05 A
  constantes (real ~3.80); angulo N-CA-C 69-90 deg (real ~111); CA-N-C ~59 deg
  (real ~121); omega ~47 deg (no trans); phi constante +57.8 deg (signo de helice
  zurda). La cadena B es un zigzag degenerado, NO una helice ni un backbone valido.
- Quiralidad del CA: NO verificable (el archivo no lleva CB); bb_check tampoco la mide.
- Procedencia: la traza CA del archivo NO coincide con los PARAMS_V3 publicados en
  siguiente_paso.md (desviacion media 8.6 A, max 15.8 A); el script exporta_scaffold3.py
  no se ha compartido. Pose origen NO verificable.
- orient_gate.py da PASS (G1 12, G2 9, G3 OK, G4 0 choques): el gate mide orientacion y
  choques, NO quimica. Un PASS del gate NO certifica un backbone valido (limitacion del
  gate, ahora documentada en su cabecera).

Consecuencia: v3 no es input valido para ProteinMPNN ni para nada. El siguiente paso real
es una reconstruccion de backbone con angulos/diedros correctos (o un constructor de
backbone helicoidal real) antes de cualquier diseno de secuencia.

## scaffold_v2_DISCARDED.pdb - descartado (25 oct)

Pose del muestreo target-only (0 choques en el modelo idealizado, N-bolsillo 12/12,
C-ECD 2.5 A; verificada: regenerada atomo a atomo desde exporta_scaffold2.py seed 23, y
receptor identico al crop de 6X18, 316 residuos). DESCARTADO como input: su backbone es
decorativo, no quimico (C_i-N_i+1 = 3.77 A sin enlace, N-CA = 1.32, N-CA-C = 180 deg,
C-O = 0.80 A, sin CB). Se conserva como evidencia y para regresiones del gate.

## Herramientas

- scripts/orient_gate.py: gate G1-G4 robusto (exit 0/1/2). Exige receptor crop completo
  (316 resseqs del mapeo P43220 verificado) y binder 1..36; error explicito en datos
  incompletos, columnas desalineadas o cadenas ausentes. Sin registros MODEL: un modelo
  con AVISO. Regresiones: tests/test_gate.py (8 casos).
- scripts/gate_pose.py: gate esterico por contactos; ahora aborta (exit 2) si no parsea
  MODELs o faltan cadenas.
- scripts/bb_check.py: longitudes de enlace + holgura CB aproximada. SOLO longitudes:
  no certifica angulos, diedros ni quiralidad.
- scripts/exporta_scaffold2.py, scripts/verifica_scaffold_v2.py: generador/verificador
  de v2 (artefacto descartado; conservados para procedencia).
