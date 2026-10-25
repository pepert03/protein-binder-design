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

## scaffold_v4.pdb - SUPERA la verificacion geometrica (con reservas; NO es un diseno)

Verificacion independiente (escritor, 25 oct; ademas crosscheck del investigador con
biopython 1.88 como segunda implementacion, reproducido en el sandbox del escritor):

- Procedencia: regenerado atomo a atomo (2719, 0 mismatches) desde
  scripts/exporta_scaffold4.py (seed 23, determinista, SOLO coordenadas del receptor;
  pose por muestreo + helice NeRF ajustada por Kabsch a la traza).
- Longitudes: 0 fuera de tolerancia (bb_check y crosscheck: N-CA 1.46, CA-C 1.52,
  C-O 1.23, C-N 1.33).
- Angulos/torsiones: dentro de tolerancia en residuos internos (N-CA-C 111, CA-C-N 116,
  C-N-CA 121, omega +-180, phi -58, psi -47; crosscheck biopython: idem). DOS defectos
  terminales PENDIENTES DE CORRECCION antes de cualquier diseno de secuencia: B1
  N-CA-C = 90 grados y B36 CA-C-O = 180 grados (O colineal). No son cosmeticos: v4 NO
  es input valido para ProteinMPNN hasta que los extremos esten corregidos.
- Choques: 0 internos no enlazados (check_scaffold) y 0 backbone-receptor (gate G4).
- Gate de orientacion: PASS (G1=3, G2=10, G3 OK, G4=0).
- Calibracion (test_nativo.py, SOLO calibracion, no base de diseno): la helice canonica
  ajustada a la traza nativa de GLP-1 da RMSD 1.44 A con 1 choque (min 1.89 A);
  RMSD hélice ideal res1-12 vs 6X18 U233-244 = 0.344 A. Sirven para calibrar el
  constructor, no para inicializar disenos.

RESERVAS (ninguna resuelta):
- QUIRALIDAD: NO establecida. El PDB no lleva CB; cualquier afirmacion "L" deriva de un
  pseudo-CB reconstruido por convencion, no del archivo. Afirmacion L retirada
  explicitamente. Se fijara al construir secuencia real.
- Sidechains: CONFLICTO ABIERTO en B6 y B28: CB APROXIMADO a <2 A del receptor (1.83 y
  1.45 A; posiciones CB estimadas por bb_check, no reales). NO esta demostrado que la
  eleccion de residuo pequeno o de rotamero lo resuelva: queda como reserva, no como
  solucion. Puede requerir revisar la pose.
- Sin secuencia, sin sidechains reales, sin plegado, sin funcion. Nada predice union ni
  activacion. Es una HIPOTESIS DE POSE con backbone geometricamente valido (salvo
  extremos), candidata a diseno de secuencia, no un diseno.

Siguiente paso propuesto (docs/siguiente_paso.md, NO ejecutado): ProteinMPNN en CPU
sobre la cadena B. Pausado hasta revision; la validacion AF2 posterior requeriria GPU
con autorizacion explicita.
