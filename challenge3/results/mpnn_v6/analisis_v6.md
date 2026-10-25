# scaffold_v6 + MPNN v6 - analisis verificado (25 oct)

## Provenance del archivo: COMPLETA
- scaffold_v6.pdb regenerado atomo a atomo desde exporta_scaffold6.py (PARAMS fijos,
  builder NeRF+Kabsch verificado en v4/v5): 2719 atomos, 0 mismatches. sha256
  ccd1db10...d09a = declarado.
- TRAYECTORIA DETERMINISTA REPRODUCIDA (no solo el comentario del exportador):
  refine_wide 1500 -> superv=23, mejor=9 (EXACTO al registro declarado); refine2 1200
  anclado en el mejor de wide -> superv=306, mejor=4, C-ECD=6.02, PARAMS finales
  EXACTOS a los exportados. Misma SEED=20261025, mismo FASPR (commit 0d55732, lib
  dun2010bbdep.bin sha256 ed3f7be5...d0e; ambos = declarados, compilados de fuente).
  Nota: la etapa global laxa del script adjunto da 0/4400 (no alimento la cadena
  final); la cadena verificada es refine_v5 -> refine_wide -> refine2. Registros de
  etapa (cands/pack_rank/best_final) no archivados; la reproduccion los sustituye.

## GATE: seleccion post-hoc, documentada SIN relajar criterios
- Corte ORIGINAL (heredado v5): n_bols>=3 y d_ecd<=6.0. La busqueda v6 RELAJO durante
  el refinamiento a n_bols>=2 y d_ecd<=8. La pose final: n_bols=5 (cumple), d_ecd=6.02
  **INCUMPLE el corte original por 0.02 A**. Se reporta tal cual: v6 es una seleccion
  post-hoc bajo un gate relajado; NO declarar cumplimiento del corte original ni
  comparar con v5 como si fuera el mismo filtro.
- Metrica distinta del repo (orient_gate v2, G2 = nº residuos 27-36 a <=6 A del ECD):
  PASS (G1=5, G2=3, G3 limpio, G4=0). gate_pose: 0 pares <2.0 A, min 2.5 A.
  Las dos metricas C-ECD no son la misma; quedan etiquetadas por separado.

## Quimica/geometria (verificacion independiente)
- crosscheck biopython: phi -57.8..-57.7, psi -47.1..-46.9, omega +-180, N-CA-C
  110.9-111.0, enlaces nominales; RMSD 0.343 A vs helice real U233-244.
- bb_check: 0/36 CB<2.0 A; 0 colisiones internas. cb_map (formula validada 301 CB
  reales): 0 conflictos <2.0 A. 0 choques backbone-receptor.
- Quiralidad: el archivo NO lleva CB; la clase L es de la CONSTRUCCION calibrada,
  no observable del PDB. No se afirma quiralidad del scaffold.

## Packing poly-ALA en el bucle (proxy declarado, NO mejor caso universal)
- REPRODUCIDO: v6=4 solapes vdW>0.4, v5=44 (mismo FASPR, mismo criterio). Control
  nativo 6X18 real sin packing: 3, dmin 2.46.
- poly-ALA = sidechain minima fija (Gly no tiene CB; Ala solo CB): NO es el mejor
  caso universal ni validacion de secuencia. Umbral <=5 exploratorio, no validado.

## MPNN v6 (16 secuencias, CPU, sin GPU)
- FASTA verificado: 16 muestras, 36 aa, unicas; cabecera v_48_020/git 8907e66/seed 42/
  T=0.1; sha256 contenido 0b520cc8...da9a = declarado. Secuencias = tabla_v6.txt.
- Scores MPNN = log-verosimilitud condicionada al scaffold: NO confianza de validacion.
- max 8-mer vs GLP-1: 3-4/8 (reproducido): proxy debil, no prueba diversidad
  estructural ni zero-shot ni cumplimiento de regla.
- Packing FASPR de las 16 (contexto receptor fijo): REPRODUCIDA la tabla completa
  (21-49 solapes, dmin 1.04-2.10; mejor sample2 21/2.10). Sigue lejos del nativo (3):
  cribado del protocolo rigido, no prueba de imposibilidad. Nada validado: sin
  plegado independiente, sin funcion. Cualquier afirmacion de plegado requiere
  estructura nueva con autorizacion de GPU.
