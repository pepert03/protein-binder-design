# Preparación verificada: auto-consistencia ESMFold (sin ejecutar modelo, sin GPU aún)

QUÉ CONTESTA: cómo pliega ESMFold cada secuencia SOLA y a qué distancia queda la salida del modelo
(CA-RMSD Kabsch + pLDDT) del backbone diseñado. Es salida de UN modelo, no evidencia de que la
secuencia "prefiera" el pliegue.
QUÉ NO CONTESTA: unión a GLP1R, calidad de interfaz, agonismo, ni pliegue en contexto del receptor.
Controles EXPLORATORIOS (no techo/suelo garantizados): GLP-1 nativo 30 aa aislado y permutación
de top02 (misma composición, seed 7).
Umbrales de referencia con fuente: scTM>=0.5 (FoldingDiff, Nat Commun 2024, PMC10844308);
criterio AF2 RMSD backbone <2 A + pAE<5 (RFdiffusion, Nature 2023, s41586-023-06415-8).

SELECCIÓN DE CANDIDATAS (criterio explícito): samples 2, 14 y 5 de las 16 secuencias MPNN sobre v6,
elegidas SOLO por menor número de choques vdW medidos (21, 23, 24). No es un ranking de mes ni de
bondad global; global_score MPNN no es confianza calibrada.

VERIFICADO EN CPU (sin modelo): API transformers 5.20.2 (EsmForProteinFolding.infer_pdb por
inspección), Kabsch CA-RMSD (identidad = 0.0), extracción de refs/FASTAs, asserts de longitud,
parser por ID de cadena sin depender de TER (regresión B36).
Reconciliación provenance: binder scaffold_v6[B] == complejo[B] == packed LIBRE6_2[B] con
RMSD 0.0000 exacto; el 1.96 A anterior fue artefacto de mi parser (ventanas sobre cadena A).
Exactitud entre entornos LIMITADA (verificador, commit 0c92c9a): misma preparada, 1 hilo,
E0 igual, pero mínimo difiere lig RMSD 0.893 vs 0.863 y contactos 553 vs 550 entre máquinas:
valores de cola NO portables entre entornos; dirección sí.
Cobertura: scaffold_v6 receptor tiene 11 sidechains recortadas a CB (resseq renumerados
46,48,51,69,75,85,86,87,94,101,108) que el complejo lleva completas (reconstruidas por FASPR):
métricas de choques entre esos archivos NO son parejas; usar siempre el mismo archivo por comparación.

PESOS: facebook/esmfold_v1, licencia MIT (tarjeta HF). Modelo = ESM-2 3B + 48 bloques
(fuente: facebookresearch/esm esmfold/v1/pretrained.py). Sin MSA, sin servicios externos:
todo local tras descargar pesos.

RECURSOS DGX — ESTIMADOS SIN VERIFICAR (hardware/VRAM/RAM del DGX desconocidos ahora; no afirmo
disponibilidad actual por estado antiguo): pesos 3B fp32 ~12 GB / fp16 ~6 GB; >=16 GB GPU esperado
para 30-36 aa. En frío: descarga de varios GB + carga + segundos-minutos por péptido corto.
Primer paso EN sesión: nvidia-smi + free -g antes de prometer tiempos.

PERMISO: este documento NO autoriza nada. La tanda autorizada hasta ahora cubre generar y
evaluar candidatos, con código y resultados en el Umbrel, sin envío al concurso. ESMFold es una
NUEVA ejecución fuera de esa tanda: requiere decisión explícita del usuario antes de lanzarse.

EJECUCIÓN (cuando se autorice): 
python3 esmfold_sctm.py --fasta v6_top02.fa v6_top14.fa v6_top05.fa aleatorio36.fa --ref ref_binder_v6.pdb --ref-n-ca 36 --out out/ --fp16
python3 esmfold_sctm.py --fasta nativo_glp1.fa --ref ref_nativo.pdb --ref-n-ca 30 --out out/ --fp16
