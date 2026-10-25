# ESMFold auto-consistencia CA-RMSD — tanda DGX 25oct (EJECUTADA)

SCOPE: aprobado por el usuario 15:47: top3 por choques + 2 controles, comprobando
recursos/coste cero antes. Sin submission, sin servicios externos, sin pagos.

ENTORNO (comprobado en sesión, no estimado): DGX H200, MIG 1g.18gb (16384 MiB), RAM nodo 2 TB,
4 cores, python 3.10.12, torch 2.14.0+cu126 CUDA True, transformers 5.17.0 (infer_pdb verificado).
Modelo facebook/esmfold_v1 (HF, MIT), fp32, low_cpu_mem_usage=True. fp16 descartado: NaN en
compute_tm (loss.py:103) — ver run2.log... (no exportado; env.txt lo registra). Primer intento
OOM-kill (pico fp32->half en CPU); corregido con carga streaming. 5 inferencias < 1 min total.

RESULTADOS (CA-RMSD Kabsch vs backbone diseñado / pLDDT medio):
- v6_top02 (21 choques): 0.695 A / 0.92
- v6_top14 (23):         0.480 A / 0.90
- v6_top05 (24):         0.640 A / 0.92
- aleatorio36 (permutación top02, misma composición, seed 7): 0.607 A / 0.89
- nativo GLP-1 30 aa aislado: 6.426 A / 0.64

LECTURA HONESTA (limitada): con ESTE control, el test no separó candidatas de ESTA permutación
(0.607 vs 0.48-0.70; pLDDT 0.89 vs 0.90-0.92). Un solo control no permite generalizar a
"cualquier secuencia helicoidal". Conclusión acotada: auto-consistencia CA-RMSD sobre esta hélice
única no aportó evidencia discriminadora a favor ni en contra de las candidatas EN ESTA corrida.
El nativo aislado dio RMSD alto (6.43 A, pLDDT 0.64): es la predicción de ESTE protocolo sobre la
secuencia aislada; NO prueba dependencia de contexto.
NOTA METRICA: se calculó CA-RMSD Kabsch + pLDDT, NO TM-score: no es scTM en sentido FoldingDiff.
pLDDT en escala 0..1 (transformers 5.17.0 escribe pLDDT/100 en B-factors; verificado en PDB).
Ninguna lectura dice nada sobre unión a GLP1R ni agonismo.

VERIFICACIÓN LOCAL POST-DESCARGA: Kabsch recalculado en sandbox sobre los PDB descargados
reproduce EXACTAMENTE los 5 pares de valores. Hashes DGX == sandbox (7/7 archivos).

SESIONES: terminal 710 y Jupyter 713 borradas tras exportar (cuota GPU liberada). Jupyter 712
cayó en el bug de plataforma "sesion pertenece a otro usuario"; borrada y relanzada (713 OK).
Lancé por error una Jupyter con GPU (711, mi select-option falló por flag --tab-id); quedó Queued
y la borré de inmediato; nunca ejecutó. Terminal 710 quedó a 8h por el mismo fallo de select;
se usó ~10 min y se borró manualmente.

ARTEFACTOS: esmfold_results/ con 5 PDB + run3.log + env.txt + este md. Hashes:
fcf5426eb5cf5c96b4e80d798f86eb4aa5f989278b1faa640eb8e8b68ff6d8b3  aleatorio36_esmfold.pdb
4c6345b0f11a9722dc4612b63b07087e9621a924cd8caa525d3d542d8250921d  nativo_glp1_esmfold.pdb
2fefee215775a95828e7926e813db3efdba41d156acec180306ff4440f4c7b33  v6_top02_esmfold.pdb
3bf8697783b3b118236a0f6c16b1c4e990e1bec2c73bb3eb61f94562a3f89247  v6_top05_esmfold.pdb
cfd85fc15d47309d75ca561480848186f3075a9c055439e75f09bbe9006fa144  v6_top14_esmfold.pdb
36d266d90a26ce6e5f8f5ff91f3d49e81da4693756d5c859ad87e78338b8b45c  run3.log
13457c84f7c058c680ce2811a5b61bb645ad81e514f119eccda557cda4452e19  env.txt
