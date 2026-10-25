# Packing FASPR de las 32 candidatas sobre scaffold_v5 - revision del metodo (25 oct)

Artefacto del investigador: pack_eval.py + mejor pose RESTR_13 (55 solapes; rango
reportado 55-103 en las 32; nativo reportado 3). CPU, sin GPU. Cribado, no validacion.

## Verificado por el verificador (reproduccion independiente del CONTEO)
- RESTR_13.pdb: backbone identico a v5 (RMSD CA 0.0000 A, 36/36; FASPR fija backbone).
- Recuento propio binder x receptor con el MISMO criterio (solape vdW >0.4 A,
  radios N 1.55 / C 1.70 / O 1.52 / S 1.80): **55 solapes = lo reportado**. dmin
  1.72 A (CZ B16 - CH2 A63).
- Choques INTERNOS del binder (excluidos enlazados y vecinos secuenciales): **0**.
  El packing interno de FASPR es limpio; TODO el conflicto es binder-receptor.
- Distribucion por posicion (mi recuento): B19=13, B16=7, B28=7, B18=5, B30=4,
  B27=3, B31=3, B6=2, B32=2, B1-B5/B8/B15/B35<=1. Se concentra en la cara al receptor.

## Limites del metodo (por que NO concluye inviabilidad)
- **Packing SIN contexto del receptor**: FASPR eligio rotameros viendo solo el binder
  aislado; el receptor entra unicamente en la evaluacion posterior. Un packer con
  contexto (o relajacion) tiene grados de libertad que este protocolo no exploró:
  los 55-103 solapes son una cota del protocolo rigido, no prueba de rotameros
  imposibles ni de que ninguna relajacion lo recupere.
- Backbone fijo: ninguna relajacion de backbone evaluada.
- Control nativo = 3: REPORTADO, NO REPRODUCIDO (falta el PDB nativo empaquetado).
- El script no evalua choques internos (los cubre la minimizacion de FASPR; verificado
  independientemente = 0 para esta pose).

## Provenance PENDIENTE del investigador
Tabla completa de las 32 secuencias; version/commit de FASPR; sha256 de bb_binder.pdb
y receptor_A.pdb; PDB del control nativo empaquetado. RESTR_13 se acepta como "mejor"
en tanto su 55 reproduce y es el minimo del rango reportado, pero la tabla falta.

## Otros
- Propuesta de busqueda sidechain-aware a 3.4 A: heuristica NO validada; sin promesa.
- No se publica render: sin verificacion visual disponible; la tabla por posicion de
  arriba es el dato verificado.
