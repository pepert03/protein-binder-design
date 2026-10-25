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

---

# v2 (25 oct 14:4x): pasada CON receptor como contexto fijo - verificada

El investigador corrigio el metodo (su primera pasada era binder aislado, como marco
el verificador): complejo_bb.pdb = receptor COMPLETO + backbone binder; secuencia del
receptor en minusculas (recseq_lower.txt, 316 aa) como contexto fijo; choques contados
contra el receptor ORIGINAL (coords 6X18 exactas), no contra la salida de FASPR.

## Verificacion del verificador (reproduccion independiente, tarball extraido con revision de rutas)
- Mecanismo: la cadena A de la salida conserva los 2575 atomos del receptor, 0 ausentes,
  desplazamiento max 0.067 A (drift <=0.07 confirmado). El contexto SI entro en el packing.
- Backbone binder intacto: drift CA = 0.000 A.
- Reconteo propio RESTR_3 contra receptor original: 46 choques, dmin 1.46 (NH1 B22 -
  CH2 A186) = EXACTAMENTE lo reportado.
- **Control nativo REPRODUCIDO directamente desde 6X18 real (cadena P vs R, SIN
  packing): 3 choques, dmin 2.46 = lo reportado.** Ya no es provenance pendiente.
- Tabla de las 32 publicada (pack_ctx_results.txt): rango 46-90 choques, dmin
  0.56-1.46 A; mejor = RESTR_3 (46, 1.46).

## Lectura honesta
- Las 32 candidatas sobre esta pose rigida quedan a 46-90 solapes del receptor frente
  a 3 del nativo: interfaz muy lejos de nativa bajo este criterio.
- dmin de 0.56-1.46 A (atomos casi coincidentes en varias) indica que el scoring
  interno de FASPR no impidió solapes severos en la interfaz: cribado, no optimizador
  fisico. Sigue siendo evidencia del FALLO del protocolo rigido ensayado (rotamer
  library + backbone fijo + sin relajacion), NO prueba de imposibilidad universal.
- El investigador RETIRO la busqueda global a 3.4 A (el nativo no la cumple): correcto.
- Pendiente menor: version/commit de FASPR y sha256 de bb_binder.pdb/receptor_A.pdb.
