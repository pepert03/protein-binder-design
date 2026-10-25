# MPNN CPU sobre scaffold_v5 - analisis verificado (25 oct)

Artefactos del investigador (ejecutados en su entorno CPU, reportados: 16/lote,
seed 42, T=0.1, modelo v_48_020). **Comandos exactos y logs NO suministrados**:
provenance parcial, recuperada de cabeceras FASTA y contenido (ver abajo).
Sin GPU. Sin plegado nuevo: son SECUENCIAS sobre el scaffold fijo; nada de esto
valida plegado ni funcion.

## Cabeceras (ambos archivos, identicas)
- modelo v_48_020, git_hash 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57, seed=42, T=0.1
- fixed_chains=['A'], designed_chains=['B']
- 16 muestras por archivo, 36 aa, todas unicas, alfabeto AA valido

## Asignacion libre/restringida (recuperada del contenido, NO de nombres de archivo)
El investigador no sabia que archivo era cual. Prueba de cumplimiento contra
omit_B3B6.jsonl (posiciones B3 y B6 restringidas a {A,G,S,T,V}):
- `scaffold_v5_LIBRE.fa` (hash adjunto 453d904d): pos3 viola en 16/16, pos6 en 15/16 -> LIBRE
- `scaffold_v5_RESTRINGIDA.fa` (hash adjunto 5893bd87): pos3 0/16, pos6 0/16 violaciones -> RESTRINGIDA
  (pos6 colapsa a G en 16/16; pos3 solo {A,V})

## Hallazgo critico: las restricciones NO cubren el conflicto esterico real
Comprobacion independiente de posiciones estrechas sobre scaffold_v5 (dos metodos):
- bb_check (CB aproximado, repo): B6=1.83 A, B28=1.45 A (los unicos <2.5 A)
- Reconstruccion CB tetraedrica independiente (N,CA,C, quiralidad L):
  B28=1.86 A, B32=1.97 A; siguiente B19=2.57 A (B6=3.23 A por este metodo)

Conclusion:
- **B28 es el unico conflicto confirmado por AMBOS metodos y quedo SIN restringir**
  (muestra {A,E,R} en ambas corridas; E/R son cadenas largas).
- B6 solo aparece con el metodo aproximado (ambiguo por orientacion de CB).
- B3/B4 (criterio del investigador) NO aparece en ningun metodo sobre v5; el valor
  B3 procedia de scaffold_v3, ya descartado. Restringir B3 no corrige nada conocido.
- Las restricciones aplicadas NO son una solucion demostrada a los conflictos CB:
  son una heuristica sobre posiciones parcialmente equivocadas.

## Novelty (proxy debil, limitaciones explicitas)
- Maxima identidad en ventana de 8 aa vs GLP-1(7-36): 3-5/8 (ambas corridas)
- vs exendina-4(1-30): 3-4/8
- Esto NO demuestra diversidad estructural, ni cumplimiento de la regla de
  diversidad de la competicion, ni zero-shot: es solo un proxy de secuencia local.

## Alcance
Gates/RMSD sobre el scaffold miden fidelidad de construccion de la pose, no
validan estas secuencias ni su plegado. Cualquier afirmacion de plegado/funcion
requiere estructura nueva (re-plegado independiente) con autorizacion de GPU.
