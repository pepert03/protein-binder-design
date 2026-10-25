# Auditoria independiente relax3 = PROTOCOLO A original (receptor masa 0) (verificador, 25 oct)

## Config verificada
relax3.py usa system.setParticleMass(i, 0.0) para todos los atomos del receptor (no del
ligando): es el protocolo A pre-declarado. Coherente con la fuente oficial OpenMM
(gradiente de masa-0 anulado en ReferenceMinimize.cpp): las tres salidas tienen
receptor CA RMSD = 0.000000 EXACTO (max 0.000000), verificado atomo a atomo. El
investigador retiro su afirmacion anterior sobre masa-0; la retirada es correcta.
Esto confirma ademas que la variante A' (k=10000) era innecesaria.

## Reproduccion independiente (mi codigo, mapeo secuencia-explicito) - EXACTA
- C1 NATIVO:    choques 3->3 | dmin 2.46->2.45 | contactos 980->944 | lig RMSD 0.441 (max 0.836) | receptor 0 EXACTO
- C2 poly-ALA:  choques 4->0 | dmin 2.36->2.94 | contactos 245->248 | lig RMSD 0.372 (max 1.066) | receptor 0 EXACTO
- C3 LIBRE6_2:  choques 21->3 | dmin 2.10->2.60 | contactos 601->577 | lig RMSD 0.918 (max 3.331) | receptor 0 EXACTO
Anclaje C-term (metrica reconciliada): 16.65 -> 17.88 (sigue ALEJANDOSE con receptor fijo exacto).

## Lectura (cribado, no validacion) - sin rankings
- Con receptor EXACTAMENTE fijo, la candidata absorbe 21->3 solapes con deriva de
  backbone 0.92 A (2.1x el nativo, 0.44) y PERDIDA neta de contactos (601->577).
  El nativo queda 3->3 con deriva menor y tambien pierde contactos (980->944).
- NO declarar contactos conservados (disminuyen en ambos). NO ordenar 3 vs 3: el
  criterio vdW>0.4 post-relajacion no esta calibrado (marca la interfaz nativa) y
  los sistemas difieren; no es un ranking de calidad ni funcion. Sin GPU.

## ANOMALIA E0 (NO RESUELTA, documentada)
E0 candidata en protocolo A = 831595 kJ/mol vs A' = 124608 kJ/mol con LA MISMA entrada
(sha256 73994371...80). Verificacion: las salidas relA y rel2 tienen EXACTAMENTE los
mismos hidrogenos por cadena (2837 total, distribucion identica por cadena): misma
protonacion y topologia. La diferencia de E0 solo puede venir de la POSICION INICIAL
de hidrogenos colocados por PDBFixer/addHydrogens (no bit-determinista entre ejecuciones)
o de otra preparacion no determinista. Consecuencia: las energias E0/E1 entre corridas
NO son comparables; solo las metricas sobre atomos pesados lo son. Pendiente: entradas
preparadas (post-PDBFixer) y logs completos para reconciliar. Sin resolver.
