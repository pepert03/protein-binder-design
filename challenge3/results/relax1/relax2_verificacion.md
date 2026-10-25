# Auditoria independiente relax2 (verificador, 25 oct)

## Desviacion de protocolo: la justificacion del cambio es FALSA (verificado en fuente oficial)
El protocolo pre-declarado fijaba el receptor con masa=0. El investigador lo sustituyo por
restriccion armonica k=10000 kJ/mol/nm2 sobre TODOS los atomos del receptor alegando que
"LocalEnergyMinimizer ignora las masas y masa 0 no inmoviliza".
Verificacion contra el codigo oficial OpenMM (la version usada: 8.6.1):
- LocalEnergyMinimizer.minimize -> ContextImpl::minimize -> MinimizeKernel de la plataforma.
- CpuPlatform HEREDA de ReferencePlatform (platforms/cpu/include/CpuPlatform.h l.51) y usa su
  MinimizeKernel registrado -> ReferenceMinimize::minimize
  (platforms/reference/src/SimTKReference/ReferenceMinimize.cpp).
- En ReferenceMinimize.cpp l.56: el gradiente L-BFGS de toda particula con masa == 0 se pone a
  CERO (g[3i..3i+2]=0.0), y la proyeccion de restricciones tambien las excluye (l.98-103).
  Gradiente cero => la particula NO se mueve en la minimizacion.
Conclusion: masa=0 SI inmoviliza en el minimizador usado. El cambio de protocolo era
innecesario; el protocolo A pre-declarado habria funcionado tal cual.
Lo ejecutado es por tanto la VARIANTE A' (receptor fuertemente restringido), con receptor CA
RMSD 0.06 A (max 0.29): pequeno pero NO cero. Retiradas por el investigador (correctamente)
las etiquetas "inmovil", "interfaz intacta" y "native-like integral".

## Reproduccion independiente (mi codigo, mapeo secuencia-explicito) de los 3 sistemas
Hashes de entradas/salidas = relax2_log.txt. TODOS los numeros reproducidos EXACTAMENTE:
- C1 NATIVO:    choques 3->2 | dmin 2.46->2.55 | contactos 980->992 | lig RMSD 0.369 (max 0.888) | receptor 0.0599
- C2 poly-ALA:  choques 4->0 | dmin 2.36->2.91 | contactos 245->257 | lig RMSD 0.388 (max 1.392) | receptor 0.0615
- C3 LIBRE6_2:  choques 21->2 | dmin 2.10->2.61 | contactos 601->582 | lig RMSD 0.868 (max 3.170) | receptor 0.0623
Anclajes con la metrica reconciliada (CA primer/ultimo residuo -> U367-391 / U120-138):
- LIBRE6_2 pre 7.35 / 16.65 -> rel2 7.13 / 17.83 (N-term se acerca 0.2; C-term SE ALEJA 1.2)
- Nativo pre: 8.19 / 15.14. ALA rel2: 7.87 / 17.31.
El "enganche C-term preservado" de relax1 queda RETIRADO tambien en esta metrica: 20.00 (relax1)
confirmado; en rel2 sigue alejandose. La interfaz se sostiene solo por conteo global de
contactos, que en C3 ademas BAJA (601->582).

## Lectura (cribado, no validacion)
- Bajo receptor casi-fijo, la candidata absorbe 21->2 solapes pero con deriva de backbone
  2.3x la del nativo (0.87 vs 0.37) y PERDIDA neta de contactos; el control nativo mejora
  con geometria mucho mas estable y gana contactos. poly-ALA 4->0 con deriva minima:
  los 4 solapes polyALA eran "de biblioteca de rotameros", no comparables a secuencia real.
- Vacio, amber14, un minimizador local, un criterio vdW estricto que marca hasta la interfaz
  nativa: estas metricas describen el protocolo, NO calidad de diseno, NO funcion.
- No se amplia a las 16 ni se ejecuta nada mas sin revision. Sin GPU.
