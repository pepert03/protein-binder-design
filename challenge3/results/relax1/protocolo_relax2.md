# Protocolo relax2 (PREPARADO, NO EJECUTADO — pendiente de reproduccion de relax1)

Objetivo: aislar la explicacion del cambio choques 21->1 en LIBRE6_2, separando contribucion
de movimiento de backbone binder vs receptor, con receptor REALMENTE fijo.

Protocolo A (receptor fijo duro): masa 0 a todos los atomos del receptor en OpenMM
(particle mass = 0 => excluidos de la integracion/minimizacion), binder libre,
misma FF (amber14), vacio, L-BFGS tol=10 maxIter=500. Sin restriccion CA (innecesaria con receptor fijo).
Variante A': ademas CA del binder con k=100 para acotar deriva de pose.

Controles (mismo protocolo):
- C1: nativo 6X18 R-P (peptido libre, receptor fijo): esperado ~0 deriva si protocolo sano.
- C2: candidata con coordenadas del binder SIN empaquetar (backbone ALA): negativo de secuencia.
- C3: LIBRE6_2 con receptor fijo: lectura principal.

Lecturas pre-declaradas: choques/dmin/contactos (audita_relax.py), CA RMSD binder,
energia (solo relativa, no afinidad), y DERIVA vs CHOQUES: si choques bajan solo con RMSD>2 A,
la explicacion es "separacion", no "acomodo local". Criterio de exito pre-declarado: choques <= nivel
del control nativo en el MISMO protocolo con RMSD binder <=1 A.

No ejecutar hasta luz verde; CPU local, ~1-2 min por sistema. Sin GPU.
