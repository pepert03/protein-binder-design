# MPNN CPU sobre scaffold_v5 - analisis verificado v2 (25 oct, reconciliado)

Artefactos del investigador ejecutados en su entorno CPU (torch 2.14.1+cpu, modelo
v_48_020, seed 42, T=0.1, cadena B disenada, A fija). Sin GPU. Sin plegado nuevo:
son SECUENCIAS sobre el scaffold fijo; nada de esto valida plegado ni funcion.

## Provenance: COMPLETA y verificada (provenance.md + logs en este directorio)
- Entrada: scaffold_v5.pdb sha256 7931f785...87da = archivo del repo (commit 6b33d50). Coincide.
- Repo dauparas/ProteinMPNN commit 8907e66, pesos v_48_020.pt sha256 c9cb4a67...f5bd (declarado).
- Comandos exactos y logs de re-ejecucion publicados; re-runs byte-identicos (determinista, seed 42).
- Hashes de CONTENIDO sha256 de los FASTA publicados aqui:
    LIBRE       62675e45...b546  = salida declarada de la corrida libre. Coincide.
    RESTRINGIDA 6299aefc...3de8  = salida declarada de la corrida restringida. Coincide.
  (Los sufijos 453d904d/5893bd87 del transporte eran hashes del attachment, no del contenido:
  no hubo corrupcion; la asignacion LIBRE/RESTRINGIDA queda confirmada por DOS vias:
  cumplimiento del omit y hashes declarados.)

## Datos basicos (verificados)
16 muestras por archivo, 36 aa, todas unicas, alfabeto valido. Cabeceras: modelo
v_48_020, git 8907e66, seed=42, T=0.1, fixed_chains=['A'], designed_chains=['B'].

## Mapa CB RECONCILIADO (sustituye a todas las listas anteriores)
Tres implementaciones independientes (canonica ProteinMPNN l.968-970; tetraedrica
109.47 grad con quiralidad por datos reales; reimplementacion propia del verificador),
todas validadas contra los 301 CB REALES del receptor en el propio archivo:
med error 0.056 / 0.064 / 0.056 A; quiralidad contraria FALSADA (med 2.44 A).
Mapa unico, distancia minima CB->receptor con par explicito:
- CONFLICTOS (<2.0 A): B2=1.17 (CD1 A356) | B6=1.71 (CD1 A113) | B31=1.96 (CA A90)
- LIMPIOS: B28=2.68 | B32=3.52 | B3=3.63 | B4=4.30 | borde: B11=2.03, B35=2.58
RETRACTADO: el mapa previo del verificador (B28/B32, construccion especular con
error med 2.35 A contra CB reales) y el CB aproximado de bb_check (B6/B28).
B3/B4 del investigador venia de scaffold_v3, descartado.

## Efecto sobre las 32 secuencias (verificado posicion a posicion)
- B2: GLICINA en 32/32 (ambas corridas). Sin CB: el conflicto no se materializa.
- B31: GLICINA en 32/32. Idem.
- B6: corrida restringida -> G en 16/16 (evitado); corrida LIBRE -> L en 15/16
  (choque CB real si se usara tal cual).
- B3 (restringida innecesariamente): limpio a 3.63 A.
Conclusion honesta: los conflictos reales B2/B31 se evitaron por glicina ESPONTANEA
de MPNN, no por diseno ni por la restriccion; la restriccion B6->AGSTV solo fue
efectiva porque el muestreo eligio G (A,S,T,V tienen CB y chocarian). Ninguna de
estas evasiones es garantia para futuras campanas: cualquier candidata con residuo
con CB en B2/B6/B31 choca; las opciones son forzar G o revisar la pose.

## Novelty (proxy debil, limitaciones explicitas)
Max. identidad ventana 8 aa: vs GLP-1(7-36) 3-5/8; vs exendina-4 3-4/8.
NO demuestra diversidad estructural, ni cumplimiento de la regla de la competicion,
ni zero-shot: es solo un proxy de secuencia local.

## Alcance
Gates/RMSD sobre el scaffold miden fidelidad de construccion de la pose, no validan
estas secuencias ni su plegado. Cualquier afirmacion de plegado/funcion requiere
estructura nueva (re-plegado independiente) con autorizacion de GPU.
