# Auditoria independiente relax1 (verificador, 25 oct)

## El error original, confirmado
Los PDB post-relajacion tienen 7 cadenas (candidata) y 4 (nativo) por cortes TER +
renumeracion de OpenMM. El selector "cadena mas grande = receptor" comparo DOS
FRAGMENTOS DEL RECEPTOR entre si: los numeros "21->0, dmin 8.24" y "3->0" eran falsos.
Retirados por el investigador; confirmado por inspeccion directa de las cadenas.

## Verificacion independiente (mi propia implementacion, mapeo secuencia-explicito)
Hashes de entradas/salidas = env_fijos.txt (sha256 coinciden los 4 archivos).
Mapeo: ligando localizado por SECUENCIA exacta de resnames; receptor por concatenacion
con secuencia verificada identica; atomos pesados comunes por nombre; OXT (capping)
excluidos. Resultados reproducidos EXACTAMENTE:
- CANDIDATA LIBRE6_2: choques 21->1 | dmin 2.10->2.58 | contactos<5A 601->634 |
  binder CA RMSD 1.140 (max 3.24) | receptor CA RMSD 0.881 (max 4.54)
- NATIVO 6X18 R-P: choques 3->6 | dmin 2.46->2.58 | contactos<5A 980->1102 |
  peptido CA RMSD 0.567 (max 2.14) | receptor CA RMSD 0.674 (max 4.38)
Protocolo reportado: OpenMM 8.6.1, amber14, VACIO, PDBFixer+H, restriccion CA k=100
(NO fija receptor), LocalEnergyMinimizer tol=10 maxIter=500, CPU. Sin GPU.
Config coherente con lo observado; NO re-ejecute la minimizacion (OpenMM no instalado
en este entorno): verifico mapeo+metricas sobre sus salidas, no la dinamica.

## Discrepancia NO resuelta: anclajes
El investigador reporta N-term-bolsillo 6.26->5.36 y C-term-ECD 5.24->6.09
("enganche preservado") sin especificar metrica/zona. Mi metrica (distancia minima
CA a CAs de zona, zonas U367-391 / U120-138 del muestreo) da 7.35->5.75 (se acerca)
y 16.65->20.00 (se ALEJA 3.3 A). El mantenimiento de interfaz lo sostiene el conteo
global de contactos (601->634), pero "C-term enganchado" NO queda confirmado con mi
metrica. Pendiente: que el investigador fije zona/metrica exactas de su anclaje.

## Limites (publicar con esto)
- Vacio, amber14, UNA candidata, UN protocolo: no valida diseno ni generaliza.
- El criterio estricto vdW>0.4 marca incluso la interfaz NATIVA tras relajacion
  (3->6): el cero no es referencia calibrada post-relajacion y 1 (candidata) vs 6
  (nativo) NO ordena a la candidata por encima del nativo (ranking retirado, correcto).
- Receptor NO fijo en relax1 (deriva 0.88 A): la absorcion mezcla acomodo del binder
  y movimiento del receptor. relax2 (masa 0) lo separa.

## Revision de protocolo_relax2 (PREPARADO, NO EJECUTADO): pertinente con condiciones
Diseno adecuado: receptor masa 0 (fijo duro), binder libre, controles C1 nativo /
C2 backbone-ALA / C3 candidata, lecturas y criterio de exito PRE-declarados
(choques <= nativo mismo protocolo con RMSD binder <=1 A). Condiciones/anotaciones:
1. Verificar receptor RMSD == 0 exacto como sanity del congelado (masa 0).
2. El criterio de exito es de CRIBADO: superarlo no valida el diseno (falta plegado
   independiente, funcion, solvente); no convertirlo en certificado.
3. Energias solo relativas, nunca afinidad (ya declarado).
4. Variante A' (CA binder k=100) util para acotar deriva; reportar ambas.
CPU local ~1-2 min/sistema, sin GPU: sin objecion tecnica. Luz verde la da main.
