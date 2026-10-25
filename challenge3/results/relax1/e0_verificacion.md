# Verificacion de la resolucion E0 (verificador, 25 oct)

## Verificado independientemente
- Hashes de las 3 preps canonicas = declarados en resolv_e0.md.
- Atomos PESADOS: comparacion por conjunto (resname, nombre, coords) raw -> prep:
  0 atomos raw movidos/ausentes en los 3 sistemas. La preparacion NO toca ningun
  atomo pesado existente.
- Anadidos por PDBFixer (documentados): OXT de capping en C3/C2 (7 cada uno); en
  NATIVO ademas ~45 atomos de cadena lateral que FALTABAN en el 6X18 recortado
  (CG/OD/OE/CD/OG/NE/CZ/NZ...): el control nativo de los relax incluye sidechains
  RECONSTRUIDOS por PDBFixer, no todos cristalograficos. Se etiqueta aqui.
- H: 2837/2837 (C3), mismos conteos por cadena entre salidas A y A'.

## Consistente pero NO demostrado aqui (limites honestos)
- La no-determinismo de addMissingHydrogens entre procesos y que el delta E0 salga
  de NonbondedForce de H: consistente con todo lo observado, pero el PAR original de
  preps (la de A y la de A') no se adjunto: no puedo demostrar la comparacion por
  pares. Pedir via main si se quiere cerrado.
- La descomposicion por grupos (CustomExternal=0 en el punto inicial) requiere OpenMM:
  no verificada aqui.

## Correccion de alcance (incorporada, clave)
Con preps distintas, A y A' NO son una comparacion aislada del protocolo: los H
distintos cambian fuerzas y el minimo heavy. Por tanto la diferencia 21->2 (A') vs
21->3 (A) NO es atribuible solo al esquema de restriccion: mezcla ruido de preparacion.
Las metricas intra-corrida (pre vs post de la MISMA corrida) usan los mismos pesados y
siguen siendo validas. El benchmark con preps canonicas fijas (autorizado, CPU,
sin re-anadir H, sin 16/GPU) es lo que cerrara la comparacion; ejecuta el investigador,
audito yo.
