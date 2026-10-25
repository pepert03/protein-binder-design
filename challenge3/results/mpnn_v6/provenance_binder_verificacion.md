# Resolucion de provenance binder scaffold_v6/complejo — verificador, 25 oct

## El flag "resseq 148-183, 1.96 A" es un error de parseo: REPRODUCIDO y explicado
- scaffold_v6.pdb tiene 2 cadenas: A = receptor (resseq 1..376, 316 res) y B = binder (resseq 1..36).
- Mi barrido de ventanas de 36 CA: la mejor ventana de la cadena A (RECEPTOR) es resseq 148-183
  con Kabsch RMSD 1.961 A — el numero exacto del investigador. La cadena B da 0.000 A.
- Causa: su extractor agrupaba por TER (scaffold_v6.pdb no escribe TER) y no vio la cadena B;
  el "mejor segmento" era receptor. Investigador lo ha retirado; queda aqui la evidencia.

## Identidad verificada (codigo propio)
- scaffold_v6 sha256 ccd1db10...d09a (repo) y cand2_complex sha256 73994371...2880 (repo): hashes declarados.
- Cadena B scaffold = 36xALA (proxy poly-ALA); cadena B cand2 = secuencia EXACTA del sample 2 de MPNN v6
  (SAESLLTEGLGEALRARNARLAREAAEAAELAERAA): cand2_complex = scaffold_v6 empaquetado con sample 2.
- Backbone binder (N/CA/C/O): 144/144 atomos IDENTICOS (max 0.0 A). Drift FASPR de backbone: 0 confirmado.
- Backbone receptor (N/CA/C/O): identico entre ambos (max 0.0 A).
- Receptor, cadenas laterales: 11 residuos (46,48,51,69,75,85,86,87,94,101,108) van RECORTADOS a CB
  en scaffold_v6 y completos en cand2_complex (33 atomos CG+ ausentes; CB desplazado <=0.067 A).
  Ninguno de esos residuos es pareja de choque en pack_v6_results (A120/A177/A182/A205), pero los
  contactos d<5 A y cualquier metrica de interfaz que use esas cadenas laterales NO son comparables
  entre los dos archivos. Etiquetado.

## Cadena CERRADA (25 oct, segunda entrega)
- complejo_v6.pdb sha256 5e8c41ec...130c VERIFICADO: receptor identico al de scaffold_v6
  (2575/2575 atomos, max 0.0 A, mismas 11 sidechains recortadas); binder poly-ALA sin CB;
  backbone binder identico (144/144, 0.0 A).
- LIBRE6_2.pdb ("packed") sha256 a32d4455...e6a84 VERIFICADO: receptor completo (2608) IDENTICO
  al de cand2_complex; binder 261/261 atomos IDENTICOS a cand2_complex (max 0.0 A).
  cand2_complex.pdb ES LIBRE6_2.pdb regrabado (solo difieren REMARK/TER/END).
- Cadena scaffold -> complejo -> packed -> relajado COMPLETA y con identidad 0.0 A en cada salto.
- Choques de LIBRE6_2 medidos con codigo propio: 21 / dmin 2.10 = fila sample 2 de pack_v6_results.
- ABIERTO menor: poly-ALA proxy mide 3 choques (dmin 2.50) sobre los artefactos comprometidos
  (con y sin CB); el 4 reportado en analisis_v6.md proviene de archivos del bucle no preservados.
  Discrepancia 3<->4 sin resolver, sin impacto: ambos <= umbral exploratorio 5 y << 21.

## Regresion anadida
challenge3/scripts/test_binder_ref.py: exige cadena B unica con 36 CA e identidad de backbone
scaffold<->complejo; FALLA explicito si falta B, si no son 36 CA, o si difiere el backbone.
Sin fallback por ventana/RMSD (el patron que produjo el falso 148-183).

## Paquete ESMFold (auditado, NO ejecutado; GPU sin autorizar)
- ref_binder_v6.pdb == cand2_complex cadena B (max 0.0 A). ref_nativo.pdb == nativo_RP cadena P
  (30 CA, resseq 7..36, max 0.0 A). nativo_glp1.fa == secuencia de la cadena P. aleatorio36 =
  permutacion de top02 con composicion identica. Referencias CORRECTAS.
- Seleccion v6_top02/05/14 = samples MPNN 2/14/5 = los 3 MENOS choques de tabla_v6 (21/23/24):
  criterio razonable pero el nombre "topNN" es el ID de sample, no el rango; el mejor por score
  MPNN (sample 10) NO esta seleccionado. Sin impacto cientifico; clarificar nombre.
- Script corregido: agrupa por ID de cadena y exige exactamente 1 cadena con n_ca CA (sin fallback
  TER/RMSD). Correcto.
- Tercera entrega: FLAG retirado y criterio de seleccion explicito (samples 2/14/5 por 21/23/24
  choques). Sanitizado por el verificador antes de publicar: eliminados identificadores de mensajes
  y referencias a coordinacion interna del parrafo PERMISO (la autorizacion no necesita detalles
  del chat; un documento nunca es evidencia de permiso). Estado: PREPARADO, NO EJECUTADO; GPU sin
  autorizar; estimaciones de recursos DGX sin verificar.
