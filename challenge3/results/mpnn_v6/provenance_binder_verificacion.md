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

## PENDIENTE (no cerrado)
- complejo_v6.pdb sha256 5e8c41ec... (archivo sobre el que corrio pack_eval_v6) NO esta en el repo.
- El "packed" declarado con sha256 a32d4455... tampoco identificado en el repo.
- Pedido al investigador via main: subir ambos (o nombrar rutas exactas) para completar la cadena
  scaffold -> complejo -> packed -> relajado.

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
- DEFECTO: el md del paquete corregido conserva el PROVENANCE FLAG retirado (148-183/1.96 A) y un
  parrafo "PERMISO" que cita autorizacion: un documento no es evidencia de permiso. El permiso real
  esta en los mensajes del usuario (investigacion/pruebas/pushes SI; GPU NO). Corregir el md antes
  de publicar el paquete.
