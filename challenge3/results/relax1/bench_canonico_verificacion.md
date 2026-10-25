# Auditoria independiente: benchmark canonico (preps fijas, masa0, CPU 1 hilo) — verificador, 25 oct

Verificado con codigo propio (parser y metricas re-escritos) y re-ejecucion con OpenMM 8.6.1
(pip, sandbox del verificador). No se ha confiado en ningun resumen.

## 1. Hashes: TODOS EXACTOS
relax_fixed.py, prep_LIBRE6_2/prep_nativo/prep_ala6 (repo), prep_LIBRE6_2_b, out_libre_t1a,
out_libre_t1b (= t1a byte a byte), out_nativo_t1, out_ala6_t1, out_libre_prepb: sha256 = declarados.

## 2. Metricas del benchmark: REPRODUCIDAS EXACTAS (codigo propio)
- C1 nativo:    3->3 choques | dmin 2.46->2.47 | contactos 980->963 | lig 0.325 | receptor 0.000000
- C2 poly-ALA:  4->0 | 2.36->2.93 | 245->249 | 0.336 | 0.000000
- C3 LIBRE6_2: 21->1 | 2.10->2.61 | 601->550 (-51) | 0.863 | 0.000000
Receptor exactamente inmovil (masa0) en los 3 sistemas: RMSD 0.000000 confirmado.

## 3. Robustez H-init: CONFIRMADA y cuantificada
- prep vs prep_b: 2876/2876 pesados con claves y posiciones IDENTICAS; 2837/2837 H con posicion distinta.
- out_libre_prepb (minimo con prep_b): 21->2 choques, lig RMSD 0.833, contactos 601->551, receptor 0. EXACTO.
- RMSD pesados entre minimos (t1a vs prepb): 0.1050 A sobre 2876 atomos. EXACTO.
- Conclusion limitada correcta: direccion robusta (21->1-2), colas NO reproducibles sin preparada fija.

## 4. Causalidad H->deltaE: RE-EJECUTADA, identica
decomp_causal.py re-ejecutado por el verificador: Bond +1838.5, Torsion -575.1, Nonbonded +1180.5,
Angle +25.5; totales 109992 vs 112461 (delta +2469). Numeros identicos al log del investigador.
Alcance: demuestra que H-init cambia E0 via terminos con H en ESTE par de preps; no es prueba general.

## 5. Determinismo CPU: CONFIRMADO, con matiz de entorno NUEVO
- Mis 2 ejecuciones con CpuThreads=1 (misma prep): BYTE-IDENTICAS entre si.
- Mis 2 ejecuciones con CpuThreads=2: DIFIEREN (contactos 553 vs 560, lig RMSD 0.871 vs 0.892).
  Las dos afirmaciones del investigador quedan confirmadas independientemente.
- MATIZ NUEVO (cross-entorno): mi ejecucion 1-hilo (pip 8.6.1, sandbox) != out_libre_t1a del
  investigador: contactos 553 vs 550, lig RMSD 0.893 vs 0.863, dmin 2.60 vs 2.61; RMSD pesados
  entre minimos 0.027 A (ligando 0.089, max 0.533). E0 identico (109992): la trayectoria diverge
  en la minimizacion, no en la entrada. Consecuencia: los valores EXACTOS de cola solo son
  reproducibles con la misma preparada Y el mismo build/entorno; la direccion (21->1) es robusta.

## 6. Cobertura nativo (raw vs prep): VERIFICADA
Emparejamiento por coordenadas (mismo nombre, <1e-3 A): 3371/3371 atomos pesados del raw presentes
en la prep, 0 movidos. Anadidos: 45 = 41 de cadena lateral (CG 10, OD1 5, OD2 4, CD 3, OG 3,
CD1 2, CD2 2, OE1 2, CE/CE1/CE2/CG2/CZ/ND2/NE2/NZ/OE2/OG1 1 c/u) + 4 OXT. Confirmada la cifra de main.
- Conteo nativo pre = 3 choques / dmin 2.46 / contactos 980: IDENTICO medido sobre raw y sobre prep
  (los 45 atomos anadidos no cambian el conteo pre en este caso).
- Etiqueta obligatoria: raw->post NO es parejo en cobertura (post incluye los 45); prep->post SI es
  parejo (mismo conjunto de atomos). No mezclar las dos lecturas. C2/C3: solo 7 OXT excluidos; parejos.

## Limites que se mantienen
Cribado/minimizacion en vacio, una candidata, un protocolo; no es ranking de calidad ni validacion
del diseno; plegado y funcion sin evaluar. Comparaciones A/A' retrospectivas: no reparables
(par original de preparadas no existe). Sin GPU.
