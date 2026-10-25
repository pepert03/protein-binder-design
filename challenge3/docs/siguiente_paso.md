> CORRECCION 2026-10-25 (verificacion independiente del escritor): scaffold_v3 FALLA la
> verificacion quimica (ver docs/scaffolds.md): longitudes correctas pero angulos/diedros
> invalidos, traza CA no coincide con PARAMS_V3. La propuesta de abajo queda BLOQUEADA
> hasta existir un backbone valido. No ejecutada; la estimacion "<30 min GPU" NO esta
> comprobada; cualquier sesion GPU requiere autorizacion explicita (no concedida).

# Siguiente paso: viabilidad concreta y coste (sin GPU, sin campaña nueva)

Fecha: 2026-10-25. Estado: scaffold_v3.pdb = hélice ideal 36 aa con geometría REAL de backbone
(0 enlaces fuera de tolerancia Engh-Huber), pose encontrada SOLO con coordenadas del receptor,
0 choques backbone-receptor, gate de orientación PASS. NO es un diseño: no tiene secuencia.

## Bloqueo real (como señaló main): método que conserve el frame Y genere secuencia Y se valide independientemente

Propuesta concreta en 2 etapas:
1. Diseño de secuencia sobre scaffold_v3: ProteinMPNN (fixbb) sobre la cadena B, fijando la
   pose. Restricción local: B3 está ajustada (CB a 1.97 A del receptor) → forzar residuo
   pequeño (G/A/S) en B3 y vigilar B1-B8 (bolsillo profundo).
2. Validación independiente: re-plegar las secuencias con ColabDesign/AF2 (modelo DISTINTO al
   que generó la secuencia) y exigir: RMSD binder vs scaffold_v3 < 2 A + orient_gate PASS +
   gate_pose 0 choques + ensemble_rmsd binder < 2 A. Solo entonces "compatible con modo de
   unión nativo" (nunca "agonista").

## Coste y recursos (todo gratuito conocido)

- Etapa 1 (ProteinMPNN): CPU local o Colab free; minutos para ~50 secuencias. Coste 0, sin GPU.
- Etapa 2 (validación AF2/ColabDesign): requiere GPU. Opciones gratuitas conocidas:
  (a) Colab free T4: ~10-30 min por diseño, factible para 5-10 diseños pero lento y con
      límites de sesión; (b) DGX (sesión nueva): <1 min por diseño en H200 MIG — pero eso ES
      una sesión GPU nueva y requiere autorización explícita; estimación total <30 min GPU
      para validar ~10 candidatos.
- Nada de esto se ejecuta sin autorización; esto es solo el plan y su coste.

## Criterio de éxito de la siguiente prueba (propuesto)

- PASS si >=1 secuencia re-plegada cumple todos los gates con RMSD al scaffold <2 A.
- Si ninguna lo cumple: el frame es geométricamente válido pero no diseñable/plegable →
  resultado negativo honesto, documentado, sin más GPU.

## Provenance scaffold_v3

Búsqueda numpy local (exporta_scaffold3.py, seed 23, determinista). PARAMS_V3:
anclaje [122.732, 126.488, 107.399], eje [-0.0849, 0.9769, 0.1960], roll -0.00129.
Verificación del artefacto pura Python: bb_check.py (química) + orient_gate.py v2 (G1-G4).
