# Viabilidad local del rediseño por extremos (Challenge 3, GLP1R) — 2026-10-25

Análisis local SIN GPU sobre artefactos verificados del repo. Cada cifra fue
reproducida de forma independiente antes de publicar (salvo lo marcado "no
reproducido"). Esto NO prueba la viabilidad del rediseño completo.

## 0. Marco de reglas (verificado 2026-10-25 en la página oficial)

Las reglas oficiales exigen diseños "de novo and zero-shot": no se puede tomar
un binder existente (p. ej. el GLP-1 de 6X18) como punto de partida ni
modificarlo; sí se permite usar binders conocidos para calibrar filtros y
métricas de diseño. Fuente:
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/glp1r (FAQ 8).
Consecuencia: todo uso del péptido nativo de 6X18 en este documento es
CALIBRACIÓN geométrica, nunca base de diseño. La vía de diseño es target-only.

## 1. Gates implementados y validados (reproducidos)

- `scripts/ensemble_rmsd.py`: RMSD entre MODELs con superposición Kabsch sobre
  CA del receptor. Reproduce: g3_1 = 0.49 A global / 1.13 A binder;
  g3b_1 = 1.62 A global / 4.93 A binder. (El 1.58 A citado antes para g3b_1
  queda supersedido: era de una variante de ajuste; diferencia 0.04 A.)
- `scripts/orient_gate.py`: gate G1-G4. G4 usa pares de ATOMOS <2.0 A con la
  misma semántica all-atom de gate_pose.py y reproduce exactamente
  g3_1 M1=27/M2=24 y g3b_1 M1=15/M2=54. Ambos diseños conocidos dan FAIL
  global, como corresponde a un gate que discrimina.

## 2. Calibración geométrica sobre 6X18 (marco nativo = SOLO calibración)

`scripts/viabilidad_scaffold.py` sobre 6X18_RP.pdb (R=receptor, P=GLP-1 7-36),
reproducido:
- Corredor nativo: 0/30 CA del péptido a <2.0 A del receptor; mínimo 2.95 A en
  el tramo insertado. El corredor admite un tubo de backbone.
- Span N->C del péptido nativo: 40.6 A en 30 res; una hélice ideal cubre ese
  span con ~28 res (un binder de 36 aa es axialmente suficiente).
- 0 choques de backbone del marco nativo contra el receptor.
Límite: esto muestra que el corredor nativo admite un backbone; NO implica que
un diseño de novo vaya a ocuparlo.

## 3. Muestreo target-only (pseudo-backbone, sin sidechains)

- `scripts/muestreo_helice.py` (4000 poses, solo coordenadas del receptor;
  anclaje N en U367-391, cono 40 grados, giro libre), reproducido: mejor pose
  = 4 choques, extremo C a 2.6 A del ECD. Ninguna pose con 0 choques.
- Hélice recta sobre eje de centroides (helice_target_only.py): 126 choques
  reportados; script no compartido en este lote — NO reproducido aquí.
- Muestreo con bisagra (dos segmentos, 2500 poses): mejor = 12 choques,
  reportado por el investigador; script aún no compartido — NO reproducido
  aquí. Tampoco demuestra pose cero.

LÍMITES DUROS de estas cifras:
- El muestreo usa un pseudo-backbone aproximado (CA + N/C estimados) contra
  átomos pesados del receptor, SIN cadenas laterales del binder. Sus conteos
  NO son comparables con el gate all-atom de los diseños (27/24, 15/54).
- Backbone sin choques no demostraría sidechains compatibles, y nada aquí
  predice unión, activación ni agonismo.

## 4. Afirmaciones RETIRADAS (no citar como vigentes)

- "Los ~4 choques residuales los absorbe una etapa de diseño con backbone
  flexible / relajación": hipótesis NO demostrada. RETIRADA.
- "15-54 choques son irrecuperables por relajación": hipótesis NO demostrada.
  RETIRADA. Además la comparación 4 vs 15-54 no es válida: miden cosas
  distintas (pseudo-backbone vs all-atom, ver §3).
- "VIABLE" como conclusión del rediseño completo: NO probado. Lo demostrado es
  solo: existe corredor geométrico en la estructura nativa (calibración), y un
  muestreo simple de varilla rígida no encuentra pose perfecta.

## 5. Siguientes pasos locales propuestos (sin GPU)

1. Muestreo variable (dos segmentos con bisagra y parámetros variables) hasta
   0 choques o demostrar lo contrario. En preparación.
2. Si aparece una pose: exportar PDB del scaffold y dejar listo el input de
   fixbb/ProteinMPNN (sin ejecutar).
3. Integrar orient_gate.py + ensemble_rmsd.py como gate post-trayectoria del
   driver cuando se autorice campaña (vía escritor único Umbrel).
