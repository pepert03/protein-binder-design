# Propuesta: condicionamiento por extremo sin copiar ligando conocido

Objetivo: que el extremo N del péptido contacte el bolsillo de activación (U367-391) y el C el ECD,
sin copiar motivos de GLP-1/exendina y sin asumir que un hotspot global produce orientación.

**Estado: PROPUESTA DE DISEÑO. Marco de reglas verificado (FAQ 8 oficial): de novo/zero-shot; el marco del ligando nativo solo vale para calibrar filtros, nunca como base de diseno.**

**Detalle: NO ejecutada. Ninguna de estas opciones ha corrido; el agonismo NO
está probado para ningún diseño de este repositorio. Requeriría una nueva campaña autorizada (GPU).**

## Opciones (de menor a mayor cambio de pipeline)

1. Scaffolding geométrico (recomendada): colocar una hélice ideal sobre coordenadas del TARGET
   (no del ligando): extremo N en la boca del bolsillo U367-391, extremo C sobre la superficie ECD
   (~U120-138). Diseñar secuencia sobre ese andamio (p.ej. fixbb/ProteinMPNN) y validar con AF2.
   La orientación queda fijada por construcción; solo se usan coordenadas del receptor.
2. Sesgo de estructura secundaria: una hélice de 30-36 aa es una varilla de ~45-54 Å; si se exige
   contenido helicoidal (loss de SS) y contacto con ECD, la varilla solo puede alcanzar el bundle
   con el extremo opuesto. Inductivo, no garantiza orientación; gate post-hoc obligatorio.
3. Hotspots dirigidos por extremo: hotspot del bolsillo de activación SOLO evaluado contra
   residuos 1-12 del binder y hotspot ECD contra 27-36 (requiere loss por segmentos, modificación
   del driver; ColabDesign estándar no lo trae). No copia ligando, pero asume geometría objetivo.
4. Sin cambio de driver: gate de orientación post-hoc (scripts/orient.py) como criterio de
   aceptación junto al gate estérico; descarta inversas sin rediseñar. Barato pero desperdicia
   GPU en trayectorias condenadas — por eso se recomienda combinarlo con (1).

## Gates de aceptación propuestos (por trayectoria, antes de cualquier relajación)

- Estérico: 0 pares intercadena <2 Å (scripts/gate_pose.py).
- Orientación: ≥3 residuos de 1-12 con contacto ≤6 Å a U367-391; ≥3 residuos de 27-36 con
  contacto ≤6 Å a ECD; 0 residuos de 27-36 con contacto ≤3 Å al core TM.
- Ensemble: binder RMSD <2 Å entre miembros (scripts/ensemble_rmsd.py).
- Presentación honesta: "compatible con el modo de unión nativo" — nunca "agonista".
