# Análisis de orientación extremo-a-extremo — g3_1 y g3b_1 (Challenge 3, GLP1R)

Fecha: 2026-10-25. Artefactos: `results/traj1/g3_1.pdb`, `results/traj2/g3b_1.pdb` (md5-verificados en el repo).
Script: `scripts/orient.py`. Uso: `python3 orient.py <complejo.pdb>`.
Método: por MODEL (2 miembros de ensemble por PDB de ColabDesign, MISMO mejor paso — una candidata
por trayectoria), cada residuo del binder (cadena B) se clasifica por el residuo de receptor (cadena A)
más próximo ≤6 Å; zonas según UniProt P43220: ECD ≤139, ECL1 202-227 / ECL2 291-305 / ECL3 371-383,
TM 140-404 resto. Mapeo resseq→UniProt idéntico al de gate_pose.py, verificado fragmento a fragmento
contra la secuencia real P43220 (6/6 exactas). Altloc solo ' '/A; sin duplicados.
Salidas reproducidas en el clon Umbrel: `results/traj1/orient_g3_1.txt`, `results/traj2/orient_g3b_1.txt`.

## 1. Orientación observada (hechos, reproducidos)

| Diseño | N-terminal (1-12) | C-terminal (27-36) |
|---|---|---|
| g3_1 (ambos MODELs) | ECD + ECL loops (U68-138, U210-214) | Enterrado en core TM: contactos U237-U387 con solapes imposibles 0.6-2.1 Å (mínimo 0.6 Å en U314; también 0.9-1.8 Å en U368-369) |
| g3b_1 (ambos MODELs) | ECD + ECL (U29-141, U202-210; también U299 en MODEL 2) | Solo 2 residuos con contactos TM débiles (3.9-5.8 Å; cola 30-36 sin contacto ≤6 Å); choques en segmento medio (B24), ver gate |

Hallazgo geométrico: en 2/2 trayectorias la orientación es la INVERSA al modo de unión nativo
(referencia: two-domain model de GPCR clase B1 — el ECD ancla el extremo C del péptido y el N se
inserta en el core TM; esa segunda interacción es la esencial para activación. Fuente textual:
https://pmc.ncbi.nlm.nih.gov/articles/PMC3415635/ ; coherente con 6X18 y Wootten 2017).

## 2. Hipótesis causal (NO probada)

Una lectura posible es que optimizar i_pae/ptm con hotspot global no induce orientación nativa.
NO está demostrado: es un hallazgo en DOS trayectorias con UN objetivo de optimización.
Alternativas no descartadas: inicialización, MSA, ausencia de membrana, contexto citoplasmático
ausente, longitud/secuencia. El receptor completo tiene diferencias de contexto reales (membrana,
ICLs, proteína G) aunque la topología verificada indica que no explican directamente los choques
extracelulares observados.

## 3. Límites (leer antes de interpretar)

- La orientación nativa es un criterio de priorización de ESTA estrategia, no un requisito
  universal probado para todo agonista de novo.
- Ninguna métrica de esta pipeline (i_pae, plddt, gates de choque u orientación) predice
  activación ni agonismo. Ambas trayectorias además FALLAN el gate estérico.
