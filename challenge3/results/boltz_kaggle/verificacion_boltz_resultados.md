# Verificación de resultados Boltz-2 (Kaggle) — verificador independiente, 25 oct

Tanda autorizada por el usuario el 25oct 18:15 CET: 1 candidata (top02) + péptido nativo GLP-1 +
control permutado (aleatorio36), en Kaggle con cuota gratuita, sin submission al concurso.
Ejecutor: agente investigador. Este documento es la auditoría independiente: todos los números
recomputados desde los artefactos, ninguno copiado de resúmenes.

## 1. Identidad de inputs — VERIFICADA
- sha1 de los YAML ejecutados: top02 58e5b22da917b798779237e5c8fdb4b972b7aca2,
  nativo a353048878351afb231c964c1fb921c6ded1f570, aleatorio36 f7f36d775a64bab4bf528e32336511ff25e2584f.
- Los YAML ejecutados difieren de los DRAFT del repo (results/boltz_prep) SOLO en las cabeceras
  de comentario (copia anterior al saneo de cabeceras). Las SECUENCIAS son idénticas:
  cadena A = GLP1R UniProt P43220 24-463 (440 aa, match exacto verificado contra el SQ de UniProt);
  cadena B = binder exacto de cada caso (top02 36 aa; nativo 30 aa = GLP-1 7-36; aleatorio36 36 aa).
- 33/33 hashes del manifiesto out_hashes.json verificados byte a byte.

## 2. Métricas de confianza — REPRODUCIDAS (escala 0..1)
| caso | ipTM | confidence | pTM | complex pLDDT |
|---|---|---|---|---|
| top02 | 0.9067 | 0.7378 | 0.8580 | 0.6956 |
| nativo | 0.3662 | 0.6513 | 0.7635 | 0.7226 |
| aleatorio36 | 0.1663 | 0.6268 | 0.7466 | 0.7420 |

complex_plddt del JSON == media del array pLDDT (npz). pLDDT medio de CA por cadena:
top02 receptor 0.702 / binder 0.615; nativo receptor 0.744 / péptido 0.406; aleatorio36 receptor 0.753 / péptido 0.607.

## 3. Log de ejecución — SIN ERRORES
3/3 predicciones "Number of failed examples: 0" (~2:45-2:57 min cada una, GPU T4, 1 muestra de
difusión, recycling 3, sampling 200). Solo avisos benignos (deprecaciones de pytorch_lightning,
nota de versión de checkpoint Lightning). Sin MSA (msa: empty) y sin templates aportados.

## 4. Análisis de pose (propio; alineamientos Kabsch sobre CA comunes)
Referencia: EM 6X18 (receptor 384 CA observados UP 29-423 + GLP-1 7-36) y complejo_v6 (diseño).
- Receptor predicho vs EM: top02 6.38 Å, nativo 8.58 Å (384 CA). Aislando el CORE TM
  (166 CA de los 7 segmentos FT): top02 2.59 Å, nativo 2.57 Å. El desvío global viene del ECD,
  loops y cola citoplasmática, no del bolsillo TM.
- Binder top02 vs pose de diseño (complejo_v6, mismo receptor de referencia): 5.20 Å CA-RMSD.
- Binder top02 (36 aa) vs péptido nativo EM (30 aa): NO se usa RMSD (longitudes distintas, sin mapa
  1:1 justificado). Medida: distancia de cada CA del binder al CA nativo más cercano, con el receptor
  alineado por 166 CA del core TM. Media global 10.5 Å, pero por mitades: res 1-15 = 3.9 Å (el
  N-terminal del binder sigue la trayectoria nativa, varios a ~2 Å), res 16-36 = 15.2 Å (la mitad
  C-terminal diverge hasta 25 Å). Mapeo residuo a residuo completo en analiza_result.json/stdout.
- Péptido nativo predicho vs su pose EM real: 34.3 Å — el modelo NO reproduce el modo de unión nativo.
- Contactos de interfaz (<4.5 Å, receptor en numeración UniProt):
  top02: 45 residuos en ECD (29-148) y bolsillo TM/ECL (197-217, 233-237, 298-310, 364, 375-388)
  — región ortostérica correcta, pose desplazada.
  nativo: 31 residuos, mayoría ECD (67-148) + 209-222 — no el bolsillo TM profundo.
  aleatorio36: 28 residuos, todos en región citoplasmática/cola (254-463) — sitio distinto.

## 5. Lectura honesta
- ipTM 0.907 de top02 es la confianza DEL MODELO en SU pose, con 1 muestra, sin MSA, sin templates,
  sin contexto de membrana/Gs. La pose favorece la región ortostérica pero desplazada ~5 Å del diseño.
- El control nativo bajo (0.366) NO autentica agonismo ni unión del nativo: el modelo tampoco
  reproduce la pose nativa experimental (34 Å).
- Esto NO es un diseño validado. Es una comprobación de coherencia con resultado favorable y acotado.
- Figura comparativa: comparativa_top02_vs_6X18.png (trazas CA, alineada por 166 CA del core TM;
  gris = receptor EM 6X18, azul = receptor predicho, rojo = binder top02 36 aa, verde = GLP-1 nativo
  EM 30 aa; inspección visual realizada).

## 6. Incidencias de esta verificación (transparencia)
- Mi primera pasada de alineamiento tenía un bug en la convención de Kabsch (RMSD receptor 52 Å,
  imposible); detectado con un test sintético, corregido y repetido todo el análisis.
- 6X18_RP.pdb cadena R incluye 20 aguas numeradas 501-520: excluidas del análisis (solo proteína 24-463).
- Escaneo de datos privados sobre todo el paquete: sin wamids, IDs de agentes/chats, cuentas, rutas de
  usuario ni secretos (los únicos matches eran texto de ayuda de boltz y nombres de hiperparámetros).
- Script de análisis incluido (analiza.py); resultados crudos en analiza_result.json.
