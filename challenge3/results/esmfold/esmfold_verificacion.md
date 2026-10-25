# Auditoria independiente ESMFold tanda 25oct — verificador

## Reproducido EXACTO (codigo propio, sobre los PDB entregados)
- 7/7 hashes = declarados (5 PDB + run3.log + env.txt).
- Secuencias de los 5 PDB = FASTAs publicadas (top02/top14/top05/aleatorio36/nativo). nCA 36/36/36/36/30.
- CA-RMSD Kabsch vs refs publicadas: 0.695 / 0.480 / 0.640 / 0.607 / 6.426 A = reportado EXACTO.
- pLDDT medio (B-factors CA): 0.923 / 0.904 / 0.918 / 0.885 / 0.640 = reportado (0.92/0.90/0.92/0.89/0.64).
- Escala pLDDT: B-factors en 0..1 (min 0.42, max 0.94 en los 5 PDB): transformers 5.17.0 escribe
  pLDDT/100. Confirmado; la nota del md es correcta.
- Metrica: CA-RMSD Kabsch + pLDDT, NO TM-score. Correcto el retiro de la etiqueta scTM.

## Provenance y entorno
- env.txt: DGX H200 MIG 1g.18gb, torch 2.14.0+cu126, transformers 5.17.0, fp32, facebook/esmfold_v1.
- DISCREPANCIA MENOR resuelta: el paquete de preparacion afirmaba "transformers 5.20.2 verificado";
  la sesion real uso 5.17.0. El md lo declara; la API infer_pdb existe en ambas. Sin impacto en
  numeros, pero el "VERIFICADO EN CPU" de la preparacion era sobre otra version.
- Log sin secretos (HF sin token, rutas relativas). UUIDs de GPU/MIG en env.txt: identificadores
  de hardware, no credenciales; se conservan como provenance.
- md saneado por el verificador antes de publicar: retirado identificador de mensaje (wamid) del
  parrafo SCOPE. La aprobacion del usuario no necesita detalles del chat en el repo publico.

## Limites (se mantienen)
- Auto-consistencia de UN modelo sobre secuencias AISLADAS: no es evidencia de union, interfaz,
  funcion ni pliegue en contexto del receptor.
- Lectura acotada: el test NO discrimino candidatas de ESTA permutacion en ESTA corrida; un solo
  control no generaliza. Nativo aislado RMSD 6.43 A: prediccion de este protocolo, no prueba de
  dependencia de contexto.
- Boltz-2 u otro modelo complejo: disponibilidad/capacidad NO verificadas; no presentado como listo.
