# Campaña zero-shot Ch3 GLP1R (409 aa) - 2026-10-26T07:37:24.051265Z
Diseños de novo (ColabDesign AF2-binder, sin binder de partida) sobre target oficial UniProt P43220 24-432 (409 aa, sha1 fasta 689647a8c8b799bdf43c9b753e1e83ff5446b1f7).
Filtro Boltz-2 2.2.1 (--recycling_steps 3 --sampling_steps 200 --diffusion_samples 1), semilla primaria + semilla 1234.
Criterios: ipTM, estabilidad entre semillas, contactos en epítopo GLP-1 (docs/epitope_contacts.tsv, 57 contactos; hotspots ECD/TM en scripts/contacts409.py).
Controles: nativo GLP-1 (30 aa), aleatorio36.
TODO ES PREDICCIÓN COMPUTACIONAL; no medición ni evidencia de agonismo.
Reproducir: scripts/ contiene generación y filtro; env/ el entorno; MANIFEST.sha256 los hashes.
