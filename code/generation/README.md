# Generation pipeline (recovered)

Recovered 2026-10-05 from Kaggle notebook `notebook6198e9562d`, **Version 2**
(`scriptVersionId=355016993`, 2026-10-04 01:00:55 CEST) — the most recent
notebook version containing the BoltzGen binder-generation code before the
notebook was switched to controls-only versions (V3–V6).

## `diii_batch_v2.py`

Single-cell Kaggle script that runs one full design batch end to end:

1. Installs `boltzgen`, downloads PDB 6ARU, extracts chain A (EGFR).
2. Writes `patchDIII.yaml`: de novo binder (id B, 50–90 aa) against a 17-residue
   EGFR patch: `369,395,397,400,402,403,406,407,408,427,428,431,446,447,449,451,452`
   (6ARU numbering = UniProt − 27 for this region, verified residue by residue).
3. Asserts the installed boltzgen uses **empty local MSA** (`msa={}`, no
   `msa_server`/`mmseqs`) — no external sequence search.
4. Runs: `boltzgen run patchDIII.yaml --output oc_DIII --protocol
   protein-anything --num_designs 16 --budget 8 --config filtering
   filter_bindingsite=true`
   Filters: `rmsd ≤ 2.5, monomer_rmsd ≤ 2.0, plddt ≥ 75, iptm ≥ 0.6,
   ipae ≤ 10, pae ≤ 10, affinity ≤ 0, complex_rmsd ≤ 2.5, ptm ≥ 0.6`.
5. Post-processes: per-design metrics table, <4 Å contacts vs the patch,
   SHA256 of every output CIF.
6. Archives everything into `DIII_all_outputs.zip` with a sha256 manifest.

Measured runtime for this batch: **2 h 48 m on 2× Tesla T4**, 229 output files.
Outcome of the original run: 16/16 completed all stages, **0/16 passed filters**
(best: patchDIII_15, patchDIII_09, patchDIII_01).

## Reusing it for a new batch

To target a new patch/hypothesis, change only:
- `LABELS` (comma-separated EGFR residue numbers in 6ARU numbering),
- the `--num_designs` / `--budget` flags,
- the `--output` directory name.

Everything else (install, MSA check, filtering config, post-processing,
archiving) is reusable as-is. Do not run two batches from the same notebook
session concurrently; Kaggle sessions can die and lose unsaved results.
