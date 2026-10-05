# pH-Selective EGFR Binder Design — Research Report

Design of a conditional (pH-selective) EGFR binder for the ProteinBase challenge
(anthropic-adaptyv-2026, Challenge 01: EGFR).
Repository: `ph-selective-egfr-binder`. Author: José Ridruejo.

This is a living document. New findings are appended to the **Running log**
(newest first) and the summary sections are updated accordingly.
Companion files: `docs/STATUS.md` (reconstruction of Oct 2–4 work from the
transcript), `code/kaggle_v6/` (recovered, hash-verified Version6 code).

**Deadline: Tuesday 2026-10-06, 23:59 AoE = 2026-10-07 13:59 CEST (Madrid).**
Nothing has been submitted. Nothing is submitted without José's explicit go-ahead.

---

## 1. The challenge

- Official page: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
- Task: design a **conditional EGFR binder** against the human EGFR extracellular
  region. Domain III is the recommended functional epitope (targeted by
  cetuximab/panitumumab).
- Ranking of objectives (in this order):
  1. **pH-selective binding**: binds at pH 6.5, no detectable binding at pH 7.4.
  2. **Mouse cross-reactivity**: the same sequence binds human and mouse EGFR.
  3. **Affinity** against human EGFR.
- Submission: CSV ordered by the submitter's own ranking; columns `name`,
  `sequence`, `molecule_class` (`protein` | `nanobody` | `scfv` | `fab_kappa` |
  `fab_lambda`). De novo designs only, 10–250 aa. Up to 40 designs (Track 1) or
  20 (Tracks 2/3). The methodology description is made **public** on submission.
  Only one submission batch counts (the last, or the one explicitly designated).
- Winners are announced per category; a weak but clearly pH-sensitive binder may
  beat a strong non-selective one. Method novelty, diversity and quality are
  valued over filling the quota.

## 2. Key discoveries so far

### 2.1 Assay conditions (Slack, Tudor/Adaptyv, Oct 4)

- Binding measured primarily by **SPR**, with **BLI** to validate hits.
  Our designs are immobilized; EGFR flows over them.
- Human EGFR tested at **pH 6.5 and 7.4**; **mouse EGFR only at pH 6.5**.
  Implication: a genuine pH-switch (off at 7.4) can still demonstrate mouse
  cross-reactivity, since mouse is measured at the acidic pH.
- Targets are full ectodomains, tethered: **Human Met 1–Ser 645**,
  **Mouse Met 1–Ser 647**, supplied by Sino Biological
  (human-egfr-10001-h08h, mouse-egfr-51091-m08h).
- Receptors are produced in **HEK293 → glycosylated**. Glycan-clash alerts from
  modeling are relevant to the real assay, not just the models.
- Cross-reactivity is framed as KD(mouse)/KD(human) ≈ 1.
- Sources: https://proteinbase.slack.com/archives/C0C4VEG57HU/p1791137831013399
  (Tudor FAQ), earlier assay thread
  https://proteinbase.slack.com/archives/C0C4VEG57HU/p1790724749287889?thread_ts=1790637669.096089

### 2.2 Novelty requirement (Slack, Oct 2–5)

- **Novelty Level 3** required, checked at submission; both sequence *and*
  structural comparison against known proteins. Local similarity checks are not
  a substitute. Reference: Adaptyv's novelty blog post (adaptyvbio.com).
- The official novelty pipeline code was promised for Oct 5. As of Oct 5 ~18:05
  CEST, **no verified publication** in the Slack workspace (checked read-only;
  thread replies could not be opened). Still pending.
- Affibodies do not get antibody-style novelty treatment in this challenge.

### 2.3 Design pipeline and batches (Oct 2–4, previous AI "Instinct")

Pipeline: **BoltzGen** (de novo generation + sequence design) → **Boltz 2**
refolding (usually 5 samples; metrics ipTM, PAE, pLDDT) → **Chai-1** orthogonal
pose check → official **DockQ** for interface comparison. Internal filters for a
candidate: complex RMSD (refolded vs generated) ≤ 2.5 Å, contact with the target
patch, no interface contacts < 2.5 Å (clashes), interface confidence judged per
case, partial novelty checks.

| Date | Batch | Designs | Passed | Note |
|---|---|---|---|---|
| Oct 2 | test + numbering fix + patch 2 | 10 | 0 | Low interface confidence; a residue-numbering bug was found and fixed |
| Oct 2 | relaunch | 4 | 0 | 1 touched the patch (ipTM 0.665) but had a too-close contact |
| Oct 2 | evening | 8 | 2 | Best ipTM 0.81 |
| Oct 3 | 16-design | 16 | 2 | 1 convincing (ipTM 0.713, patch contacted, no clashes) |
| Oct 3 | focused 16-design | 16 | 3 | 1 with better signals, 2 weaker |
| Oct 3 | 2 fixed histidines | 16 | 0 | All failed RMSD; His conserved in sequence but not constrained to interface |
| Oct 4 | new DIII patch | 16 | 0 | 15/16 passed the binding-site filter; all failed complex RMSD (80 refold samples checked) |

**Total: 86 designs, 7 passed internal filters (~8%).** Only 4 of the 7 have the
best interface signals. Naming in the transcript is inconsistent (cC_6/cC6);
two passing designs are never named [uncertain].

### 2.4 Candidate assessment (all still predictions, no wet-lab evidence)

- **cC_6** (batch C2): most pose-consistent across Chai-1/Boltz-2; places a Glu
  near receptor **H433** (pH hypothesis). **Risks: high structural-novelty risk**
  (local match to part of a PDB protein); glycan clash is
  conformation-dependent, unresolved.
- **cC_12** (batch C3): Glu near H433 in Chai-1 models; pose consistent in only
  1 of 3 Boltz-2 models.
- **cC4_00**: better novelty signal than cC6; position varies between Chai-1 and
  Boltz-2.
- **cC4_15**: avoids the N352 glycan in the tested context; issues when full
  EGFR domains are added.
- **cC5_00**: strong clash with the first EGFR-attached glycan in two alignments.
- **2 unnamed**: passed filters with weaker signals [uncertain].

Cross-model caveats: Boltz-2 scores **shuffled-sequence controls** highly too
(its score is not discriminative here); Chai-1 separates designs from shuffled
controls in score and position, but that is not independent validation.
**No candidate has any affinity, pH-selectivity, or mouse-binding evidence.**

### 2.5 Controls: what the pipeline can and cannot do

- **Fab (cetuximab mutant, PDB 6ARU)**: recovers 88–92% of native interface
  contacts, interface RMSD ~0.61 Å, **DockQ 0.72–0.84**. The pipeline *can*
  recover a real binder pose. (Full Fab + 5 simultaneous samples OOMs; 1 sample
  fits.)
- **EGF (PDB 1NQL)**: DockQ 0.048/0.014, complex RMSD 11.98–24.06 Å — but the
  input prep had **dropped the 3 disulfide connectivities**, and 1NQL is an
  autoinhibited low-affinity state crystallized at very acidic pH, not the
  active DIII complex. Not a valid negative control of the pipeline.
- EGF replicas show **high seed-to-seed variability** (2 of 5 seeds reach DockQ
  0.43–0.46); one sample is not enough to judge. Recorded seeds do not control
  all internal randomness.

### 2.6 pH-selectivity hypothesis research

- Literature lead: an engineered anti-EGFR antibody placing **Asp/Glu facing
  receptor H433 (domain III)** gained binding in acidic medium —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC9703009/ — but the final construct
  still bound at pH 7.4. Lead, not a solution.
- Second lead: binder histidine facing a receptor Glu for acidic binding —
  https://www.tandfonline.com/doi/full/10.1080/19420862.2019.1682866
- In Chai-1 models, **cC_6 and cC_12 place a Glu near H433** (checked against the
  active-state structure 1IVO). Weakened: shuffled controls also contact H433,
  and contact positions differ between Chai-1 and Boltz-2.
- **Proton-PottsMPNN** (Oct 4, flagged by Amir/organization as highly relevant
  for pH-sensitive design): paper
  https://www.biorxiv.org/content/10.64898/2026.09.30.755438v1.full.pdf, code
  https://github.com/christian-creator/ProtonPottsMPNN. Under active evaluation
  (see running log).
- Caution from literature review: adding histidines can produce the *opposite*
  pH effect to the one wanted; the mechanism needs grounding before more compute.

### 2.7 Compute characterization

- **Kaggle** (jridruejotunon@gmail.com): 2× Tesla T4 (15 GB each); design jobs
  use ~9–10 GB/GPU; ~50 min per 4-design batch. Sessions can die and lose
  unsaved results. Quota on Oct 5: 19:21/30 h used (~10.5 h left).
- **Comillas DGX** (https://hpc.comillas.edu, Open OnDemand, no VPN):
  student profile = **1× 16 GiB MIG slice of an H200, sessions up to 8 h**
  (email said 24 h; portal shows 8 h). Measured: 4 CPU cores, **~5.86 GiB RAM
  per job (binding limit — not GPU memory)**, Python 3.10.12, PyTorch 2.14.0,
  CUDA 12.6, internet OK. EGF prediction 2:56 (DGX) vs 6:53 (Kaggle). Storage
  is **not backed up (RAID 0)** — always copy results off.

### 2.8 Code recovery and validation (Oct 5, Muse)

- Recovered the private Kaggle dataset `proteinbase-codigo-version6`
  (https://www.kaggle.com/datasets/peperidruejo/proteinbase-codigo-version6):
  9 files, in `code/kaggle_v6/`.
- **Version6 is hash-verified** (SHA256 match `expected_source_hashes.json`),
  syntax-checked, synthetic assertion regression passes, embedded payload CRC
  valid (Fab_6ARU + EGF_1NQL references, configs, driver).
- **Limitation: Version6 covers only the corrected reference-folding controls**
  (4 arms: Fab/EGF × disulfide on/off, 1 sample each). It does **not** include
  the BoltzGen design-generation pipeline that produced the 86 designs —
  recovery of the generation code from notebook history is in progress.

## 3. Open questions

1. No pH-selectivity evidence for any candidate (the #1 criterion).
2. Official novelty pipeline code: promised Oct 5, not yet verified.
3. Submission needs a Proteinbase login (account/method unconfirmed); track
   unconfirmed (affects 40 vs 20 design cap).
4. Glycan clash on cC6 unresolved; cC5_00 clashes.
5. Mouse cross-reactivity cannot be confirmed computationally (but mouse is
   tested at pH 6.5, so a pH-switch design is compatible with it).
6. Design-generation code not yet recovered (notebook history in progress).

## 4. Running log (newest first)

### 2026-10-05 ~19:05 CEST — Structural pre-analysis of the 9 candidates
- Binder chains: 54–74 aa (all within the 10–250 aa challenge range);
  sequences in `data/candidates/binder_sequences.fasta`.
- Numbering resolved: model His406 = UniProt P00533 His433 (model numbering
  = UniProt − 27, verified via the LKEISDGDVI motif).
- Closest binder-acid (Asp/Glu CA) to a target His (CA-CA):
  rank03 5.7 A (Asp37–His406), rank06 6.7 A (Asp52–His406),
  rank07 8.2 A (Glu12–His406), rank01 8.5 A (Glu19–His331).
  rank03/rank06 already place a binder Asp near His433: prime starting
  points for the Proton-PottsMPNN pH redesign (Path B).

### 2026-10-05 ~18:45 CEST — All 9 candidate complexes recovered
- `data/candidates/`: rank01..rank08 CIFs (482–498 kB each) + cC5_00.cif.
  All verified as binder–target complexes (chains A+B, ~5.0–5.3k atoms).
  Ready for Proton-PottsMPNN Path A (placement scan / selectivity gaps).

### 2026-10-05 ~18:40 CEST — GitHub push fixed and verified
- Fixed `gh_push.py` technical debt: commits now preserve history
  (parent = previous tip, base_tree reuse) and accept `--message`.
  Also excludes `__pycache__`/`.pyc` (a 551 kB `.pyc` was causing HTTP 400).
- Pushed commit `83fe8e62` (parent `77eb506d`) with the Proton-PottsMPNN
  brief, DGX survey, candidate recovery and the pH-switch batch draft.

### 2026-10-05 ~18:30 CEST — First candidate structures recovered
- Downloaded via time-limited signed GCS URLs (minted from the Kaggle output
  file preview): `data/candidates/rank01_cC5_03.cif` (492 kB) and
  `data/candidates/cC5_00.cif` (444 kB). Both verified as genuine mmCIF
  binder–target complexes (chains A+B, ~5k atoms).
- Remaining 6 final-8 CIFs being fetched the same way.

### 2026-10-05 ~18:30 CEST — Kaggle batch archives located
- Found in notebook `notebook6198e9562d`: `C5_all_outputs.zip` (19.11 MB,
  version 1, scriptVersionId=354935690) containing
  `oc_C5/final_ranked_designs/final_8_designs/` with 8 ranked candidate CIFs
  (rank01_cC5_03 … rank08_cC5_07) plus `intermediate_designs/` with
  cC5_00.cif … cC5_14.cif; and `DIII_all_outputs.zip` (19.07 MB, version 2).
- **Not found**: no C2/C3/C4 archives in any of the 6 notebook versions, so
  the cC_6 / cC_12 / cC4_00 / cC4_15 complex backbones may be unrecoverable
  from Kaggle. If so, Proton-PottsMPNN Path B will start from C5 backbones
  (cC5_00 passed filters; the final-8 are batch-ranked).
- Browser zip-download clicks did not land on the VM; retrying with
  individual small CIF downloads.
- DGX driver-script downloads blocked: the portal's Files app triggers native
  "Save As" dialogs the browser automation cannot capture. Parked: not on the
  critical path. The RAM-workaround logic is documented from the survey; the
  originals can be re-fetched via a quick user takeover when the Boltz-2
  folding step needs them.

### 2026-10-05 ~18:20 CEST — DGX history surveyed (user tip)
- Portal survey (read-only): 4 batch-connect logs from Oct 4 evening
  (jobs 426–429: inputs upload, EGF runs, Fab runs). All 5 driver scripts
  intact, including the RAM-workaround variant `driver_dgx_phase.py`
  (RSS guard at 5.5 GiB, malloc_trim, phase markers) — downloading it plus
  the Fab provenance/log for reuse in the Boltz-2 folding step.
- Run outputs all present: egf_dgx_out, egf_phase_notrim, egf_phase_trim
  (with successful .npz/.cif predictions), fab_dgx_out, fab_phase_trim.
- Important: the DGX never ran candidate generation — no candidate binder
  structures there. Candidates live only in Kaggle notebook outputs.

### 2026-10-05 ~18:15 CEST — Proton-PottsMPNN brief done; plan unlocked
- Brief at `docs/proton_potts_mpnn.md` (code-grounded; paper PDF unfetchable).
  Key result: **design is CPU-only** (21 MB checkpoint), no GPU needed.
- Plan: **Path A** — rescore the 7 candidates with the placement scan
  (CPU-minutes per backbone); **Path B** — pH-redesign 2–3 backbones
  (prefer cC4_00/cC_12, avoid cC6 novelty risk) with λ sweep, fold top ~10 each
  in Boltz-2 on the DGX, re-apply filters.
- His433 mechanism confirmed exploitable: binder ASP-D/GLU-D vs target HIS-P =
  salt bridge at pH 6.5, lost at 7.4 (matches PMC9703009 precedent).
- Installing Proton-PottsMPNN locally now (uv + Python 3.12 + foundry);
  recovering candidate complex structures from Kaggle in parallel.

### 2026-10-05 ~18:15 CEST — Draft of next (pH-switch) batch prepared
- `code/generation/ph_switch_batch_DRAFT.py`: parameterized copy of the
  validated pipeline, targeting the H433 region. Patch residues and pH
  mechanism are explicit TODOs awaiting the Proton-PottsMPNN brief — not to be
  run as-is. Goal: launch within minutes of the brief arriving.

### 2026-10-05 ~18:10 CEST — Generation pipeline recovered
- Recovered the BoltzGen **design-generation code** from Kaggle notebook
  Version 2 (`scriptVersionId=355016993`, 2026-10-04): full pipeline in
  `code/generation/diii_batch_v2.py` + `code/generation/README.md`.
- Key parameters: 16 designs, binder 50–90 aa, 17-residue DIII patch
  (369,395,397,400,402,403,406,407,408,427,428,431,446,447,449,451,452 in 6ARU
  numbering), empty local MSA, filters rmsd ≤2.5 / iptm ≥0.6 / plddt ≥75 etc.
  Original run: 2h48m on 2×T4, 0/16 passed filters.
- This unblocks new design batches: retargeting = changing `LABELS` + output dir.
- Still in progress: Proton-PottsMPNN deep dive (paper + code).

### 2026-10-05 ~18:10 CEST — Report created; Slack assay details confirmed
- Created this report; repo `ph-selective-egfr-binder` live on GitHub with the
  workspace, STATUS.md, and the verified Version6 code.
- Slack (read-only) findings: mouse EGFR tested **only at pH 6.5**; exact
  constructs Human Met 1–Ser 645 / Mouse Met 1–Ser 647 (Sino Biological);
  Proton-PottsMPNN flagged by the organization for pH-sensitive design;
  novelty pipeline code still unverified.
- 30-minute status updates scheduled until the deadline.
- In progress: Proton-PottsMPNN deep dive (paper + code); recovery of the
  BoltzGen generation code from Kaggle notebook history.

### 2026-10-05 ~17:30–18:00 CEST — Access and recovery
- Signed in (with José's approvals): Kaggle (dataset downloaded, code
  validated), Comillas DGX (dashboard verified, 16 GiB / 8 h profile, nothing
  launched), Proteinbase Slack (read-only; CAPTCHA solved by José).
- GitHub repo created by José; connected via fine-grained token; SSH blocked by
  the sandbox proxy, pushed via the GitHub git-database API.
- Official challenge page re-read: deadline and ranking criteria confirmed.

### 2026-10-02 → 2026-10-04 — Prior work by "Instinct" (from transcript)
- 86 designs generated with BoltzGen across 9 batches; 7 passed internal
  filters; extensive cross-checks (Boltz-2, Chai-1, DockQ, disulfide controls,
  glycan and novelty analyses); pH-hypothesis literature research.
- Full detail in `docs/STATUS.md`.

### 2026-10-05 ~20:15 CEST — Path A complete: placement scan of all nine recovered candidates

**Status: DONE. All nine candidates scanned; `data/candidates/placement_scan.json` saved and committed.**

After a chain of environment fixes (below), the Proton-PottsMPNN placement
scan ran end to end on the recovered C5 backbones. For each candidate it
reports where the extended-vocab protonated states (HIS-P, ASP-P, GLU-P)
fit best on the binder, split by structural class (interface / core /
surface). Lower score = better fit. This is a **model prediction, not
experimental evidence**: it ranks which backbones are most amenable to
pH-switch redesign, it does not prove binding or pH selectivity.

Results, best interface HIS-P placement per candidate (lower is better):

| candidate | best iface HIS-P | residue | best iface ASP-P | residue |
|---|---|---|---|---|
| cC5_00 | -3.26 | 6 | -2.06 | 28 |
| rank01_cC5_03 | -2.56 | 13 | -0.47 | 16 |
| rank02_cC5_13 | -1.93 | 5 | -2.58 | 36 |
| rank05_cC5_06 | -1.79 | 38 | +4.62 | 44 |
| rank03_cC5_01 | -1.70 | 33 | -0.46 | 16 |
| rank04_cC5_02 | -1.43 | 7 | +2.26 | 24 |
| rank07_cC5_10 | -1.10 | 26 | +1.50 | 26 |
| rank08_cC5_07 | -0.95 | 7 | +0.33 | 48 |
| rank06_cC5_09 | -0.72 | 49 | -2.11 | 2 |

Reading (two independent pH-switch mechanisms):
- **Receptor-His mechanism** (from the earlier geometry): a binder acidic
  residue facing EGFR His433 forms a salt bridge only when His433 is
  protonated (pH 6.5), losing it at pH 7.4. rank03 leads here: Asp37 sits
  5.7 A (C-alpha) from His433, followed by rank06 (Asp52, 6.7 A).
- **Binder-His mechanism** (this scan): a binder histidine at the interface
  that prefers protonation could pair with a receptor acid at pH 6.5.
  cC5_00 leads here with the strongest interface HIS-P signal (res 6: -3.26,
  res 28: -2.16), but it carries a strong glycan-clash warning from two
  alignments and needs that filter before any redesign.
- **rank03 combines both**: best H433 geometry (Asp37, 5.7 A) plus a
  moderate HIS-P signal (-1.70 at res 33) and a mild GLU-P fit (-1.16 at
  res 16). It is the leading Path B backbone.
- rank06, second-best on H433 geometry, has the weakest binder HIS-P
  signal (-0.72); the model does not favor protonated His at its interface.

Path B plan: redesign the top 2-3 backbones with
`inference/design_ph.py` (block descent, `combined_lambda` sweep across
0.1 / 0.3 / 0.5 / 0.7, multiple seeds), starting with rank03 and cC5_00
(pending its glycan recheck), then fold selected sequences with Boltz-2.

### 2026-10-05 ~19:00–20:10 CEST — Proton-PottsMPNN install and Path A debugging

- Compiled HBPLUS from the official EBI/GitHub source
  (`~/workspace/ppmnn/hbplus/hbplus/hbplus`); the v6 featurization calls it
  via `HBPLUS_PATH`. The earlier brief's claim that design needs no external
  HBPLUS was wrong for this code snapshot; corrected in
  `docs/proton_potts_mpnn.md`.
- Reinstalled the venv CPU-only (torch 2.14.1+cpu, CUDA stack removed after
  José asked about acceleration).
- Installed undeclared runtime deps one by one as the scan's import chain
  surfaced them: `ipdb`, `propka`, `flaml`, `scikit-learn`, `xgboost`,
  `lightgbm`. Note: an early quiet `pip install xgboost lightgbm` failed
  silently; verified with explicit output on retry.
- `code/placement_scan_A.py`: Biotite API fix (PDBxFile -> CIFFile), Atomworks
  parsing for charge annotations, finite-coordinate SASA split, HBPLUS_PATH
  export. All nine candidates now complete without errors.

### 2026-10-05 ~18:45 CEST — Candidate structures recovered

- Recovered 9 structures from Kaggle notebook `notebook6198e9562d`
  (C5 batch): `data/candidates/` rank01-08 plus cC5_00. Binder 54-74 aa.
- C2/C3/C4 archives not found in any notebook version; cC_6, cC_12,
  cC4_00, cC4_15 structures remain unrecovered.
- Binder sequences in `data/candidates/binder_sequences.fasta`.
- Model numbering offset +27 vs UniProt (model His406 = EGFR His433,
  verified via LKEISDGDVI motif).
