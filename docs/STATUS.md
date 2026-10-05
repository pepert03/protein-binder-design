# ProteinBase EGFR Challenge, Status

Source: WhatsApp export of Jose's conversation with a previous AI assistant ("Instinct"),
Oct 2-4 2026, file `~/workspace/user/files/2026-10-05_0_l0kg.md` (605 lines).
Work resumed by Muse on Oct 5 2026. This file is a faithful reconstruction; anything not
in the transcript is marked [uncertain].

## 1. Challenge summary

- Competition: ProteinBase, anthropic-adaptyv-2026, Challenge 01: EGFR.
  Official page: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
- Task: design a conditional binder against human EGFR extracellular region.
  Domain III recommended (functional epitope, targeted by cetuximab/panitumumab).
- Ranking of objectives (this order matters):
  1. pH-selective binding: binds at pH 6.5, no detectable binding at pH 7.4.
  2. Mouse cross-reactivity: same sequence binds human and mouse EGFR.
  3. Affinity against human EGFR.
- Deadline: Tuesday 2026-10-06, 23:59 AoE (= 2026-10-07 11:59 UTC = 13:59 CEST Madrid).
  Extended from Oct 4; announced in Slack Oct 2, later confirmed on the official page.
- Submission rules (from official FAQ):
  - De novo designs only, single chain / nanobody / scFv / Fab, 10-250 aa.
  - CSV ordered by the submitter's ranking, columns: name, sequence, molecule_class
    (protein | nanobody | scfv | fab_kappa | fab_lambda).
  - Up to 40 designs for Track 1, up to 20 for Tracks 2 and 3.
    [uncertain: which track Jose is in was never pinned down in the transcript]
  - Encouraged: metrics (ipTM, ipSAE, self-consistency, physics-based, liability
    scores), structure models, methodology write-up (methodology will be made public).
  - Only one submission batch counts: the last one sent or the one explicitly designated.
  - Selection for synthesis via a Claude-based workflow; all data incl. failures
    published on Proteinbase afterwards.

## 2. Timeline of work so far

- Oct 2 13:45 - Oct 4 21:50: "Instinct" (another AI) did the research and compute,
  with Jose approving steps. No money spent, nothing submitted.
- Oct 5: Muse resumes the work. Repo scaffolded at `~/workspace/proteinbase-egfr/`.

## 3. Design pipeline used

1. BoltzGen: de novo binder generation against a chosen receptor patch, then sequence design.
2. Boltz 2: refolding of the designed complex (5 samples in most batches), confidence
   metrics (ipTM, PAE, pLDDT).
3. Chai-1: orthogonal structural check (poses reproduced for 5 of 7 candidates;
   note: Chai-1 runs used the truncated relevant EGFR domain, Boltz 2 used the full
   receptor, so the comparison is not equivalent).
4. DockQ (official implementation): interface comparison against reference complexes.
5. Internal filters (all must pass):
   - Complex RMSD (refolded vs generated complex) <= 2.5 A.
   - Binder contacts the requested target patch.
   - No clashes: no interface contacts below 2.5 A.
   - Interface confidence (ipTM) judged case by case, not a fixed cutoff.
   - Novelty: partial checks only, local sequence similarity plus structural
     comparison (TM-align style) against PDB entries through July 2025.
     NOT equivalent to the official novelty pipeline.

## 4. Batches run

| Date (Oct) | Batch | Designs | Passed | Key outcome |
|---|---|---|---|---|
| 2, ~15:36 | Test batch | 2 | 0 | Low interface confidence; kept as exploratory only |
| 2, 16:12 | Numbering fixed | 4 | 0 | None passed filters |
| 2, 17:39 | Patch 2 | 4 | 0 | Interface confidence slightly higher; top 2 missed target patch |
| 2, ~19:56 | Relaunch (session lost) | 4 | 0 | 1 touched patch: ipTM 0.665, min PAE 5.6 A, but one too-close contact |
| 2, 21:53 | Evening batch | 8 | 2 | Best ipTM 0.81, better refold self-consistency; clashes not yet checked |
| 3, 04:06 | 16-design batch (from 01:05) | 16 | 2 | 1 convincing: ipTM 0.713, patch contacted, no contacts < 2.5 A |
| 3, 07:44 | Focused 16-design batch | 16 | 3 | 1 with better interface signals, 2 weaker (one design initially miscounted, corrected to 3/16) |
| 3, 20:12 | 2 fixed histidines | 16 | 0 | All failed complex RMSD; His were conserved in sequence (initial misread corrected) but not constrained to the interface |
| 4, 04:04 | New DIII patch | 16 | 0 | 15/16 passed binding-site filter; all failed complex RMSD; 80 refold samples (5/design) checked, none passed |

Total: 86 designs generated, 7 passed internal filters (~8%).
Running judgment at the time: 7 pass, but only 4 have the best interface signals.
## 5. The 7 candidates that passed internal filters

IDs as they appear in the transcript (naming is inconsistent: cC_6/cC6, cC_12/cC12).
Two of the seven are never named [uncertain].

- **cC_6** (batch C2): most pose-consistent across Chai-1 and Boltz-2 contrasts.
  Places a Glu near receptor H433 (pH hypothesis, weakened by controls).
  Issues: HIGH structural-novelty risk (local comparison found high similarity to
  part of a PDB protein); glycan clash risk is conformation-dependent and unresolved.
- **cC_12** (batch C3): Glu near H433 in Chai-1 models. Issues: pose consistent
  in only 1 of 3 Boltz-2 models.
- **cC4_00**: more favorable structural-novelty signal than cC6. Issues: position
  varies strongly between Chai-1 and Boltz-2.
- **cC4_15**: avoids the N352 sugar in the tested context. Issues: problems when
  full-length EGFR domains are added.
- **cC5_00**: issues: strong clash with the first EGFR-attached sugar in two
  different alignments.
- **2 unnamed designs**: passed filters but with weaker interface signals
  [uncertain: IDs and batch of origin not in transcript].

General caveats on all candidates: Boltz-2 scores high even for shuffled-sequence
controls (its score is not discriminative here); Chai-1 separates designs from
shuffled controls in both score and position, but is not an independent validation.
No candidate has any affinity, pH-selectivity, or mouse-binding evidence.

## 6. Controls

- **EGF (PDB 1NQL)**: 5 samples (1 per disulfide mode). Complex RMSD 11.98-24.06 A
  (threshold 2.5), DockQ 0.048 / 0.014. Later finding: the input prep had dropped
  explicit disulfide connectivity (3 disulfides), so this tested a simplified input,
  not a chemically equivalent reference. Additionally, 1NQL is an autoinhibited,
  low-affinity state crystallized at very acidic pH, not the active DIII complex,
  so it is not a general negative control for the pipeline.
- **EGF replicas (Kaggle, 5 seeds x 2 disulfide modes, 1 sample each)**: high
  seed-to-seed variability. 2 of 5 seeds reach DockQ 0.43-0.46 in both modes.
  48 resolved residues sit ~5.4 A from the reference fold in both samples;
  mean pLDDT 0.50 (EGF) vs 0.74 (receptor). DGX vs Kaggle gave a practically
  identical pose by DockQ (timing 2:56 vs 6:53 for one prediction).
  Recording seeds does not control all internal randomness: disulfide-mode pairs
  are not perfectly paired comparisons.
- **Fab (PDB 6ARU, cetuximab Fab mutant)**: 1 sample per disulfide mode.
  Recovers 88-92% of native interface contacts, interface RMSD ~0.61 A,
  DockQ 0.72-0.84 (both receptor interfaces, both modes). Full Fab + 5 simultaneous
  samples OOMed on Kaggle; 1 sample fits. DGX run (after RAM workaround): DockQ
  0.78 / 0.82 for the two interfaces.

## 7. pH-selectivity hypothesis research

- Literature lead: an engineered (not de novo) anti-EGFR antibody placed Asp/Glu
  facing receptor H433 (domain III) and gained binding in acidic medium
  (PMC9703009). Its final construct still bound at pH 7.4. Lead, not a solution.
- Second lead: binder histidine facing a receptor Glu for acidic binding
  (10.1080/19420862.2019.1682866).
- Geometric scan: in Chai-1 models, cC_6 and cC_12 place a Glu near receptor H433
  (checked against structure 1IVO, the active state with the two receptor
  histidines). Weakened: shuffled-sequence controls also contact H433, and the
  contact positions differ between Chai-1 and Boltz-2.
- A rigid-geometry His-placement scan completed but used an unvalidated approximate
  histidine geometry; its distances are not evidence of feasibility.
- Protonation-state joint-sampling audit: checking sites independently can hide
  differences between joint states (tested on synthetic data; useful for auditing
  any future pH calculation, not a molecular simulation).
- Proton PottsMPNN (shared in Slack for protonation-aware design) is a PD-L1
  example, not validated for this EGFR setup.

## 8. Compute setup

### Kaggle (account: jridruejotunon@gmail.com, private notebooks)
- Free GPU quota; notebooks used: notebook6198e9562d (main, "Version6"),
  notebookf13ce1b745 (EGF replicas).
- Measured: 2x Tesla T4, 15 GB each; design jobs use ~9-10 GB per GPU.
  Roughly 50 min per 4-design batch (gen + seq design ~19 min, folding + eval ~29 min).
- Sessions can shut down and lose unsaved results (happened once; fixed with an
  auto-collect script afterwards).

### Comillas DGX (https://hpc.comillas.edu, Open OnDemand, no VPN, SLURM)
- Default student profile: 1 GPU of 16 GiB (a MIG 1g.18gb slice of an H200),
  sessions up to 8 h (announcement email said 24 h; portal shows 8 h).
- Measured inside a session: 4 CPU cores, ~5.86 GiB RAM per job (binding limit),
  Python 3.10.12, PyTorch 2.14.0, CUDA 12.6. Internet OK (PyPI, GitHub, HF).
  GPU memory was NOT the limit; RAM was (Fab 1-sample OOMed on RAM first, then
  succeeded after releasing temporary memory after model load, with no change to
  weights or scientific parameters).
- EGF prediction: 2 min 56 s on DGX vs 6 min 53 s on Kaggle (single comparison).
- Storage is NOT backed up (RAID 0): always copy results off. Working dir:
  clusters/dgx in the home folder. Delete sessions via "My Interactive Sessions".

## 9. Slack findings (Proteinbase workspace, read by Instinct Oct 2-4)

- Deadline extension to Oct 6 23:59 AoE announced Oct 2 (page confirmed later).
- Novelty: Level 3 required, checked at submission; sequence comparison alone is
  not enough, structural comparison included; broader official pipeline, code
  promised for publication on Oct 5.
- Affibodies: novelty filter has no antibody-specific treatment yet for EGFR
  (planned for a later challenge, not this one).
- Assay (Tudor, Adaptyv): binding measured mainly by SPR, BLI to validate hits;
  our designs are immobilized and EGFR is passed over them. Human and mouse EGFR
  produced in HEK293, i.e. glycosylated: the glycan alerts are relevant to the
  real assay, not just modeling.
- Team with Joaquin: put team name in profile, both names in bio and methodology.
- Can target another domain/epitope; no need to fill 20 designs; method novelty,
  diversity and quality valued over quantity.
- Unanswered or unclear at last check: exact human/mouse EGFR fragments used,
  assay conditions, pH-selectivity threshold; novelty pipeline code (due Oct 5).

## 10. Open gaps (as of Oct 4 21:50)

1. No evidence of pH selectivity for any candidate (the #1 ranking criterion).
2. Official novelty pipeline code (promised Oct 5) not yet reviewed.
3. No submission made; submission form needs a Proteinbase login
   (account/method never confirmed).
4. Glycan clash risk on cC6 unresolved (conformation-dependent); cC5_00 clashes.
5. Mouse cross-reactivity: candidates touch residues that differ human vs mouse;
   cannot be confirmed or ruled out computationally.
6. Track not confirmed (affects max designs and selection mechanics).
7. Full Fab 5-sample control never run (RAM); only 1-sample done.
8. Code/results recovery: Kaggle dataset "proteinbase-codigo-version6" (Version6
   code + execution examples); the tmpfiles ZIP is expired.

## 11. Notes on transcript reliability

- Candidate naming is inconsistent (cC_6 vs cC6); two passing designs unnamed.
- "C2"/"C3" batch labels are used but never fully defined.
- The disulfide option in the BoltzGen version used contradicted its own docs;
  it was left unchanged and noted, not silently flipped.
- Several interpretations were corrected in-conversation (histidine conservation,
  DIII CIF/sample mismatch, EGF disulfide prep): the corrected versions are used above.
