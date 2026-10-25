# Challenge 3 strategy (GLP1R agonist) - honest version, 24 Oct 2026

## What makes this challenge different

Challenges 1-2 asked for binders. This one asks for an **agonist**: binding is necessary
but not sufficient, and the winner is ranked by *activation* relative to a reference
agonist. No metric in our current pipeline (ipTM / ipAE / pLDDT) predicts activation;
they only proxy binding. Any claim that a design "should activate" is a hypothesis, and
the method write-up must say exactly that.

## What the reference structure shows (computed, not recalled)

From 6X18 (GLP-1 bound to GLP1R, cryo-EM 2.1 A), `scripts/prepare_target.py` output:
- 52 receptor residues contact the peptide (heavy atoms <= 4.5 A): **22 in the ECD,
  30 in the TM/ECL core**. Contacts span UniProt 30-391.
- The peptide touches the receptor along its whole length: N-terminal residues
  (His7, Glu9, Thr11, Phe12, Thr13, Ser14, Asp15, Val16, Ser17, Tyr19 in GLP-1
  numbering) sit in the TM bundle; C-terminal residues (Glu21, Lys26, Phe28, Ile29,
  Trp31, Leu32, Val33, Arg36) occupy the ECD groove.
- This matches the established class-B GPCR two-site picture (peptide C-terminus docks
  on the ECD, N-terminus engages the TM core where activation happens) - stated here as
  literature context, with 6X18 itself as the concrete reference in this repo.
- The ECD alone is ~120 residues. Modeling binder+ECD is a proxy for the binding
  objective only; it says nothing about activation and ignores the TM core. Treat any
  ECD-only result as incomplete evidence, per our own review rule.

## Design routes, with honest odds

1. **De novo micropeptide agonist (10-40 aa, "microbinder" class)** that engages both
   the ECD groove and the TM pocket with a novel sequence. Mechanistically the most
   direct route to agonism, and GLP-1 itself proves a 30-aa peptide can do it. Hard
   parts: (a) must be sequence-diverse from GLP-1/exendin-4 (de novo zero-shot rule -
   our `validate_submission.py` flags exact >=8-aa windows); (b) getting the N-terminal
   activation chemistry right without copying the natural sequence; (c) peptides in
   detergent-solubilised receptor assays can be conformationally unstable.
2. **De novo minibinder (40-100 aa) against the ECD groove.** Binding objective is
   plausible with our existing pipeline style. Activation via the ECD alone is not
   established pharmacology for GLP1R; a minibinder would most likely win the affinity
   mention but not the challenge. Still legitimate: ranking counts binding too, and
   ~375 Track 2+3 designs get tested, so a strong binder has real expected value.
3. **Nanobody vs ECD.** Binding plausible; agonism via an ECD nanobody is a long shot.
   Only worth it if we want class diversity in the 20 slots.

Recommended mix (if generation capacity appears): mostly route 1 micropeptides with a
few route 2 minibinders, ranked honestly. One strong submission of a few designs beats
20 weak ones: the selection prompt rewards method novelty and diversity, not volume.

## Compute reality (verified 24 Oct)

Umbrel: 4 cores, 15 GB RAM, no NVIDIA GPU. BoltzGen / Boltz-2 class inference does not
run here. The competition FAQ allows any tools and resources we can access, but our
active limits are: **no spending, no credits, no sequences to external services, no
submission** without Pepe's explicit go. So right now this pack supports preparation,
review, and method writing; actual generation needs a compute decision from Pepe
(Kaggle weekly quota, Comillas DGX, or accept a non-generated entry is impossible and
aim at challenge 4/5 with proper time). Nothing here promises otherwise.

## What we CAN do on the Umbrel today (CPU-only)

- Target prep (done here, re-runnable): `scripts/prepare_target.py`.
- Submission self-checks: `scripts/validate_submission.py`.
- Literature/structure analysis and the method write-up skeleton.
- Sequence-level design tools that are CPU-feasible (e.g. small ProteinMPNN runs on a
  fixed scaffold) - only if a scaffold exists, which today it does not.

## Review checklist before anything is ever submitted

- [ ] Pepe explicitly approves submission (channel evidence), CSV and method text reviewed by him
- [ ] `validate_submission.py` passes with zero FAILs and every WARN reviewed by a human
- [ ] Every design sequence generated de novo; provenance per design documented
- [ ] Method write-up states: binding metrics only, activation is a hypothesis
- [ ] Paratope not at the C-terminus (Twin-Strep capture tag in the assay)
- [ ] No design targets FLAG/His tags; no tag fusions in designs
- [ ] Team name on profile; both names in bio/methodology (standing competition rule)
