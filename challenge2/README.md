> Historical planning/results text below is superseded by [LATEST_STATUS.md](LATEST_STATUS.md), updated11October2026. No supported pH-switch binder or submission.

# Challenge 2 - TNF-alpha (Anthropic x Adaptyv 2026, Proteinbase)

Goal: design a de novo binder to the soluble human TNF-alpha trimer
(residues 77-233, UniProt P01375-1, PDB 1TNF) that binds at pH 7.4 and shows
no detectable binding at pH 6.0. This is the reverse pH direction of Challenge 1
(EGFR: bind 6.5, release 7.4).

Ranking of objectives: (1) pH selectivity (bind 7.4, not at 6.0),
(2) mouse TNF-alpha cross-reactivity at 7.4 (about 79% identity, 124/156),
(3) affinity.

Official page: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha

**Deadline: Sunday 2026-10-11 23:59 AoE (= Monday 2026-10-12 13:59 CEST).**

Submission: up to 20 de novo designs (10-250 aa; protein, nanobody, scFv or Fab),
CSV ranked by preference (name, sequence, molecule_class). The method description
becomes public on submission.

## Status

Research done, nothing designed yet. Nothing submitted to Proteinbase and no money
spent. Submission needs Jose's explicit go-ahead.

## Layout

- `docs/RESEARCH.md` - problem summary and methods survey.
- `code/` - design and analysis scripts.
- `data/` - target structures (1TNF and asymmetric ligand-bound structures), configs.
- `submissions/` - final CSVs. Empty.

## Plan sketch (hypotheses, not results)

1. Pick an epitope between protomers, conserved in mouse.
2. Generate with BindCraft or BoltzGen / ColabDesign hallucination.
3. Place histidines at the interface so binding weakens when protonated at pH 6.0.
4. Re-predict every design independently; filter by confidence and novelty.
5. Rank 20 by pH-selectivity rationale first.

Compute caveat: a trimer is a large complex; a 16 GB Kaggle T4 may be tight for BindCraft.
