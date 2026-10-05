# ProteinBase EGFR Challenge

Design of a conditional EGFR binder for the ProteinBase challenge
(anthropic-adaptyv-2026, Challenge 01).

Ranking of objectives: (1) pH-selective binding (binds at pH 6.5, no detectable
binding at pH 7.4), (2) mouse cross-reactivity, (3) affinity to human EGFR.

**Deadline: 2026-10-06 23:59 AoE (= 2026-10-07 13:59 CEST).**

## Where things stand

Read [`docs/STATUS.md`](docs/STATUS.md) first: it reconstructs all prior work
(Oct 2-4, 2026) from the conversation transcript: batches run, the 7 candidates
that passed internal filters, controls, compute setup, Slack findings, and open gaps.

Short version: 86 designs generated with BoltzGen, 7 passed internal filters,
none has any evidence of pH selectivity (the #1 criterion). Nothing has been
submitted.

## Layout

- `code/` - design and analysis scripts (recovered from Kaggle).
- `data/` - inputs: receptor structures, reference complexes, configs.
- `results/` - batch outputs: metrics, structures, diagnostics.
- `notebooks/` - exploration notebooks.
- `submissions/` - final submission CSVs (name, sequence, molecule_class).
  Empty: nothing submitted yet.
- `docs/` - STATUS.md and notes.

## Rules

- Nothing is submitted to Proteinbase and no money is spent without Jose's
  explicit go-ahead.
- The methodology description is made public on submission: review it before
  sending anything.
- Results that live only on the Comillas DGX or Kaggle must be copied here;
  those machines are not backed up.

## Quickstart (resume)

1. Code recovered from the Kaggle dataset `proteinbase-codigo-version6` into
   `code/kaggle_v6/` (Version6, hash-verified, syntax + regression + payload CRC
   checked). This version covers the corrected reference-folding controls only,
   not the design-generation pipeline.
2. Access confirmed (2026-10-05): Kaggle (free GPU quota available), Comillas
   DGX (student profile: 1x16 GiB MIG slice, 8 h max per session), Proteinbase
   Slack (read-only).
3. Check the official novelty pipeline code (promised in Slack for Oct 5).
4. Pick the final candidate set (see open gaps in `docs/STATUS.md`) and prepare
   the submission CSV in `submissions/`.
