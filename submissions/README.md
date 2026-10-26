# Submissions

## challenge3_glp1r_designs_PREPARED.csv (26 Oct 2026, 12:09 CET)

Ten de novo designs PREPARED for the GLP1R challenge 3 form (columns as the form
requests: name, sequence, molecule_class; single_chain per the form's own wording).

Status: PREPARED, upload NOT confirmed from this side. The user submitted through
the official form with his own login before the deadline (26 Oct 12:59 CET); which
of these designs (or how many) were actually uploaded has not been confirmed back.
Until it is, this file records what was prepared, not what was submitted.

Evidence strength per design (all computational prediction, no experimental
validation; see challenge3/results/zeroshot_409_v2/ and boltz_kaggle_409/):

- gA_3, gB_2, gB_1: pass the two-seed Boltz filter on the official 409 aa target
  (metrics only; poses not reproducible between seeds; native control scores
  equal or better).
- gA_2, gA_1: rejected by that same filter (wandering site / low ipTM).
- g3_1: fails the local steric gate (see challenge3/results/README.md, traj1).
- g3b_1: weak / rejected.
- v6_top02: scores below the native control in the 409 aa filter.
- v6_top05, v6_top14: ESMFold consistency check only, no independent filter.

## Challenge 1 and 2

Nothing was ever submitted for challenges 1 or 2 (see challenge2/results/
submission_format_DRAFT_NOT_SUBMITTED.csv).
