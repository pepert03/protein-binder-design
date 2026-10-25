# Challenge 3 results (GLP1R agonist)

## traj1 (first DGX trajectory, 24-25 Oct 2026) - FAILS the steric gate

- `traj1/g3_1.pdb` - best-model export of trajectory 1 (ColabDesign binder, 36 aa
  sequence PPSPKIWEAHEAQARYDYLRIAIEDHVSVGMFGGAM, 2 MODELs).
- `traj1/g3_1_gate.txt` - per-MODEL pose gate (`scripts/gate_pose.py`, reproduced
  independently on the Umbrel clone): MODEL 1: 27 atom pairs < 2.0 A, minimum
  0.56 A; MODEL 2: 24 pairs, minimum 0.93 A. Both models show severe steric
  clashes (mainly binder Gly33-35 backbone vs Phe369 / Leu314 of the receptor).
- `traj1/g3_run_sanitized.log` - run log with personal paths replaced; ends in
  the clear_mem defect (RuntimeError: Array has been deleted) fixed in driver v4.

Verdict: NOT a valid candidate. The steric clashes mean this pose cannot
represent a physically plausible bound state. No binding or activation is
claimed: the pose gate is an in silico sanity check, and AF2/ColabDesign
metrics (i_pae, plddt) proxy binding confidence only - they say nothing about
GLP1R activation, which is the actual challenge objective.

What traj1 does establish (engineering, not science): the v3 crop (316 res,
helix-end boundaries) runs on the 16 GiB MIG without OOM at 12,798 MiB, the
pipeline exports best-consistent artifacts, and the resseq->UniProt mapping
used by the gate was verified fragment-by-fragment against the real P43220
sequence (6/6 exact).

Earlier contact analyses of this pose used a circular mapping and are
superseded; they are not published here.
