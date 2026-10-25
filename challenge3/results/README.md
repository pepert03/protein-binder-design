# Challenge 3 results (GLP1R agonist)

## traj1 (first DGX trajectory, 24-25 Oct 2026) - FAILS the steric gate

- `traj1/g3_1.pdb` - best-model export of trajectory 1 (ColabDesign binder, 36 aa
  sequence PPSPKIWEAHEAQARYDYLRIAIEDHVSVGMFGGAM; two MODELs from the ensemble of the same best step - one candidate).
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
helix-end boundaries) runs without OOM (but see the VRAM correction under traj2 - that 12,798 MiB figure was a whole-card nvidia-smi reading, not the slice), the
pipeline exports best-consistent artifacts, and the resseq->UniProt mapping
used by the gate was verified fragment-by-fragment against the real P43220
sequence (6/6 exact).

Earlier contact analyses of this pose used a circular mapping and are
superseded; they are not published here.

## Notation: two MODELs = one candidate

The PDB exports carry two MODEL records that come from the ensemble of the
SAME best optimization step (confirmed against the colabdesign source in the
DGX venv). They are not two independent candidates. The run metrics in the
log are averages over the ensemble members (e.g. plddt 65.7 / 64.6 -> mean
65.1, matching the log). The steric gate is evaluated per MODEL; both models
of traj1 fail it, so the verdict is unchanged.

## DGX-side tooling

`scripts/supervisor_dgx.sh` - the actual supervisor that ran on the DGX for
run g3b: it waits for the composite condition (TRAJ traj=1 line in the log +
g3b_results.json containing "traj": 1 + json mtime stable for 20 s +
non-empty g3b_1.pdb), copies the artifacts to export/ BEFORE stopping the
process, and leaves a g3b_SUPERVISOR.txt record. DGX-only; not for the
Umbrel.

## traj2 (run g3b traj 1, 25 Oct 2026) - FAILS the steric gate; NOT validated

- `traj2/g3b_1.pdb` - best-model export of trajectory 1 of run g3b (driver
  v4), binder 36 aa TEYAKKQMYDYMLNDLEMLTGCVPAGFFFGGFFGGG; two MODELs from
  the ensemble of the same best step (one candidate, see the Notation
  section above).
- `traj2/g3b_1_gate.txt` - per-MODEL pose gate, reproduced identically on
  the Umbrel clone. MODEL 1: 15 pairs < 2.0 A, min 0.71 A; MODEL 2: 54
  pairs, min 0.49 A. Both models clash severely in the same receptor region
  (around U367 PHE / U369 PHE / U368 ALA).
- `traj2/g3b_results.json` - run metrics (plddt 0.445, i_pae 0.543, i_ptm
  0.461). Sequence verified identical to PDB chain B and to the TRAJ line
  in the log.
- `traj2/g3b_run_sanitized.log` - run log with personal paths replaced. It
  also shows a second trajectory had just started (one line) when the
  supervisor stopped the process after its composite export check.

Verdict: NOT a valid candidate. Besides the clashes, the ensemble is
unstable: the two models of the same best step disagree sharply (15 vs 54
clash pairs), so the pose has not converged. Do not read traj2 as an
improvement over traj1 on any single axis: fewer clash pairs in MODEL 1
(15 vs 27) does not imply a globally better pose - the models disagree
with each other and both fail the gate. No binding or activation is
claimed.

VRAM correction (applies to traj1 notes too): the nvidia-smi figures
reported so far (12,798 MiB for traj1; 33,746 MiB in this log, of
143,771 MiB total) are whole-card readings, not the MIG slice and not
this process. Slice/process memory use is currently unmeasured; "runs
without OOM" is the only supported claim. Two consecutive trajectories
without OOM do not prove the clear_mem fix is stable.

## Orientation analysis (25 Oct 2026)

`../docs/orientacion_gate.md` (observed orientation, reproduced with
`../scripts/orient.py`; per-run outputs in `traj1/orient_g3_1.txt` and
`traj2/orient_g3b_1.txt`): in 2/2 trajectories the geometry is INVERSE to the
native class-B1 binding mode (N-terminus sits on ECD/ECL, C-terminus in the
TM core). This is an observed geometry plus an UNPROVEN causal hypothesis,
not evidence about activation. `../docs/propuesta_end_conditioning.md` lists
candidate design changes (end-conditioning); none has been executed and no
agonism is claimed or proven for anything in this repository.

## Local redesign viability (25 Oct 2026)

`../docs/viabilidad_rediseno.md`: gates ensemble_rmsd.py and orient_gate.py
implemented and reproduced (0.49/1.13 A and 1.62/4.93 A; G4 = 27/24 and
15/54), 6X18 corridor calibration (native frame = calibration only under the
official de novo/zero-shot rule), target-only sampling. No zero-clash pose
demonstrated; recoverability claims withdrawn; nothing here proves redesign
viability, sidechains, binding or agonism. Run outputs in `viability/`.
