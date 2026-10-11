# Whole Challenge2 review and recommended plan
11October2026,19:37CEST. Source snapshot: main e73e9cf5a20914898e0cc55a5a1fbb27bbcc923a. Reviewed79 exported text/source files,19 Python files (all parse), git log, current live entry docs and locally preserved structural artifacts. Snapshot excludes PDB/images/weights; actual local PDBs fill that evidence gap for recent branches. Earlier lost outputs cannot be reconstructed faithfully.

## Answer

Our biggest problem is not too few histidines. The scaffold/interface evidence is weak, crop geometry can mislead, and historical code can score one prediction while inspecting another structure. Fix export/provenance first, then test whether a no-His scaffold can fold and form a stable full-trimer interface. Only after that should we engineer pH response. Continuing small His substitutions or a stagnant hard stage is not the recommended use of compute.

## Confirmed defects and fixes

1. **Historical stale-PDB export path.** ColabDesign1.1.1 save_pdb defaults get_best=True. run3 optimizes designs, then calls predict for MPNN variants using the same object, then saves tmp.pdb without disabling saved-best selection. Its current log and exported coordinates can diverge. hisCat on that file is therefore not trustworthy without recovery/revalidation. Patch explicitly exports current prediction with get_best=False. run5 and repredict are made explicit too as a guard. This identifies a code path, not proof every historical value is wrong. Fresh validators start clean and their all-model coordinate audit remains useful.
2. **Best sequence vs last log.** run3/run4 pair get_seqs default saved-best with aux last-step log. Patch selects the matching saved-best aux log. Future bundles must also match PDB sequence and exact model-specific metrics.366 hard logs were all identical at print precision, so no demonstrated large366 discrepancy.
3. **Multi-MODEL flattening.** Historical hiscat manually collects every ATOM line while ignoring MODEL boundaries. With two-model output, atoms from different models can be pooled and a minimum distance taken across incompatible structures. Require one contact table perMODEL, no pooled proximity. The supplied explicit-chain geometry auditor does this for4/5/8Å contacts; a protonation-aware His analysis remains to implement after a valid interface exists.
4. **Nonstandard JSON.** Archived compare/predictions.json contains NaN for absent alone iPAE. Patch emits null for missing metrics and allow_nan=False. Keep historical files unchanged; new data follows a strict schema.
5. **Stale entry docs.** README says nothing designed, method lists early.498 scores, results ends at seeded hit and mentions old running jobs. Add LATEST_STATUS.md and prominent supersession notices instead of rewriting history.
6. **Fragile setup/checkpoints.** Historical shell helper returns failure but callers continue; pip/download versions and hashes are not consistently recorded. Future runner must fail fast, pin versions, record parameter hashes/model names, verify GPU memory before first gradient, and checkpoint sequence/PDB/log after each stage. Do not rely on session-local files or long shell-injected base64 arguments.

## Confirmed scientific weaknesses

- Original20 have low diversity and17 share VVVVLLCCLLLLLLLILLL. This is a heuristic liability, not measured aggregation. Do not score it away with AF confidence.
- Regenerated110crop poses of control/repairs show466-481 sub2Å atom pairs against omitted full-trimer atoms after target fitting. Fit RMSD~2.2Å limits fine interpretation, but extensive overlaps make crop-only evidence unreliable.
- Fresh all-model re-audit: control, repair02,277 and366 each have two intact-human-trimer models, all zero<4Å binder contacts. No severe clash does not mean binding.
- Repair02/03 binder-alone vs complex CA RMSD~4.36/5.31Å, control4.12Å. Repair reduced hydrophobic run but produced low-complexity E/R sequences and did not improve shape consistency.
- Crops contain internal gaps.366 retains full152aa chainA but truncated B/C. Source prep preserves original residue-index gaps; no confirmed false immediate-neighbour joining. Still target completeness and artificial exposed surfaces remain risks. Assert index mapping and never draw a CA line across a gap in the preview.
- PAE is normalized by31 and hotspot-masked values differ from whole-target values. Old<.4 is not a validated success gate. Model-average logs also hide per-model variation.
- Mouse screen is threading; even its documented alignment leaves a human residue at a mouse gap. It is not an exact native mouse construct test.
- Histidine proximity at8Å or ring-N proximity at5Å does not establish a hydrogen bond, protonation state, pKa shift, or the required6.0 release. The cited TNF preprint result was at5.4, not6.0. A weaker effect at6.0 is plausible but not guaranteed. Known literature binders may calibrate methods, not serve as competition starting scaffolds: official rules disallow modifying an existing binder.

## One recommended plan: rescue the pipeline before the molecule

###A. Now, noGPU needed
Apply proposed export-state patches and current-doc notices; keep historical numbers as historical. Require sequence SHA256, target SHA256, modelID, seed, recycle count, crop manifest, metric masks and normalized/Å values for every future candidate. Run the new perMODEL auditor on every PDB and verify exact sequence before ranking. Revalidate old leading hypotheses rather than delete them based on another control.

###B. One controlled pilot, only after fresh account/compute approval
Use one documented geometry-preserving target representation and identical seed/binder length for baseline and one folding-focused objective. No His bias. Retain native trimer placement and omitted-atom audit. If larger-memory compute is available, multimer ensemble design is the preferred longer-term branch; do not pretend full456 gradients fit16GB because366 did. If constrained to the current slice, cropped design is a pilot with explicit limitations, not a substitute for intact-target validation.

Keep the pilot small: baseline plus one folding-focused objective, same seed; more seeds only if either produces a stable discrete fold. Log early-stage confidence/contact trends and stop clearly weak trajectories before costly later stages. Suggested stop rule pLDDT<.7 at the end of a discrete stage is a prespecified engineering filter, not a law of binding. Do not compare weights alone without considering loss definitions/scales. No evidence that pLDDT weight.1 itself is a software bug.

###C. Fixed-sequence validation gate
Separate fresh model object, no binder coordinate template, two independent models and seeds, binder alone, intact human trimer. Export current aux, inspect allMODEL pixels and distances, record same hotspot and whole-target iPAE separately. Require a stable folded shape and physically credible multi-residue interface. Confidence gates are selection heuristics, never affinity claims. If the pilot still collapses, abandon that scaffold/epitope branch rather than relaxing filters.

###D. Only afterC: species and pH
Validate exact mouse sequence and construct mapping independently. Select interface contacts conserved across species. Then design a few His networks on our own de novo scaffold, evaluate atom-level geometry/protonation-aware alternatives at7.4 and6.0 with stated assumptions, and retain no-His controls. Core-destabilization and interface-repulsion are separate hypotheses. Experimental binding/off-rate measurements are required to claim a switch.

## Alternative choices

1. Recommended: repair pipeline + one bounded pilot, without promising an improved submission before deadline.
2. Lower-value fallback: computationally audit exact original top20 full-trimer sequences for potential salvage, which needs new authorized compute and may simply eliminate unsupported candidates.
3. Not recommended: more His combinatorial scanning of the same crop pose, lengthening the stagnant366 hard stage, calling a threaded mouse score cross-reactivity, or presenting confidence as pH affinity.

## What changed and what did not

Produced source patch, offline geometry checks, CPU mock export regression and strict JSON changes. All modified Python files compile. Real366 allMODEL audit and three negative chain/sequence tests passed. GPU inference on patched scripts has NOT run. No new candidate, ranking change, account re-entry, spending or submission. Original20 remain preserved.

## Sources

Official target/conditions/de-novo restriction:
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha

Versioned source for export/best state and masks:
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/utils.py
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/design.py
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/prep.py
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/loss.py
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/mpnn/model.py

Method/model diversity and hardware guidance, not a TNF success guarantee:
https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft
https://github.com/martinpacesa/BindCraft/blob/main/README.md

pH preprint, version1; TNF effects at5.4 vs7.4, not competition6.0:
https://www.biorxiv.org/content/10.1101/2025.09.29.678932v1.full-text

Actual coordinate source:
https://files.rcsb.org/download/1TNF.pdb

Repository source:
https://github.com/pepert03/protein-binder-design
