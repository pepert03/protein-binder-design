# Controlled next-pilot specification,11October20:37

Proposed experiment, not launched or approved compute entry. No promise of a deadline-ready binder.

## Metric correction
The existing ColabDesign binder "pae" loss with mask_1d=binder and no mask_1b averages binder rows against ALL columns, including interface. It is not a pure intrabinder-folding loss. Merely increasing its weight under a "folding-first" label would still mix fold and interface. A real intrabinder term needs binder-binder masks; compare it separately with hotspot and entire-target terms. Supplied pureNumPy helper reports explicit scopes/units/pair counts from one raw Angstrom matrix at a time; it is NOT wired into AF optimization and raw matrices from earlier runs were not exported. Existing aggregate numbers cannot recover these scopes.

##Fixed experiment factors
-50aa binder and366crop as baseline only if limited to the verified16GBslice; geometry and gaps remain explicit limitations. No new length change in the same weight comparison.
-Same seed32101, model_1_ptm, zero optimization recycles, subbatch1/chunk32/bf16/remat as prior memory-fitting baseline. Do not silently add multimer/recycles and claim an objective-only experiment.
-No histidine bias; cysteine excluded as before. Record complete weights and exact source/parameter hashes.

##Intervention
Baseline vs one arm that adds a pure binder-binder PAE term and raises binder pLDDT weight. These are a joint "folding emphasis" intervention, not independent attribution of each term. Prespecify weight values after checking differentiable-loss wiring against installed source. Alternatively test each change separately if resources permit. pLDDT.1 by itself is not a software bug. There is no evidence for a guaranteed optimum weight.

##Run and stop
At most one seed perarm initially. Export after soft,temperature and hard stages with identical chosen aux for sequence/PDB/log. Retain raw perMODEL PAE matrices to report binder-intra/hotspot/all-target scopes. Current get_seqs best vs aux last semantics must not be mixed.
Use early folding/contact checks to avoid blindly spending the full140steps on a clear weak trajectory. Discrete pLDDT<.7 is a proposed engineering stop rule, not experimental binding law. If baseline weak but intervention becomes stable, proceed to fixed-sequence validation; otherwise stop this pilot and change scaffold/target branch rather than relax acceptance metrics.

##Validation
Fresh prediction object, no binder template, independent model(s)/seed(s), binder-alone and intact456human-trimer predictions. Use explicit get_best=False, verify PDB exact sequence, allMODEL contact/clash analysis and pixels. Scores are filters, not affinity. Audit omitted full-target atoms even for a valid-looking crop pose. Mouse exact native construct and pH7.4/6.0 work come only after fold/interface gate.

##Resource alternative
Larger-memory multimer ensemble design is a separate plan, not part of the paired objective test. Verify resource access/availableVRAM before proposing launch. BindCraft32GB/~550residues guidance is approximate. No admin request, spending, new account or DGX authentication started.

Source definitions:
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/loss.py
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/utils.py
