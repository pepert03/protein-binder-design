# TNF restart preparation,12October2026

Status: local preparation only. No new resource access or GPU run. Official deadline freshly read:18October23:59AoE =19October13:59CEST.

Distinct route: ColabDesign AF2-multimer, local277-residue neighbourhood around functional site,60-70aa binders, multi-seed early rejection. Not full BindCraft. Maintainer recommends>=32GB for full workflow, so16GiB fit is unknown. Do not repeat failed monomer366+50 trajectory.

## Preflight, before investing a long trajectory
- Inspect actual device memory and available multimer_v3 parameters. Do not buy resources, request admin, use another person's account or silently replace missing models with monomer.
- Use preserved neighbouring protomer atoms in neighbour24.pdb;277 residues(A117/B56/C104), same original receptor-site hotspots A65,A67,A113,A115,A141,A143,A146. Crop is a generation expedient, not intact-trimer validation.
- After prep_inputs, set target protomer asym/entity/sym IDs using chain_ids.py. Baseline ColabDesign1.1.1 _prep_binder uses[target_total,binder_length] and collapses target asym IDs. CPU mapping tests passed; propagation into multimer features remains GPU-unverified.
- Try one60aa binder, one multimer_v3 model, zero-gradient-recycle memory preflight first. If fit, test1recycle before a full trajectory. Cap first fit attempt; stop on OOM rather than repeatedly waste allocation. Saved params must contain required multimer models.
- Verify native target residue mapping, chain IDs, pure binder metric mask and callback gradients. Save exact config/source/input hashes.

## Generation if preflight passes
Use fresh seeds rather than resume either failed sequence. Small initial cohort, stop early when best selected state's binderpLDDT<=.65 at50logit steps. This.65 is a BindCraft stage-screening heuristic, not proof of binding. Then remaining25logit/45temp/5hard steps only for survivors, with renewed fold gate at eachstage. Do not call soft sequence confidence a discrete-sequence validation. Saved sequence/PDB/metrics must come from the same selected aux. Sample available trained models if verified in memory; do not count deterministic duplicate seeds as independent evidence.

Provisional model/contact settings informed by BindCraft source: pLDDT.1, pure binder-intraPAE.4, interfacePAE.1, con1,i_con1; intra contacts2/cutoff14Å/seqsep9, inter contacts2/cutoff20Å. These are an explicitly modified ColabDesign objective, not an exact reproduction of BindCraft: installed stock pae includes binder-interface pairs, callback isolates binder-intra. No proof these settings work for TNF. Consider radius-of-gyration only after implementing and testing it, not implicit. Per-term scalar losses do not establish gradient dominance.

## Validation before recommending any candidate
1. Fixed discrete sequence binder-alone predictions, distinct trained models,3recycles. Check fold confidence, shape and sequence/backbone match. Reject extended/non-globular lowconfidence shapes.
2. Independent intact1TNF456-residue trimer predictions with properly distinct protomers and binder. Retain per-model rawPAE/pLDDT plus binder-interface geometry, native target-fit/omitted-atom audit. Do not accept crop-only contacts or pooledMODEL contacts.
3. Visually inspect actual full/alonePDB pixels. Save evidence regardless of success.
4. Require fold and physical interface together. No thresholds or AF2 scores establish affinity. No pH response or mouse binding claim from histidine counts. Only after a supported human fold/interface candidate should mouse conservation/structural checks and bounded pH hypotheses proceed.

No submission, shortlist change, spending or resource launch authorized by this plan. New DGX sign-in/session needs explicit user approval relayed by main. Check/export logged-in results at most hourly unless user agrees to faster checks and risk. Report progress from retained state every30minutes while an actual experiment is active.

Sources freshly fetched:
- https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha
- https://github.com/martinpacesa/BindCraft/blob/main/README.md
- https://github.com/martinpacesa/BindCraft/blob/main/functions/colabdesign_utils.py
- https://github.com/martinpacesa/BindCraft/blob/main/settings_advanced/default_4stage_multimer.json
- https://github.com/sokrypton/ColabDesign/blob/main/af/README.md
- Installed version source previously retrieved: https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/prep.py
