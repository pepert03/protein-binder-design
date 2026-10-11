# Paired folding pilot: complete, no successful binder
11October2026,23:57CEST

Same target366/binder50/seed32101/model1ptm/zero design recycles. Baseline pure-intra weight0 and pLDDT.1; folding pure-intra.4 and pLDDT.5. No histidine bias, but histidines allowed. Joint intervention, one seed, not causal isolation of either weight.

| Metric | Baseline | Folding |
|---|---:|---:|
| Discrete design binderpLDDT | .30155 | .28223 |
| Discrete binder-intraPAE normalized | .58613 | .58640 |
| Independent full-trimer binderpLDDT mean | .30935 | .31138 |
| Independent full-trimer iPAE mean normalized | .91329 | .91513 |
| Independent full-trimer ipTM | .15971 | .15803 |
| Binder-alonepLDDT mean | .58910 | .51287 |
| Binder-alonepTM | .26865 | .13668 |
| FullMODEL1/2 minheavy distanceÅ |4.364/6.036|5.854/11.571|
| FullMODEL1/2 binderresidues<4Å |0/0|0/0|
| Alone/fullCA shapeRMSDÅ |13.28/10.00|8.99/8.77|

Baseline exactly reproduces earlier366 sequence. Folding sequence WWEDQIYGHDMGIWNLSIFDNEANQWKRDPNDSTGNAANRTHHGMYFHYG. Relative shapeRMSD improvement is insufficient: still large, lowconfidence, extended shapes and no<4Åinterface. Actual PDB pixels inspected for both full/alone models. Optimized target-fit~.637Å; no omitted native target atom pairs<2Å for eitherarm.

Fixed-sequence predictions used2models and2seeds; seed32102/32103PDBs byteidentical with dropout off. Count two distinct models, not four independent stochastic trials. Optimization hotspot iPAE differs in scope from fulltarget mean. Exact rawPAE scopes/arrays exported. No affinity, mouse orpH evidence.

Conclusion: the callback is mathematically/GPU operational but this particular folding-emphasis intervention does not yield a usable binder. Do not add either result to the final shortlist. It does not rule out different weights, scaffolds, multimer/ensemble design or larger-memory resources. Next useful branch should change generation/model/scaffold/epitope strategy, not continue this weak hard-stage sequence or premature His screening. Needs a new plan and fresh resource access approval, not automatically launched.

Both arms completed, no pilot process at final export. Actual inputs/stagePDB/rawPAE/full+alonePDB/log/config preserved in folding_pilot_final.tar.gz. Final local perMODEL/omitted geometry and shape analysis in pilot_final_analysis.json. Browser released, rescue wake cancelled. GPU allocation itself left to expire at01:14:49, no newjob. Original20 unchanged, no submission or spending.
