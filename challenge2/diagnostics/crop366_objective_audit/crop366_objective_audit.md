#366 objective and trajectory audit

The actual366 launch weights were con.5, i_con1, i_pae.5, pae.1, plddt.1. Optimization used one model, zero recycles, seven chainA hotspots. Independent full-trimer prediction used two models, three recycles and no hotspot mask.

Version1.1.1 binder loss assigns target_id to hotspots when present. It uses that same mask for both interface contact loss and iPAE. Consequently the optimized.8415 is a binder-to-seven-hotspots metric, while independent.9133 averages binder-to-entire-target pairs. They are NOT identical metrics or a clean numeric measure of worsening on rescore. Both remain poor interface-confidence evidence alongside the coordinate/contact audit.

All140 logged steps show pLDDT.30-.36 and iPAE.82-.88. The ten discrete hard-stage steps131-140 are identical at log precision: loss7.90, iPAE.84, ipTM.16, pLDDT.30. No strong candidate appears hidden in the rounded trajectory. Best loss5.12 occurred at soft step35 with pLDDT.30 and ipTM.21. This is not a validated alternative discrete sequence.

Hypotheses for a future controlled test: interface/contact losses dominate the weak folding weight; zero recycles and single-model optimization may limit stable folding; the discretization stage stagnates. These are hypotheses, not demonstrated causes. Change one factor or use a balanced factorial pilot, stop early when binder confidence remains low, and independently validate a fixed discrete sequence on the intact trimer. Do not treat more steps of the same stagnant hard sequence as new evidence.

Sequence extraction defaults to saved best aux; design_3stage saves best in hard stage. Here all hard steps have the same printed scores, so no evidence of a strong hard-stage design missed by only reporting the last log.

Sources freshly fetched:
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/loss.py
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/design.py
Actual fresh_design.py and crop366.log in final366 archive.
