# Independent prediction and full-trimer steric audit (11 October 2026)

The preserved seeded optimization sequence does not retain its 0.365 normalized i_pae under independent 2-model, 3-recycle prediction. Two repeats give approximately 0.432, equivalent to 13.4 A under the installed ColabDesign /31 normalization.

| sequence | complex normalized i_pae | complex i_ptm | alone pLDDT | alone/complex binder CA RMSD (A) |
|---|---:|---:|---:|---:|
| control | 0.432335 | 0.614841 | 0.744533 | 4.11671 |
| repair03 | 0.413969 | 0.586119 | 0.853668 | 5.31077 |
| repair02 | 0.410512 | 0.583711 | 0.845388 | 4.36431 |

Sequence repair reduces the explicit hydrophobic-window heuristic but creates low-complexity E/R-rich tails. Higher predicted confidence alone is not measured folding, solubility, binding or pH selectivity. Whole-binder self-consistency is not improved. Actual six prediction ensembles and logs are preserved separately in the comparison archive.

## Full-trimer compatibility warning

Predicted target-crop CA atoms were aligned to the exact original 110-residue crop from 1TNF. Exact target sequence identity and ordering were verified. Then transformed binder heavy atoms were checked against full-trimer atoms omitted from the crop (including protomer C).

Control has 466 binder/omitted-target heavy-atom pairs below 2 A, with 30 binder residues affected and a minimum 0.20 A distance. Repairs have 481 and 479 such pairs, minima 0.33 and 0.29 A. CA target-fit RMSD is 2.17-2.26 A, so exact counts depend on alignment, and the screen is not relaxed binding energy. Nevertheless, widespread deep overlap makes these particular predicted poses unsuitable as evidence of viable full-trimer binding.

Visual inspection: full_trimer_pose.png shows the crop-derived control binder running into the intact trimer. compare_ca_traces.png shows all three as extended helices with context-dependent bends. These are CA traces, not full side-chain validation. This audit concerns the independently regenerated control and its two repair poses, not the lost original optimized trajectory PDB. It does not establish that every top20 candidate is impossible.

Next independent test should use the intact trimer. Continuing to optimize only this cropped pose risks improving an artifact. No designs have been submitted and the existing top20 remains unchanged, not endorsed as ready.

Structure source: https://files.rcsb.org/download/1TNF.pdb
Official soluble trimer requirement: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha
Trimming caveats: https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft
