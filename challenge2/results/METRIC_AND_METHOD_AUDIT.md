# Metric and method audit (10 October 2026, 22:02 CEST)

## The apparent hit is a normalized loss, not sub-Angstrom PAE

Upstream ColabDesign get_pae_loss computes get_pae(outputs)/31.0, symmetrizes the matrix, and averages over the specified masks. In the binder protocol i_pae uses binder-target masks. Thus, if the installed Kaggle implementation is the same, i_pae 0.365 represents approximately 11.315 Angstrom averaged interchain PAE, and the previous arbitrary gate i_pae<0.4 is approximately PAE<12.4 Angstrom. The original parent 0.391 corresponds to 12.121 Angstrom. These numbers do not establish binding, and this gate is not a published stringent acceptance filter.

The exact installed function has not yet been inspected. read_seeded_status_and_metric_units.py records it on the next permitted hourly Kaggle check. All conversions above are conditional on matching that function.

Source inspected: https://github.com/sokrypton/ColabDesign/blob/main/colabdesign/af/loss.py

## Comparison with published pH design work

The 2025 preprint reports interface-histidine design filtering at pae_interaction<8 and pLDDT>88. Those are not automatically interchangeable with this cropped ColabDesign aggregate, but they show why a normalized 0.4 should not be described as a strong hit.

The preprint tested TNF-alpha pH-sensitive designs through buried His-His and His-Arg networks in the binder core. For TNF-alpha, 48/72 retained binding, 7 had affinity tighter than 25 nM and 3 of those were pH-sensitive; the top design had 79-fold weaker binding at pH 5.4 than 7.4. This is not evidence that our interface-His candidates work, nor that the same effect holds at the competition's pH 6.0.

For interface-mediated switches, the work uses protonation-aware hydrogen bonding and electrostatic binding-energy comparisons, not His counts alone. Global His bias yielded little pH sensitivity in initial tests. Its pH-sensitive EphA2 and TNFR2 designs had two or more His-cationic contacts. Our hisCat proximity count is a limited heuristic, and target TNFR2 is not TNF-alpha.

Source inspected: https://www.biorxiv.org/content/10.1101/2025.09.29.678932v1.full-text
Official assay conditions: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha

## Missing self-consistency checks

BindCraft explains that hallucinated AF2 sequences can be difficult to purify. It uses ProteinMPNN sequence optimization, complex reprediction, binder-alone reprediction and interface/biophysical assessment. The soluble MPNN weights are relevant for repairing the present hydrophobic surface, but preservation of specific electrostatic contacts must be checked after sequence optimization.

It also warns that overly trimmed target structures can expose core hydrophobics and should preserve the realistic binding mode. Our 110-residue crop from two TNF protomers leaves the full trimer and many residues untested.

Sources inspected:
- https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft
- https://github.com/martinpacesa/BindCraft/blob/main/README.md
- https://github.com/sokrypton/ColabDesign/blob/main/af/examples/peptide_binder_design.ipynb
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9312799/

## Next bounded experiment

1. Let the current seeded run finish; retain completed sequences and structures, not only logs.
2. Confirm installed PAE normalization and design weights; report units explicitly.
3. Inspect parent and new-candidate structures to identify interface versus solvent-exposed tail positions. Compute solvent accessibility and contacts, and check full-trimer clashes after target alignment.
4. Redesign only eligible noninterface positions with soluble ProteinMPNN weights, preserving contact residues and selected His sites. Repredict candidate complex and binder alone. Keep unmodified controls.
5. Independently predict human and mouse full-trimer complexes for the best repaired candidates. Threading mouse sequence onto a human crop remains only a proxy.
6. Examine whether the His sites plausibly change electrostatic interactions at pH 6.0, and quantify the limitations. No model-only calculation should be presented as measured pH-selective affinity.

These are planned tests. No new independent validation or experimental result has been completed. No designs have been submitted.
