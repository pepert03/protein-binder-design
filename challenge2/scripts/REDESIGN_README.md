# Bounded sequence repair and reprediction

Prepared scripts, syntax-checked locally, NOT YET RUN on Kaggle. No claim of successful generation or validation.

At the next allowed check, inspect the completed seeded run's chain identities, duplicate atoms, source sequences and installed APIs before choosing a PDB. Do not run GPU predictions while the seeded process is still active.

mpnn_surface_repair.py uses exact structural chain IDs. It freezes the target, all binder residues with any heavy atom within 8 A of target, binder-only SASA below 30 A^2, and existing histidines. It samples eight sequences. This conservative mask changes only putative noninterface surface residues. The 8 A and 30 A^2 cutoffs are working choices, not validated acceptance thresholds. It prohibits new cysteines but keeps fixed ones. It asserts every fixed residue is preserved.

The constructor is inspected: when soluble weights are supported they are requested, otherwise the original packaged model is used and reported explicitly. ColabDesign v1.1.1 has original weights only. Do not assume a current upstream API on that installation.

repredict_repair.py takes the unmodified control and at most two distinct repair samples selected first by reduced longest-window hydrophobic burden and then MPNN score. Neither criterion is a binding score. It predicts binder alone without structural template and a complex against the given target structure, with two models and three recycles. Binder CA RMSD after Kabsch alignment is recorded, with PAE implementation source for unit verification. The target input can be full human trimer, once current memory/runtime constraints are checked. Mouse validation remains a separate step.

Sources inspected:
- https://github.com/sokrypton/ColabDesign/blob/main/mpnn/README.md
- https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/mpnn/model.py
- https://github.com/sokrypton/ColabDesign/blob/main/colabdesign/mpnn/model.py
- https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/design.py
- https://github.com/sokrypton/ColabDesign/blob/main/colabdesign/af/prep.py
- https://biopython.org/docs/latest/api/Bio.PDB.SASA.html

The scripts have not been integration-tested. Inspect generated coordinates visually before accepting any structure. No automatic submission or top20 replacement occurs.

Preparation update 23:16: raw duplicate-ATOM identifiers are rejected before Bio.PDB parsing because the previous run had duplicate binder atoms that caused NaNs. No automatic deduplication is performed in the repair script. A duplicated file must be recovered and inspected first. Browser daily time limit deferred the 23:15 live check until midnight. Scripts still not executed.
