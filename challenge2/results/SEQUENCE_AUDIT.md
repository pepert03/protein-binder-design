# Sequence audit (10 October 2026)

The existing top20_draft.csv is a research table, not a ready submission: it lacks the required name and molecule_class columns. A separate submission_format_DRAFT_NOT_SUBMITTED.csv preserves the exact 20 unique sequences and ranking, with stable draft names and molecule_class=protein. Nothing has been submitted.

## Reproducible sequence-only warnings

All sequences are 65 residues. Ranks 1-17 contain a 19-residue window made entirely of residues in AVILMFWYC: VVVVLLCCLLLLLLLILLL at positions 40-58 (one-based). These sequences have 14 leucines, 3 cysteines and a longest identical-residue run of 7. Ranks 18-20 peak at 16/19 residues from that set. Exact computed values are in sequence_audit.json.

This count is an explicit screening heuristic, not a validated membrane, aggregation or solubility predictor. It does not establish that the tail is exposed in the folded structure. Still, the unusually long hydrophobic block makes chasing lower AF self-scores without sequence repair a weak strategy.

## Consequence and next tests

1. Keep the originals as controls, not evidence of therapeutic readiness.
2. Before ranking a new seeded candidate, require independent full-trimer complex prediction, binder-only fold consistency and surface/packing inspection. Cropped epitope success may hide clashes with omitted target residues.
3. For solubility repair, redesign noninterface surface positions with ProteinMPNN while preserving interface contacts and selected His sites. Compare against the original. Do not blindly replace buried hydrophobic residues, and do not call an untested edited sequence an improved binder.
4. A pH switch needs a separate protonation-sensitive assessment. Count His near target basic residues only as a hypothesis, not pH evidence.
5. The present mouse screen threads mouse sequence onto the human crop. It is not an independently predicted mouse complex or measured cross-reactivity.

## Sources inspected

- Official target, pH direction, deadline, format and de novo rules: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha
- ColabDesign's example warns that AfDesign-only optimization can find adversarial sequences, and describes ProteinMPNN backbone redesign: https://github.com/sokrypton/ColabDesign/blob/main/af/examples/peptide_binder_design.ipynb
- Study of solubility-aware AfDesign design reports hydrophobicity/solubility problems and investigates sequence-level loss terms: https://pmc.ncbi.nlm.nih.gov/articles/PMC9312799/
- Private source sequences and recorded results: https://github.com/pepert03/protein-binder-design

The last three next-test items remain plans, not completed validation.
