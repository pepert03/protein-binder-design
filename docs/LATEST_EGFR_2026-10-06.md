# Latest verified EGFR results, 6 October 2026

## Completed
- Human20: 20 original CSV sequences against human EGFR P00533 25-645 (621 aa). Single-sample Boltz2.2.1, empty MSA.
- Mouse top5: 5 original sequences against Q01279 25-647 (623 aa), saved run 44m37s.
- Mouse diverse4: r02_l05, rank08_native, rank01_native, r03_l05, same mouse target, saved run37m02s.
- Exact binder/target identities, archive integrity and recorded file hashes were checked for all29 models.
- Sequence lineage, local SwissProt screen, narrow EGFR structural similarity, human/mouse contacts, experimental human receptor/glycan context and top5 human PROPKA analyses are recorded in the result directories.

## Current interpretation
There is no validated combined human/mouse winner. Scores are interface confidence, not binding probabilities, affinity or pH selectivity. r02_l07 remains a changed-epitope hypothesis. r01_l05/r02_l03 retain closer poses but known human clash/context risks. rank01_native/r02_l05 are conditional diverse reserves with substantial pose changes. rank08_native has mouse ipTM0.759 but moves37.8A and shares only2 receptor contacts. r03_l07/r03_l05 lose combined priority. Eleven of the human20 are not modeled in mouse. The new mouse poses have not been checked against experimental receptor/glycan context.

All20 remain officially novelty-unverified. No experiment, KD or acid-selective binding is established. No final list approved or submitted. Audit CSV numbers are identifiers of the earlier human review, NOT a submission ranking.

## Not completed
Lightning robustness was abandoned because the platform was unreliable; no3-sample robustness outputs exist. The earlier statement that this run was pending is superseded. No restart is planned without a new request. DGX was unsuitable under its host-RAM limit. Alternative-account GPU use was dropped. These attempts produced no valid additional folds.

## Repository artifacts
- `results/boltz2_csv20/structures/`: all20 human CIFs, inputs, provenance and records; existing pH analysis remains in this directory.
- `results/mouse_top5/` and `results/mouse_diverse4/`: all9 mouse CIFs, exact inputs, provenance and run records.
- `results/output_archives/`: three complete recovered output archives, with SHA256SUMS. No weights or model cache included.
- `code/kaggle_mouse/`: exact notebooks used for the completed mouse runs.
- `results/egfr_audit_mouse9/`: current combined table, source notes, methods and corrected audit CSV. Earlier audit directories are historical snapshots.

## Reports and saved run sources
Private current audit: https://files.instinct.com/file-01M48RS05JMG9Q33TTR9B297QT
Human20 output: https://www.kaggle.com/code/peperidruejo/proteinbase-csv20-full-human-egfr/output?scriptVersionId=355767895
Mouse top5: https://www.kaggle.com/code/peperidruejo/proteinbase-mouse-egfr-top5-audit
Mouse diverse4 output: https://www.kaggle.com/code/peperidruejo/proteinbase-mouse-egfr-diverse4-audit/output?scriptVersionId=355797295
Official target/rules: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/egfr
Novelty scale: https://www.adaptyvbio.com/blog/novelty

Historical REPORT.md and STATUS.md are preserved. This note and mouse9 data supersede their earlier target, shortlist and monitoring statements. Repo remains private. Nothing submitted or publicly shared.
