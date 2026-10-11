# Current Challenge2 status,11October2026

Read this before historical RESULTS.md, method.md or run scripts.

No experimentally validated binder, pH switch or mouse cross-reactivity. No submission. Original20 sequences/ranking remain unchanged and are exploratory hypotheses, not accepted designs.

New independent full-human-trimer tests of an unmodified65aa control, MPNN repair02, fresh277crop50aa and fresh366crop50aa do not support a convincing interface. Offline re-audit of all8 PDB models: zero binder residues within4 Å. Scores and distances are in diagnostics/whole_review/. This is computational negative evidence for those exact tested sequences, not proof every original20 fails.

Original20:20 unique65aa sequences;17 share a19aa hydrophobic/Cys-rich window. pH inference was His proximity, not protonation-aware calculation. Mouse screen was sequence threading onto human crop, not independent native mouse-trimer validation.

Installed ColabDesign1.1.1 normalizes PAE by31. Optimization with hotspot masks and full-target validation iPAE are different aggregates. Report scope and units. Do not use old<.4 gate as a validated binding threshold.

Full456+50/40 gradient attempts failed on16GBT4/MIG despite subbatch1/bf16/remat.366+50 fits but140-step binder pLDDT stays.30-.36 and hard stage stagnates. No active GPU run at review time. Bigger memory/multimer designs are proposed, not executed.

Confirmed code risk: save_pdb defaults to saved best; after predict in a model with design history, this can export stale design coordinates. Explicit get_best=False after prediction is required. Proposed patch includes a CPU mock regression and per-model offline geometry auditor. Syntax/mock/offline tests passed; modified AF scripts have not been rerun onGPU.

Recommended path: repair reproducibility and export semantics first; one no-His, folding-first controlled pilot second if fresh compute access is approved; intact-trimer and binder-alone validation before mouse/pH engineering. See diagnostics/whole_review/REVIEW_AND_PLAN.md.

Official target/deadline:
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha
Closes11October23:59AoE =12October13:59CEST. No submission authorized by this report.
