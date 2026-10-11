# Per-model His proximity repair,11October20:06

Historical hiscat in run3/run5 pooled all ATOM records across MODEL blocks. New his_contacts.model_counts computes each MODEL independently and exposes contacts separately. Legacy integer hisCat is the minimum count across all models, with hisCat_rule recorded. This is a changed heuristic, not comparable to pooled historical hisCat. Old outputs remain untouched.

Explicit binder/target chains can be supplied. Legacy wrapper infers last chain as binder in AfDesign exported layout and validates binder length. Future general use should pass roles explicitly. Incomplete histidine rings fail instead of silently counting.

CPU regression creates two models with no within-model contacts but spurious cross-model proximity. New helper returns zero/zero. One-contact/zero-contact case gives conservative0; single-contact model gives1; wrong binder length fails. Real control,repair02 and366 intactPDBs were inspected by helper. These are proximity counts, NOT hydrogen bonds, pKa, protonation, binding affinity or pH selectivity. No GPU execution; no ranking/sequence change.

Run3/run5 import this helper alongside their scripts. Compile and regression passed locally. Source patch assumes the preceding whole-review patch (790fdccb) is applied. Do not execute historical scripts as an automatic design plan.
