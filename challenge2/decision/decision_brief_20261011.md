# TNF-alpha: decision evidence, 11 October 2026

Recommendation: do not treat the existing top20 or the new candidates as validated TNF pH-switch binders. No tested new candidate currently has convincing intact-trimer interface support. Do not replace the current ranking with the366 candidate. This recommendation is about the evidence, not proof that every original sequence fails experimentally. No submission made.

## Observed results

| Test | Independent intact-human-trimer iPAE, normalized | ipTM | Binder residues within4 Å |
|---|---:|---:|---:|
| Unmodified65aa control | .913904 | .144537 | 0 |
| MPNN repair02,65aa | .913075 | .156276 | 0 |
| Fresh277crop,50aa | .900203 | .171677 | 0 in both models |
| Fresh366crop,50aa | .913287 | .159711 | 0 in both models |

The installed ColabDesign code divides PAE by31. Values near.913 therefore mean about28.3 Å, not.913 Å. Confidence scores are not measured affinity or pH response. The366 optimization itself was weak: pLDDT.30155, ipTM.15636, normalized iPAE.84153. Its intact predictions change binder pose, lack a stable folded binder interface, and have minimum target-binder distances4.36 and6.04 Å. No severe<2 Å clash is a necessary geometry check, not positive binding evidence.

Earlier crop poses could overlap omitted trimer atoms. The366 neighbourhood avoids that particular omitted-atom problem, but does not rescue interface quality. Full456-residue gradient design failed on KaggleT4 and the actual16GiB DGX MIG slice, including40aa binder and patched subbatch1/bf16/rematerialization.366+50 fits; fitting in memory is not design success.

## Original top20

The exact20 have not collectively been independently validated against intact human trimer. The tested control is not interchangeable with those20. Their human/mouse columns are crop model proxies, and mouse was threaded onto human structure. Histidine proximity is a pH-switch hypothesis, not calculated or measured selectivity. Sequence-only liabilities and diversity should be considered separately from binding.

The attached deterministic static audit checks exact sequences in the preserved top20 draft. It does not predict solubility, expression, cross-reactivity or affinity. The ranking and sequences remain unchanged.

## What remains unproven

- Binding to soluble human TNF-alpha trimer.
- Binding to mouse TNF-alpha at pH7.4.
- Binding at pH7.4 with no detectable binding at pH6.0.
- Expression, aggregation behaviour and experimental folding.

The official challenge ranks pH selectivity first, mouse cross-reactivity second, affinity third. Its recommended target is the soluble trimer, especially the adjacent-protomer receptor-binding epitope. CSV requirements include name, sequence and molecule_class. A format-correct file does not fix unsupported biology.

## Next useful branch

A credible restart needs a new folded scaffold/interface, not further His substitutions of a poorly supported pose. Preserve whole-trimer geometry in final validation, independently rescore fixed sequences, reject omitted-atom clashes, and only then investigate histidine/protonation changes and mouse interface conservation. Larger-memory gradient compute might enable full-target optimization, but no such resource is currently verified or authorized. More crop optimization on the available slice is technically possible, with low confidence that one short run changes the decision.

## Evidence

Actual PDB/json/log archives for all tested branches were exported to the parent. Final366 visual and coordinate audit were inspected. Static audit is generated from the unchanged top20_draft.csv.

Official requirements re-read11October18:03:
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha

Deadline: Sunday11October23:59AoE, calculated as Monday12October2026 13:59 Europe/Berlin. Page mentions open submissions while also listing later generic timeline stages; deadline text is clear, but no logged-in submission state was inspected.
