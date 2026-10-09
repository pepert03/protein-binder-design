# Challenge 2 research summary (checked 2026-10-09)

Full explainer (Spanish, from zero):
https://files.instinct.com/sptoiyocdtrq-challenge-2-tnf-desde-cero

## Verified against official sources
- Challenge page: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tnf-alpha
  (target, pH directions, ranking, deadline, 20-design limit).
- Previous EGFR winners: Martin Pacesa and Lennart Nickel with BindCraft
  (https://www.adaptyvbio.com/blog/po102/ and /po104/): 91.5 nM best de novo in round 2,
  6 of 7 de novo leaderboard spots; hallucination beat diffusion.
- BindCraft: https://github.com/martinpacesa/BindCraft (CUDA GPU needed; 32 GB or more
  recommended for large complexes; Colab available).

## Other participants (GitHub search sample, not exhaustive)
Their READMEs are claims, not verified results.
- wwwtalalahmed98-collab/ph-switch-binders: most developed Challenge 2 method
  (epitope first, AfDesign/ColabDesign hallucination, independent re-prediction,
  calibrated histidines).
- harishrajaram-svg/anthropic-adaptyv-2026, SyzenShen/adaptyv-binder-challenges,
  jimjmccann/CCNY-McCann-Koder: multi-challenge repos.

## OpenFold Slack (technical-updates, 2026-10-05)
Deterministic inference, seed fix, lower memory use. Note: the connected workspace is
OpenFold, not the Proteinbase competition Slack.

## Unverified
Wet-lab outcomes for Challenge 2 do not exist yet. Mouse accession P06804 comes from a
participant repo, not confirmed on the official page.
