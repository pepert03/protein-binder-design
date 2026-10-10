# Draft top 20 for TNF-alpha challenge (DRAFT, NOT SUBMITTED)

All scores are AlphaFold2 self-scores from ColabDesign hallucination (i_pae lower is better, bar used 0.40). They are NOT binding, pH 7.4 vs 6.0, or mouse affinity data. Mouse column = mouse TNF (UniProt P06804) threaded onto the human 1TNF crop, fold-survival screen only.

Ranking rule: ranks 1-16 are traj6 His variants ordered by the worse of the human and mouse i_pae (the challenge scores pH selectivity first, then mouse cross-reactivity, then human affinity; none of those are modelled, so this is a proxy). Ranks 17-20 are the parent traj6 and the next-best other backbone family for diversity. Every sequence was checked: 65 aa, His variants differ from traj6 only at the listed His positions, all 20 are unique.

Caveats: (1) 16 of 20 share one scaffold, so diversity is low. (2) Ends in a leucine-rich hydrophobic tail, likely sticky or aggregation-prone. (3) His placement is by proximity to target Lys/Arg, not a computed pH response. (4) No expression, binding or pH data exist for any of these. (5) Not submitted; submission needs José's go.

| rank | His pos (0-based) | human i_pae | mouse i_pae | sequence |
|---|---|---|---|---|
| 1 | [8, 22] | 0.381 | 0.394 | `CKDVIIIKHPDSFQPDPRWHQTHQRYRARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 2 | [6, 22] | 0.388 | 0.398 | `CKDVIIHKVPDSFQPDPRWHQTHQRYRARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 3 | [6, 26, 30] | 0.388 | 0.4 | `CKDVIIHKVPDSFQPDPRWHQTYQRYHARAHGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 4 | [6, 8, 26] | 0.392 | 0.4 | `CKDVIIHKHPDSFQPDPRWHQTYQRYHARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 5 | [6, 22, 30] | 0.389 | 0.403 | `CKDVIIHKVPDSFQPDPRWHQTHQRYRARAHGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 6 | [22, 26] | 0.397 | 0.406 | `CKDVIIIKVPDSFQPDPRWHQTHQRYHARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 7 | [6, 29] | 0.407 | 0.371 | `CKDVIIHKVPDSFQPDPRWHQTYQRYRARHEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 8 | [5, 30] | 0.387 | 0.408 | `CKDVIHIKVPDSFQPDPRWHQTYQRYRARAHGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 9 | [6, 18] | 0.39 | 0.41 | `CKDVIIHKVPDSFQPDPRHHQTYQRYRARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 10 | [17, 22, 26] | 0.405 | 0.412 | `CKDVIIIKVPDSFQPDPHWHQTHQRYHARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 11 | [1, 6, 26] | 0.409 | 0.417 | `CHDVIIHKVPDSFQPDPRWHQTYQRYHARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 12 | [8, 29, 30] | 0.421 | 0.408 | `CKDVIIIKHPDSFQPDPRWHQTYQRYRARHHGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 13 | [1, 6, 25] | 0.405 | 0.425 | `CHDVIIHKVPDSFQPDPRWHQTYQRHRARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 14 | [18, 30] | 0.395 | 0.429 | `CKDVIIIKVPDSFQPDPRHHQTYQRYRARAHGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 15 | [17, 18] | 0.402 | 0.429 | `CKDVIIIKVPDSFQPDPHHHQTYQRYRARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 16 | [5, 6, 26] | 0.403 | 0.431 | `CKDVIHHKVPDSFQPDPRWHQTYQRYHARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 17 | run6 traj6 parent | 0.391 | not scored | `CKDVIIIKVPDSFQPDPRWHQTYQRYRARAEGQDGGGDGVVVVLLCCLLLLLLLILLLKLLKKFD` |
| 18 | [28] | 0.473 | not scored | `DDDKIAAILLCLLLLLVLLTPAGPQPSEHDVADYRKYREKQLKAPRKAGKACESCIKEIIEKKDD` |
| 19 | [26] | 0.503 | not scored | `DDDKIAAILLCLLLLLVLLTPAGPQPHETDVADYRKYREKQLKAPRKAGKACESCIKEIIEKKDD` |
| 20 | run4 traj8 parent | 0.493 | not scored | `DDDKIAAILLCLLLLLVLLTPAGPQPSETDVADYRKYREKQLKAPRKAGKACESCIKEIIEKKDD` |
