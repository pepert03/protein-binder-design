> Historical planning/results text below is superseded by [LATEST_STATUS.md](LATEST_STATUS.md), updated11October2026. No supported pH-switch binder or submission.

# TNF-alpha pH-switch binder: method text (fallback, draft)

Goal: bind TNF-alpha trimer at pH 7.4, lose binding at pH 6.0, human and mouse.

## Mechanism (sourced: Baker lab preprint, bioRxiv 10.1101/2025.09.29.678932v1, fetched 10 Oct 2026)
1. Interface repulsion. Place His on the binder next to Lys/Arg/His on the target. At pH 7.4 His is neutral and can accept an H-bond; below its pKa (~6.5 per the paper) it is protonated, positive, and repels the nearby cation. In the paper's EphA2 and TNFR2 sets, every pH-dependent binder had 2 or more protonated-His to cation contacts.
2. Core destabilisation. Buried His-Arg / His-Lys / His-His H-bond network in the binder core, so protonation destabilises the fold at low pH. The paper reports TNF-alpha binders made this way (Fig 3B, His-His and His-Arg networks).
3. Caution from the paper: simply raising global His bias in ProteinMPNN gave more His but little pH sensitivity. His must be positioned.
Numbers in the paper: TNFR2 binder set, 43 designs tested, 42 expressed, 6 had >2-fold weaker binding at pH 5.4, best 122-fold. Tested at pH 5.4, not 6.0, so effect at 6.0 will be smaller (inference, not measured).

## Our pipeline
- Target: TNF-alpha trimer (PDB 1TNF/3WD5), epitope crop around the receptor-binding groove, hotspots A65,A67,A113,A115,A141,A143,A146.
- Backbones: AF2 hallucination (ColabDesign binder protocol).
- Sequence: ProteinMPNN with His bias on binder, target fixed.
- Filter: AF2 i_pae, i_ptm, pLDDT; require at least 2 binder His within 5 A of target K/R/H; reject His-rich low-complexity sequences.
- Rank up to 20 by: His-cation contact count, i_pae, pLDDT; mouse conservation of contacted target residues (to check).

## Status and limits (honest)
- Best AF2 self-scores so far: i_pae 0.498 / i_ptm 0.389 (hallucination, no His); His-MPNN redesigns i_pae ~0.80, i_ptm ~0.1. Not yet at a usable bar (i_pae <0.4 target).
- AF2 scores are not binding or pH measurements. No pKa or electrostatic (ddG_elec) calculation has been run.
- Mouse cross-reactivity not yet checked.
- Not submitted.
