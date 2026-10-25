# Provenance scaffold_v6 + MPNN v6 (2026-10-25)

## Busqueda de pose (muestreo_v6.py)
- Receptor: /tmp/ch3pack/challenge3/target/6X18_RP.pdb (6X18 cadena R; 384 residuos backbone-completo de 404; 3138 atomos)
- Seeds: SEED=20261025 en todas las etapas; etapas: global laxo 4400 (0 supervivientes duros), refine_v5 800 (4 superv, mejor 19), refine_wide 1500 (23 superv, mejor 9), refine2 1200 (306 superv, mejor 4)
- Gate geometrico: 0 choques bb <2.0 A; n_bols>=2 y d_ecd<=8 en etapas de refinamiento (RELAJADO respecto al corte original v5: n_bols>=3 y d_ecd<=6.0 — la pose final tiene n_bols=5 OK pero d_ecd=6.02 > 6.0: INCUMPLE el corte original por 0.02 A; se reporta tal cual, sin mover el umbral retroactivamente)
- Packing en bucle: FASPR commit 0d55732fd6307f373018c6bddd842291c355c5f7, lib dun2010bbdep.bin sha256 ed3f7be5f33b5fa947ac5e83cb024c6a6af6440bb50a1c8073aacabe6d792d0e; proxy poly-ALA (NO mejor caso universal: Gly no tiene CB); receptor fijo por minusculas; choques contados vs coords receptor originales; criterio solape vdW >0.4 A (Bondi); control nativo observado (sin packing): 3 choques, dmin 2.46 A. Umbral <=5 EXPLORATORIO.
- Pose v6: anclaje [109.981,134.129,111.779], eje [0.2222,0.7393,-0.6357], roll 2.8347

## Artefactos
- scaffold_v6.pdb sha256 ccd1db10f90f5295e0265295d0e81817048700669c7e2ea226d46bbe9206d09a (exporta_scaffold6.py, mismo builder verificado que v5; geometria re-verificada: phi -57.8, psi -47.0, omega +-180, RMSD 0.343 A vs 6X18 TM3; quiralidad: clase de CONSTRUCCION L, no observable sin CB en el PDB)
- complejo_v6.pdb sha256 5e8c41ecef3b1ede7d2337d8f659df856db235b1975b92ebe029a1514ceb130c
- MPNN v6: repo 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57, pesos v_48_020.pt sha256 c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd, torch 2.14.1+cpu, python 3.10.12; comandos identicos a provenance.md (parse+assign+run seed 42 temp 0.1, 16 seqs)
- mpnn6/out/seqs/scaffold_v6.fa sha256 0b520cc8baa80b3d03ebe26aa31c4b93e7d73e9b06d5f76d6277de901c57da9a
- Cribado candidatas: pack_eval_v6.py (mismo criterio que pack_eval_ctx.py sobre complejo_v6); salida completa en pack_v6_results.txt
