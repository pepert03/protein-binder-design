# Provenance ProteinMPNN CPU sobre scaffold_v5 (2026-10-25)

## Entrada
- scaffold_v5.pdb  sha256 7931f785b46a38c7794d4d746aafd33c0b6dc4bb4f7cf82393c9f59ed61087da
  (publicada por verificador en repo, commit 6b33d50)
- parsed.jsonl     sha256 59e2008973280f146fc02f8b5e73934bd37eefc3837ba0af038006c248fbf0e0
- assigned.jsonl   sha256 04c675dcae3d11256a80bc452d107678484697c21f76b21646a8c7d657f83969
- omit.jsonl       sha256 a8d73f0fd08cd76fa5cfd7a65af46d0c9805a430410fbcdd88d52a10a56d8269

## Codigo y pesos
- Repo: github.com/dauparas/ProteinMPNN, commit 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57 (arbol limpio)
- Pesos: vanilla_model_weights/v_48_020.pt sha256 c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd

## Comandos exactos
cp scaffold_v5.pdb /tmp/pose/mpnn/in/
python3 helper_scripts/parse_multiple_chains.py --input_path=/tmp/pose/mpnn/in --output_path=/tmp/pose/mpnn/parsed.jsonl
python3 helper_scripts/assign_fixed_chains.py --input_path=/tmp/pose/mpnn/parsed.jsonl --output_path=/tmp/pose/mpnn/assigned.jsonl --chain_list "B"
# corrida LIBRE:
python3 protein_mpnn_run.py --jsonl_path /tmp/pose/mpnn/parsed.jsonl --chain_id_jsonl /tmp/pose/mpnn/assigned.jsonl --out_folder /tmp/pose/mpnn/out --num_seq_per_target 16 --sampling_temp "0.1" --seed 42 --batch_size 1 --model_name "v_48_020"
# corrida RESTRINGIDA (B3,B6 solo AGSTV):
python3 protein_mpnn_run.py --jsonl_path /tmp/pose/mpnn/parsed.jsonl --chain_id_jsonl /tmp/pose/mpnn/assigned.jsonl --omit_AA_jsonl /tmp/pose/mpnn/omit.jsonl --out_folder /tmp/pose/mpnn/out_restr --num_seq_per_target 16 --sampling_temp "0.1" --seed 42 --batch_size 1 --model_name "v_48_020"

## Entorno
python 3.10.12 | torch 2.14.1+cpu | numpy 2.2.6 | CUDA no | 2 nucleos | ~46 s por corrida

## Salidas y determinismo (re-ejecucion 14:31-14:32 CET, logs adjuntos)
- LIBRE       sha256 62675e45f8bf6aa65d82921a791a20e77e13e6b1bd0f0593b13f406f5ee8b546  (re-run byte-identico)
- RESTRINGIDA sha256 6299aefc6dd98ac8571a26330890a17d85e367b3e53893e2b08b1e2377113de8  (re-run byte-identico)
- sha1 / git-blob / md5 de ambas en cb_map reporte: NINGUNO coincide con los hashes del verificador
  (453d904d libre / 5893bd87 restringida) -> discrepancia byte-nivel a reconciliar.

## Mapa CB v5 DEFINITIVO (formula canonica validada)
- Dos construcciones canonicas independientes: coeficientes ProteinMPNN (protein_mpnn_utils.py l.968-970)
  y tetraedrica exacta (angulo 109.47, quiralidad decidida por 301 CB reales de 6X18: signo +, med 0.064 A).
- Ambas validadas contra 301 CB reales: med error 0.056 / 0.064 A, max 0.52 / 0.55 A.
- Divergencia mutua por posicion: ~0.046 A. Mapa UNICO:
    B2 = 1.17 A (CD1 A356)   << conflicto
    B6 = 1.71 A (CD1 A113)   << conflicto
    B31 = 1.96 A (CA A90)    << conflicto
    B28 = 2.68 A LIMPIO | B32 = 3.52 LIMPIO | B3 = 3.63 LIMPIO | B4 = 4.30 LIMPIO
- Retractado: formula en plano anterior (med error 1.34 A) y sus listas (B6/B28, y B3/B4 que venia de v3).
- No reproducible: tetraedrico verificador B28 1.86 / B32 1.97 / B6 3.23. Nota: su B32=1.97 == mi B31=1.97
  (sospecha de desfase de indice en su pipeline). Su bb_check (B6 1.83/B28 1.45) = mi formula simple retractada.
