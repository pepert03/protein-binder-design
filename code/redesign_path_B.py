#!/usr/bin/env python3
"""Path B: pH-switch redesign of the best Path-A backbones with Proton-PottsMPNN.

Adapts the packaged inference/design_ph.py pipeline to our recovered EGFR
complexes (mmCIF, binder = chain A). For each backbone it runs a
combined_lambda sweep (block descent, production optimiser) and writes the
designs + Pareto front to data/designs/pathB/<backbone>/.

Protonated centres are restricted to the interface (placement_region=
["interface"]) because the pH switch must act at the binding site.

Usage:
    python code/redesign_path_B.py --pilot        # rank03 only, 1 design, timing check
    python code/redesign_path_B.py                 # full sweep on 3 backbones

Requires HBPLUS_PATH set to the compiled hbplus binary.
"""
import argparse
import json
import os
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.pop("DEBUG", None)
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

PKG = Path.home() / "workspace" / "ppmnn" / "ProtonPottsMPNN"
VENV_SITE = PKG / ".venv" / "lib" / "python3.12" / "site-packages"
sys.path.insert(0, str(PKG / "foundry" / "models" / "mpnn" / "src"))
sys.path.insert(0, str(PKG / "inference"))

from atomworks.io import parse  # noqa: E402

from mpnn.inference_engines.potts_mpnn_ph import (  # noqa: E402
    PottsMPNNPHEngine,
    PHDesignCriteria,
)

WORK = Path.home() / "workspace" / "proteinbase-egfr"
CAND = WORK / "data" / "candidates"
OUT = WORK / "data" / "designs" / "pathB"
CKPT = PKG / "checkpoints" / "potts_v6_afdb_edge_his0.3_acid0.06" / "epoch-0125.ckpt"

# Backbones chosen from Path A (placement_scan.json) + H433 geometry:
# rank03: best H433 geometry (Asp37 5.7 A) + moderate binder HIS-P (-1.70)
# cC5_00: strongest binder HIS-P interface signal (-3.26); glycan-clash risk noted
# rank01: second-best binder HIS-P (-2.56)
BACKBONES = {
    "rank03_cC5_01": "rank03_cC5_01.cif",
    "rank01_cC5_03": "rank01_cC5_03.cif",
    "rank02_cC5_13": "rank02_cC5_13.cif",
    # cC5_00 excluded from the first sweep: strongest binder HIS-P signal
    # (-3.26) but documented strong clash with the first EGFR-attached sugar
    # in two alignments (docs/STATUS.md); Proton-PottsMPNN is fixed-backbone
    # so redesign cannot move the binder away from the glycan. Revisit only
    # if the sweep on the three above fails to produce a Pareto front.
}
BINDER_CHAIN = "A"
LAMBDAS = [0.1, 0.3, 0.5, 0.7]


def base_criteria(trajectory: bool = False) -> PHDesignCriteria:
    return PHDesignCriteria(
        method="block_descent",
        backend="potts",
        temperature=0.05,
        samples_per_site=1,
        block_size=3,
        seed_source="native",
        center_types=["HIS-P", "ASP-P", "GLU-P"],
        dep_map={"HIS-P": ["HIS-S"], "ASP-P": ["ASP-D"], "GLU-P": ["GLU-D"]},
        forbidden_tokens=["HIS-A", "ASP-A", "GLU-A", "UNK"],
        placement_by="scan_potts",
        placement_region=["interface"],  # pH switch must act at the binding site
        repetitive_window_parents=["ARG", "LYS", "HIS", "ASP", "GLU"],
        repetitive_window_radius=2,
        repetitive_window_weight=1.0,
        neighbour_k=16,
        max_mutations=20,
        record_trajectory=trajectory,
    )


def run_backbone(name: str, lambdas, trajectory: bool = False):
    outdir = OUT / name
    outdir.mkdir(parents=True, exist_ok=True)
    cif = CAND / BACKBONES[name]
    print(f"== {name} ({cif.name}) ==", flush=True)

    engine = PottsMPNNPHEngine(
        checkpoint_path=str(CKPT), extended_vocab="v6",
        out_directory=None, write_fasta=False, write_structures=False,
    )
    atom_array = parse(str(cif))["asym_unit"][0]  # atomworks AtomArray (has charge)
    ctx = engine._build_context(atom_array, BINDER_CHAIN)
    n_free = len(ctx.chA_free_idx)
    print(f"  {n_free} designable binder positions", flush=True)

    crit = base_criteria(trajectory=trajectory)
    sweep_criteria = [
        replace(crit, combined_lambda=float(lam)) for lam in lambdas
    ]
    n_jobs = min(len(sweep_criteria), os.cpu_count() or 1)
    design_set = engine.run_ph_redesign(
        atom_array=atom_array, binder_chain=BINDER_CHAIN,
        criteria_list=sweep_criteria, seed=0, n_jobs=n_jobs,
    )
    print(f"  produced {len(design_set)} designs", flush=True)

    rows = []
    for d in design_set:
        rows.append({
            "design_id": d.design_id(),
            "backbone": name,
            "combined_lambda": d.combined_lambda,
            "potts_energy": d.final_potts_energy,
            "selective_energy": d.selective_energy,
            "centers": ";".join(
                f"{r}:{t}" for r, t in
                zip(d.center_res_ids or [], d.center_protonation_types or [])),
            "canonical_sequence": d.canonical_sequence,
            "extended_tokens": " ".join(d.extended_tokens),
        })
    df = pd.DataFrame(rows).sort_values("combined_lambda").reset_index(drop=True)
    df.to_csv(outdir / "sweep_designs.tsv", sep="\t", index=False)

    # Pareto front (minimise both energies)
    P = df[["potts_energy", "selective_energy"]].to_numpy()
    dominated = np.zeros(len(P), bool)
    for i in range(len(P)):
        for j in range(len(P)):
            if i != j and (P[j] <= P[i]).all() and (P[j] < P[i]).any():
                dominated[i] = True
                break
    df["pareto_optimal"] = ~dominated
    df.to_csv(outdir / "sweep_designs.tsv", sep="\t", index=False)

    with open(outdir / "designs.json", "w") as fh:
        json.dump([d.to_metadata() for d in design_set], fh, indent=2, default=float)

    with open(outdir / "designs.fasta", "w") as fh:
        for d in design_set:
            fh.write(f">{d.design_id()} H={d.final_potts_energy:.3f} "
                     f"sel={d.selective_energy}\n{d.canonical_sequence}\n")

    front = df[~dominated]
    print(f"  pareto-optimal: {len(front)}/{len(df)}", flush=True)
    for _, r in front.iterrows():
        print(f"    lam={r.combined_lambda} H={r.potts_energy:.2f} "
              f"sel={r.selective_energy} centers={r.centers}", flush=True)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true",
                    help="rank03 only, single design at lambda=0.3 with trajectory")
    args = ap.parse_args()

    if args.pilot:
        df = run_backbone("rank03_cC5_01", [0.3], trajectory=True)
    else:
        for name in BACKBONES:
            run_backbone(name, LAMBDAS)
    print("done")


if __name__ == "__main__":
    main()
