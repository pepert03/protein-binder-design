"""Path A: placement-scan rescoring of the 9 recovered candidate complexes.

For each candidate CIF: featurize with Proton-PottsMPNN (v6 vocab), compute the
selective placement score of each protonation state at every free binder
position, and report the top placements, especially interface positions near
target His406 (= UniProt His433, the pH sensor).
"""
import os, sys, glob, json
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

import numpy as np
import torch
from atomworks.io import parse as aw_parse
from biotite.structure.io.pdbx import CIFFile, get_structure

PKG = os.path.expanduser("~/workspace/ppmnn/ProtonPottsMPNN")
sys.path.insert(0, os.path.join(PKG, "foundry", "models", "mpnn", "src"))
sys.path.insert(0, os.path.join(PKG, "inference"))
from mpnn.inference_engines.potts_mpnn_ph import PottsMPNNPHEngine
from design_placement_scan import classify_binder, STATES, STATE_ORDER

CKPT = os.path.join(PKG, "checkpoints", "potts_v6_afdb_edge_his0.3_acid0.06", "epoch-0125.ckpt")
CAND = os.path.expanduser("~/workspace/proteinbase-egfr/data/candidates")

engine = PottsMPNNPHEngine(checkpoint_path=CKPT, extended_vocab="v6",
                           out_directory=None, write_fasta=False, write_structures=False)
t2i = None
results = {}

for fn in sorted(glob.glob(os.path.join(CAND, "rank*.cif"))) + [os.path.join(CAND, "cC5_00.cif")]:
    name = os.path.basename(fn).replace(".cif", "")
    aa = aw_parse(fn)["asym_unit"][0]  # atomworks AtomArray (has charge annotation)
    aa_bio = get_structure(CIFFile.read(fn), model=1)  # clean coords for SASA classification
    binder_chain = "A"
    cls = classify_binder(aa_bio, binder_chain)
    ctx = engine._build_context(aa, binder_chain)
    if t2i is None:
        t2i = ctx.encoding.token_to_idx
    region = ctx.chA_free_idx
    res_ids = np.array([int(ctx.token_aa.res_id[int(p)]) for p in region.tolist()])
    pos_cls = np.array([cls.get(int(r), "surface") for r in res_ids])
    seq = ctx.S_native.clone()
    with torch.no_grad():
        ef = ctx.field_potts(seq)
    out = {}
    for stt, deps in STATES.items():
        prot = ef[region, t2i[stt]]
        dep = torch.stack([ef[region, t2i[d]] for d in deps], 0).min(0).values
        scores = (prot - dep).cpu().numpy()
        k = int(np.argmin(scores))
        out[stt] = {"res": int(res_ids[k]), "class": str(pos_cls[k]), "score": float(scores[k]),
                    "all": [(int(r), str(c), float(s)) for r, c, s in zip(res_ids, pos_cls, scores)]}
    # interface positions sorted by HIS-P selective score (most preferred first)
    hip = sorted([(r, c, s) for r, c, s in out["HIS-P"]["all"] if c == "interface"], key=lambda x: x[2])[:8]
    asp = sorted([(r, c, s) for r, c, s in out["ASP-P"]["all"] if c == "interface"], key=lambda x: x[2])[:8]
    glu = sorted([(r, c, s) for r, c, s in out["GLU-P"]["all"] if c == "interface"], key=lambda x: x[2])[:8]
    results[name] = {"top": {s: {k2: v for k2, v in out[s].items() if k2 != "all"} for s in STATE_ORDER},
                     "iface_HIS-P": hip, "iface_ASP-P": asp, "iface_GLU-P": glu}
    print(f"== {name} ==")
    for s in STATE_ORDER:
        t = out[s]
        print(f"  {s:6s} top: binder res {t['res']:3d} class={t['class']:9s} score={t['score']:+.3f}")
    print(f"  iface HIS-P top3: {[(r, round(s,2)) for r, c, s in hip[:3]]}")
    print(f"  iface ASP-P top3: {[(r, round(s,2)) for r, c, s in asp[:3]]}")
    print(f"  iface GLU-P top3: {[(r, round(s,2)) for r, c, s in glu[:3]]}", flush=True)

json.dump(results, open(os.path.join(CAND, "placement_scan.json"), "w"), indent=1)
print("saved placement_scan.json")
