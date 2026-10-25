#!/usr/bin/env python3
"""
run_glp1r_colabdesign.py -- Challenge 3 (GLP1R agonist) de novo design campaign.
Adapted from challenge2/scripts/run4.py (same ColabDesign binder protocol, same crop
logic, same metrics), aimed at GLP1R using the prepared 6X18 chain R.

HONEST LIMITS (read docs/estrategia.md):
- AF2/ColabDesign optimizes BINDING CONFIDENCE (i_pae, plddt). It says nothing about
  receptor ACTIVATION. Any agonism claim for these designs is a hypothesis to be
  tested in the wet lab, not an in silico result.
- The official de novo / novelty filter is the organizers'. Our local check
  (validate_submission.py) only flags exact >=8-aa windows vs GLP-1/exendin-4.
- save_pdb / metrics / sequence are all read with get_best=True (ch2 self-audit found
  a stale-export defect when mixing best/last; keep all reads best-consistent).

Layout on the DGX (created by setup_dgx_glp1r.sh):
  ~/glp1r/            work root (override with env GLP1R_WORK)
  ~/glp1r/scripts/    this file + prepare_target.py
  ~/glp1r/target/     6X18.cif, 6X18_R.pdb (from prepare_target.py)
  ~/glp1r/docs/       epitope_contacts.tsv

Modes:
  --crop-only   build glp1r_crop.pdb and print fragment statistics, then exit.
                Needs no GPU, no jax, no AF2 params download. Run this FIRST and
                review the fragment report (see docs/estrategia.md checklist).
  (default)     full campaign on GPU (MIG 16 GiB), nohup-friendly.
"""
import os, subprocess, sys, time

# ---- config -----------------------------------------------------------------
WORK       = os.environ.get('GLP1R_WORK', os.path.expanduser('~/glp1r'))
TARGET_PDB = os.path.join(WORK, 'target', '6X18_R.pdb')   # chain R, auth==UniProt numbering
OUT_PREFIX = os.environ.get('G3_PREFIX', 'g3')   # batch tag; g3b for the continuation run
BINDER_LEN = 36          # microbinder class (<40 aa). Sizing evidence (crop_review.py on 6X18):
                         # GLP-1 spans 40.6 A over 30 aa (helical, ~1.4 A/aa); the most extreme
                         # ECD<->TM hotspot pair is 51.7 A. 32 aa reaches ~43 A in the same
                         # geometry - marginal for dual engagement; 36 aa buys margin (~52 A)
                         # while staying in the microbinder class.
CROP_RADIUS = 18.0       # Angstrom from hotspot CAs (20.0 in ch2; 18 keeps the 7TM crop smaller)
BUDGET_MIN = int(os.environ.get('G3_BUDGET_MIN', '80'))
MAX_TRAJS = int(os.environ.get('G3_MAX_TRAJS', '0'))   # 0 = sin tope; 1 = una trayectoria por proceso (gate de revision entre trayectorias)   # wall-clock minutes; checked BETWEEN trajectories
# Hotspots: strongest computed GLP-1 contacts (min dist <= ~3.4 A) spanning the ECD
# groove AND the TM pocket, from docs/epitope_contacts.tsv (auth == UniProt numbering).
# HYPOTHESIS, not a guarantee: dual-site engagement is the agonist geometry GLP-1
# uses in 6X18; a hallucinated 32-mer may satisfy the loss binding the ECD groove
# alone. Post-run pose analysis decides what each trajectory actually engaged.
HOT_ECD = "R32,R39,R68,R91,R121,R138,R141"
HOT_TM  = "R152,R190,R197,R205,R214,R234,R298,R299,R300,R310,R372,R380"
HOT = HOT_ECD + "," + HOT_TM
# UniProt P43220 TM helices (rest.uniprot.org, 24 Oct 2026). Used to avoid
# mid-helix crop cuts: a boundary falling inside a helix is extended to the
# helix end so no truncated helix presents an artificial surface.
TM_HELICES = [(140,164),(176,201),(228,251),(266,290),(306,328),(349,370),(384,404)]
EXTEND_TO_HELIX_ENDS = True   # cost measured 24 Oct: +26 res -> 315+36 = 351 total
                              # (~0.48x the 506-res system that OOM'd on this MIG)
# -----------------------------------------------------------------------------

def log(*a): print(*a, flush=True)

def _d2(p, q):
    return (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2

def build_crop():
    """Crop target chain R around hotspots. Stdlib-only (no numpy): runs before the
    venv exists on the DGX. Returns (n_kept, fragments)."""
    _L = [l for l in open(TARGET_PDB) if l.startswith('ATOM')]
    _hot = {('R', int(h[1:])) for h in HOT.split(',')}
    _ca = [(l[21], int(l[22:26]), (float(l[30:38]), float(l[38:46]), float(l[46:54])))
           for l in _L if l[12:16] == ' CA ']
    _hp = [p for c, r, p in _ca if (c, r) in _hot]
    if not _hp: sys.exit('no hotspot CAs found - check chain letter/numbering in ' + TARGET_PDB)
    _r2 = CROP_RADIUS ** 2
    _keep = {(c, r) for c, r, p in _ca if min(_d2(p, q) for q in _hp) < _r2}
    _k2 = set(_keep)
    _rs = sorted(r for c, r in _keep if c == 'R')
    for a, b in zip(_rs, _rs[1:]):
        if 1 < b - a <= 6: _k2 |= {('R', x) for x in range(a, b)}
    _rs = sorted(r for c, r in _k2 if c == 'R'); frag = []; cur = []
    for r in _rs:
        if cur and r - cur[-1] > 1: frag.append(cur); cur = []
        cur.append(r)
    if cur: frag.append(cur)
    dropped = [f for f in frag if len(f) < 4]
    for f in dropped: _k2 -= {('R', x) for x in f}
    kept_frag = [f for f in frag if len(f) >= 4]
    if EXTEND_TO_HELIX_ENDS:
        ext = []
        for f in kept_frag:
            a, b = f[0], f[-1]
            for h0, h1 in TM_HELICES:
                if h0 < a < h1: a = h0
                if h0 < b < h1: b = h1
            ext.append(list(range(a, b + 1)))
        # merge fragments that now overlap or touch
        ext.sort(key=lambda f: f[0]); merged = []
        for f in ext:
            if merged and f[0] <= merged[-1][-1] + 1:
                merged[-1] = list(range(merged[-1][0], max(merged[-1][-1], f[-1]) + 1))
            else:
                merged.append(f)
        if merged != kept_frag:
            log('EXTEND mid-helix cuts to helix ends:',
                [f'{f[0]}-{f[-1]}' for f in kept_frag], '->', [f'{f[0]}-{f[-1]}' for f in merged])
        kept_frag = merged
        _k2 = {('R', x) for f in kept_frag for x in f}
    out = os.path.join(WORK, 'glp1r_crop.pdb')
    open(out, 'w').write(''.join(l for l in _L if (l[21], int(l[22:26])) in _k2))
    log('CROP', out, '| residues', len(_k2), '| fragments', len(kept_frag),
        '| dropped(<4)', [f'{f[0]}-{f[-1]}' for f in dropped])
    for i, f in enumerate(kept_frag, 1):
        log(f'  frag{i}: UniProt {f[0]}-{f[-1]} ({len(f)} res)')
    return len(_k2), kept_frag

def main():
    t0 = time.time()
    def sh(c):
        r = subprocess.run(c, shell=True, capture_output=True, text=True)
        print(c[:80], '->', r.returncode, (r.stdout + r.stderr)[-300:]); return r
    os.chdir(WORK)
    crop_only = '--crop-only' in sys.argv
    n, frags = build_crop()
    if crop_only:
        log('CROP-ONLY done; review fragments against docs/estrategia.md before any GPU run')
        return

    sh("mkdir -p params && ls params | head -3")
    if not os.path.exists('params/params_model_1_ptm.npz'):
        sh("curl -fsSL https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar | tar x -C params")

    import jax; log('jax', jax.__version__, jax.devices())
    import jax.numpy as jnp
    _c = jnp.clip
    def _clip(a, *args, a_min=None, a_max=None, **k):   # same compat patch as ch2
        if a_min is not None: k['min'] = a_min
        if a_max is not None: k['max'] = a_max
        return _c(a, *args, **k)
    jnp.clip = _clip
    import numpy as np, json

    DRY = (jax.devices()[0].platform == 'cpu')
    BUDGET = 60 if DRY else BUDGET_MIN * 60
    log('DRY' if DRY else 'GPU', 'budget_s', BUDGET)

    from colabdesign import mk_afdesign_model, clear_mem
    from colabdesign.shared.utils import copy_dict
    af = mk_afdesign_model(protocol="binder", use_multimer=False, data_dir=".", best_metric="i_pae")

    # Memory patch proven necessary on this MIG (ch2 diagnostics 11 Oct,
    # diagnostics/mig/full456_mig/prep_patch_after.txt): force subbatch_size=1 so
    # attention is sub-batched. Note: even WITH this patch, a 506-res system
    # (456 target + 50 binder) OOM'd on the 16 GiB slice; our crop is 351 total.
    def _prep_model_subbatch(self, **kwargs):
        if not hasattr(self, "_model") or self._cfg != self._model["runner"].config:
            self._cfg.model.global_config.subbatch_size = 1
            self._model = self._get_model(self._cfg)
            if sum(self._lengths) > 384:
                self._cfg.model.global_config.subbatch_size = 1
                self._model["fn"] = self._get_model(self._cfg)["fn"]
        self._opt = copy_dict(self.opt)
        self.restart(**kwargs)
    af._prep_model = _prep_model_subbatch.__get__(af, type(af))
    af.prep_inputs(pdb_filename="glp1r_crop.pdb", chain="R", hotspot=HOT,
                   binder_len=BINDER_LEN, rm_target_seq=False)
    # NOTE: hotspot chain letters refer to the chain letters IN glp1r_crop.pdb (R),
    # matching how ch2 passed hotspot letters from its cropped PDB.

    rows = []; k = 0
    while (time.time() - t0 < BUDGET) and (MAX_TRAJS <= 0 or k < MAX_TRAJS):   # v5: tope de trayectorias para review-gate
        k += 1
        af.restart(seed=k + 300)
        af.design_3stage(30 if DRY else 100, 30 if DRY else 100, 5 if DRY else 10)
        chosen = af._tmp['best'].get('aux', af.aux); lg = chosen['log']
        seq = af.get_seqs(get_best=True)[0]
        af.save_pdb(f'{OUT_PREFIX}_{k}.pdb', get_best=True)   # best-consistent export
        rec = dict(traj=k, seq=seq, len=len(seq), plddt=round(lg['plddt'], 3),
                   i_pae=round(lg['i_pae'], 3), i_ptm=round(lg['i_ptm'], 3))
        rows.append(rec)
        log('TRAJ', rec, 'elapsed', int(time.time() - t0))
        json.dump(rows, open(f'{OUT_PREFIX}_results.json', 'w'), indent=1)   # incremental: survives a mid-batch crash
        if k == 1:
            sh("nvidia-smi --query-gpu=memory.used,memory.total --format=csv")  # real VRAM footprint of one trajectory
        if DRY: break
        # NO clear_mem(): it frees arrays still referenced by the model's grad_fn and
        # kills the NEXT trajectory's first step (RuntimeError: Array has been deleted).
        # Measured footprint is stable (12,798 MiB of 16 GiB across trajectory 1).

    rows.sort(key=lambda r: r['i_pae'])
    json.dump(rows, open(f'{OUT_PREFIX}_results.json', 'w'), indent=1)
    with open(f'{OUT_PREFIX}_ranked.csv', 'w') as fh:
        fh.write('rank,traj,seq,plddt,i_pae,i_ptm\n')
        for i, r in enumerate(rows, 1):
            fh.write(f"{i},{r['traj']},{r['seq']},{r['plddt']},{r['i_pae']},{r['i_ptm']}\n")
    log('DONE trajectories', len(rows), '->', f'{OUT_PREFIX}_ranked.csv')

if __name__ == '__main__':
    main()
