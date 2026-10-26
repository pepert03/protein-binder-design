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
OUT_PREFIX = os.environ.get('OUT_PREFIX', 'g3')
BINDER_LEN = 32          # microbinder class (<40 aa); GLP-1 itself is 30 aa
CROP_RADIUS = 18.0       # Angstrom from hotspot CAs (20.0 in ch2; 18 keeps the 7TM crop smaller)
BUDGET_MIN  = int(os.environ.get('BUDGET_MIN', '80'))         # wall-clock minutes of design (8h session cap -> several runs per session)
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
EXTEND_TO_HELIX_ENDS = True
# -----------------------------------------------------------------------------

def log(*a): print(*a, flush=True)

def build_crop():
    """Crop target chain R around hotspots. Returns (n_kept, fragments)."""
    import numpy as np
    _L = [l for l in open(TARGET_PDB) if l.startswith('ATOM')]
    _hot = {('R', int(h[1:])) for h in HOT.split(',')}
    _ca = [(l[21], int(l[22:26]), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
           for l in _L if l[12:16] == ' CA ']
    _hp = [p for c, r, p in _ca if (c, r) in _hot]
    if not _hp: sys.exit('no hotspot CAs found - check chain letter/numbering in ' + TARGET_PDB)
    _keep = {(c, r) for c, r, p in _ca if min(np.linalg.norm(p - q) for q in _hp) < CROP_RADIUS}
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
    af = mk_afdesign_model(protocol="binder", use_multimer=False, data_dir=".", best_metric="i_pae")
    af.prep_inputs(pdb_filename="glp1r_crop.pdb", chain="R", hotspot=HOT,
                   binder_len=BINDER_LEN, rm_target_seq=False)
    # NOTE: hotspot chain letters refer to the chain letters IN glp1r_crop.pdb (R),
    # matching how ch2 passed hotspot letters from its cropped PDB.

    import gc as _gc
    _incr = f'{OUT_PREFIX}_results_incr.json'
    if os.path.exists(_incr):
        rows = json.load(open(_incr))
    else:
        rows = []
    k = len(rows); fails = 0
    def _new_model():
        m = mk_afdesign_model(protocol="binder", use_multimer=False, data_dir=".", best_metric="i_pae")
        m.prep_inputs(pdb_filename="glp1r_crop.pdb", chain="R", hotspot=HOT,
                      binder_len=BINDER_LEN, rm_target_seq=False)
        return m
    while time.time() - t0 < BUDGET:
        k += 1
        try:
            af.restart(seed=k + int(os.environ.get('SEED0', '300')))
            af.design_3stage(30 if DRY else 100, 30 if DRY else 100, 5 if DRY else 10)
        except Exception as e:
            log('TRAJ_FAIL', k, type(e).__name__, str(e)[:200])
            fails += 1
            k -= 1
            if fails >= 3:
                log('TOO_MANY_FAILS'); break
            try:
                af = _new_model()
            except Exception as e2:
                log('REBUILD_FAIL', type(e2).__name__, str(e2)[:200]); break
            continue
        fails = 0
        chosen = af._tmp['best'].get('aux', af.aux); lg = chosen['log']
        seq = af.get_seqs(get_best=True)[0]
        af.save_pdb(f'{OUT_PREFIX}_{k}.pdb', get_best=True)   # best-consistent export
        rec = dict(traj=k, seq=seq, len=len(seq), plddt=round(lg['plddt'], 3),
                   i_pae=round(lg['i_pae'], 3), i_ptm=round(lg['i_ptm'], 3))
        rows.append(rec)
        log('TRAJ', rec, 'elapsed', int(time.time() - t0))
        json.dump(rows, open(f'{OUT_PREFIX}_results_incr.json', 'w'), indent=1)
        if DRY: break
        _gc.collect()

    rows.sort(key=lambda r: r['i_pae'])
    json.dump(rows, open(f'{OUT_PREFIX}_results.json', 'w'), indent=1)
    with open(f'{OUT_PREFIX}_ranked.csv', 'w') as fh:
        fh.write('rank,traj,seq,plddt,i_pae,i_ptm\n')
        for i, r in enumerate(rows, 1):
            fh.write(f"{i},{r['traj']},{r['seq']},{r['plddt']},{r['i_pae']},{r['i_ptm']}\n")
    log('DONE trajectories', len(rows), '->', f'{OUT_PREFIX}_ranked.csv')

if __name__ == '__main__':
    main()