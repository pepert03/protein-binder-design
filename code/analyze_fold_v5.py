"""Geometry analysis of Boltz-2 v5 folds vs template backbones.

Checks per design:
  - ipTM from confidence json (already known)
  - min heavy-atom distance from each designed protonation center to EGFR H433 (model His406)
  - binder-target interface contact residues (< 4.5 A)
  - distance from binder to model Asn325 (glycan proxy for UniProt N352)
  - pose shift vs template: RMSD of binder C-alpha after superposing template on fold via EGFR
"""
import gemmi, json, math, os, sys

BASE = os.path.expanduser('~/workspace/proteinbase-egfr')
PREDS = os.path.join(BASE, 'hidden_runs/fold_v5_out/preds')
CAND = os.path.join(BASE, 'data/candidates')

DESIGNS = {
    'r03_l05': ('rank03_cC5_01', [17, 18, 36]),
    'r03_l07': ('rank03_cC5_01', [17, 18, 36]),
    'r01_l07': ('rank01_cC5_03', [3, 5, 16]),
    'r02_l07': ('rank02_cC5_13', [4, 35, 39]),
}

def chains(st):
    out = []
    for ch in st[0]:
        n = sum(1 for r in ch if r.name != 'HOH')
        out.append((ch.name, n))
    return out

def atoms_of(x):
    # x: gemmi.Residue, Chain, or Model-like iterable
    if isinstance(x, gemmi.Residue):
        return list(x)
    return [a for r in x for a in r]

def min_dist(res_a, res_b):
    dmin = 1e9
    for a in atoms_of(res_a):
        for b in atoms_of(res_b):
            d = a.pos.dist(b.pos)
            if d < dmin:
                dmin = d
    return dmin

def ca_atoms(ch):
    return [r['CA'][0].pos for r in ch if 'CA' in r]

def kabsch_rmsd(pa, pb):
    # pa, pb: lists of gemmi.Position, same length; returns RMSD after optimal superposition of pb onto pa
    import numpy as np
    A = np.array([[p.x, p.y, p.z] for p in pa])
    B = np.array([[p.x, p.y, p.z] for p in pb])
    ca = A.mean(0); cb = B.mean(0)
    A -= ca; B -= cb
    H = B.T @ A
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = ca - R @ cb
    B2 = (R @ B.T).T + t
    return float(np.sqrt(((A - B2) ** 2).sum(1).mean()))

def analyze(name, backbone, centers):
    pred_cif = os.path.join(PREDS, f'boltz_results_{name}', 'predictions', name, f'{name}_model_0.cif')
    tmpl_cif = os.path.join(CAND, f'{backbone}.cif')
    st = gemmi.read_structure(pred_cif)
    chinfo = chains(st)
    # binder = short chain, target = long chain
    chinfo.sort(key=lambda x: x[1])
    bname, tname = chinfo[0][0], chinfo[-1][0]
    binder = st[0][bname]; target = st[0][tname]
    bres = [r for r in binder]
    tres = [r for r in target]
    print(f'=== {name} (backbone {backbone}) chains binder={bname}({len(bres)}) target={tname}({len(tres)})')
    # H433 = model His406
    h433 = target['406'][0]
    print(f'    target res 406 is {h433.name}')
    asn325 = None
    for r in tres:
        if r.seqid.num == 325 and r.name == 'ASN':
            asn325 = r; break
    for c in centers:
        # centers are 1-indexed positions into binder sequence
        cres = bres[c - 1]
        d = min_dist(cres, h433)
        print(f'    center {c} ({cres.name}): min dist to H433 = {d:.2f} A')
    # interface contacts
    dmin_if = []
    for r in bres:
        d = min_dist(r, target)
        if d < 4.5:
            dmin_if.append((r.seqid.num, r.name, d))
    print(f'    binder interface residues (<4.5A, n={len(dmin_if)}): ' +
          ', '.join(f'{n}{nm}({d:.1f})' for n, nm, d in dmin_if))
    if asn325 is not None:
        dn = min_dist(binder, asn325)
        print(f'    min dist binder -> Asn325 (glycan proxy) = {dn:.2f} A')
    # pose shift vs template
    ts = gemmi.read_structure(tmpl_cif)
    tch = None
    for ch in ts[0]:
        if sum(1 for r in ch) > 100:
            tch = ch
    # find binder chain in template = the shorter one
    bch = min((ch for ch in ts[0]), key=lambda c: sum(1 for r in c))
    tf_ca = ca_atoms(tch); tmpl_b = ca_atoms(bch)
    fold_t_ca = ca_atoms(target); fold_b = ca_atoms(binder)
    n = min(len(tf_ca), len(fold_t_ca), len(tmpl_b), len(fold_b))
    # superpose template on fold via target, then binder RMSD
    import numpy as np
    A = np.array([[p.x, p.y, p.z] for p in fold_t_ca[:n]])
    B = np.array([[p.x, p.y, p.z] for p in tf_ca[:n]])
    ca = A.mean(0); cb = B.mean(0)
    A0 = A - ca; B0 = B - cb
    H = B0.T @ A0
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T; t = ca - R @ cb
    nb = min(len(tmpl_b), len(fold_b))
    Bb = (R @ np.array([[p.x, p.y, p.z] for p in tmpl_b[:nb]]).T).T + t
    Ab = np.array([[p.x, p.y, p.z] for p in fold_b[:nb]])
    rmsd = float(np.sqrt(((Ab - Bb) ** 2).sum(1).mean()))
    print(f'    binder pose RMSD vs template (target-superposed) = {rmsd:.2f} A')
    print()

for name, (bb, centers) in DESIGNS.items():
    analyze(name, bb, centers)
