#!/usr/bin/env python3
"""
crop_review.py -- scientific review of the GLP1R campaign crop, CPU, stdlib-only.
Runs the SAME crop parameters as run_glp1r_colabdesign.py against target/6X18_R.pdb
and reports what could invalidate designs BEFORE any GPU time is spent:

  1. fragments kept + their UniProt ranges
  2. fragment boundaries vs UniProt P43220 topology (a boundary mid-TM-helix means
     an artificially truncated helix -> hallucinated binders may target cut surfaces)
  3. disulfide integrity: Cys pairs (SG-SG <= 2.6 A) with one partner cropped out
  4. artificial exposure: kept residues that lost many contacts (cut neighbours),
     i.e. surfaces that do not exist in the intact receptor
  5. scale check: hotspot span (ECD<->TM) vs GLP-1's own His7->Arg36 CA span in 6X18,
     to size whether a 32-aa peptide can plausibly engage both sites at once

UniProt topology source: rest.uniprot.org/uniprotkb/P43220.json, fetched 24 Oct 2026.
"""
import os, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, '..'))
os.environ.setdefault('GLP1R_WORK', WORK)
R_PDB  = os.path.join(WORK, 'target', '6X18_R.pdb')
RP_PDB = os.path.join(WORK, 'target', '6X18_RP.pdb')
ECD = (24,139)
# Crop logic and config imported from the campaign driver (single source of truth)
sys.path.insert(0, HERE)
from run_glp1r_colabdesign import build_crop, TM_HELICES

def d2(p, q): return (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2
def dist(p, q): return math.sqrt(d2(p, q))

def parse(path):
    atoms = []
    for l in open(path):
        if l.startswith('ATOM'):
            atoms.append(dict(res=l[17:20].strip(), chain=l[21], seq=int(l[22:26]),
                              atom=l[12:16].strip(),
                              p=(float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return atoms

def main():
    atoms = parse(R_PDB)
    n_kept, kept_frag = build_crop()   # writes glp1r_crop.pdb, logs fragments
    kept = {r for f in kept_frag for r in f}
    ca = [(a['chain'], a['seq'], a['p']) for a in atoms if a['atom'] == 'CA']

    print('== 1. fragments (from driver build_crop) ==')
    for i, f in enumerate(kept_frag, 1):
        print(f'  frag{i}: UniProt {f[0]}-{f[-1]} ({len(f)} res)')

    print('== 2. boundaries vs topology ==')
    for f in kept_frag:
        for b in (f[0], f[-1]):
            in_tm = [f'TM{n+1}' for n, (a, z) in enumerate(TM_HELICES) if a + 2 < b < z - 2]
            loc = in_tm[0] if in_tm else ('ECD' if ECD[0] <= b <= ECD[1] else 'loop/other')
            flag = '  <-- MID-HELIX CUT' if in_tm else ''
            print(f'  boundary {b}: {loc}{flag}')

    print('== 3. disulfides ==')
    sg = {a['seq']: a['p'] for a in atoms if a['atom'] == 'SG' and a['res'] == 'CYS'}
    ks = sorted(sg); pairs = []
    for i, r1 in enumerate(ks):
        for r2 in ks[i+1:]:
            if dist(sg[r1], sg[r2]) <= 2.6: pairs.append((r1, r2))
    for r1, r2 in pairs:
        s = 'both kept' if r1 in kept and r2 in kept else ('SEVERED' if (r1 in kept) != (r2 in kept) else 'both cut')
        print(f'  Cys{r1}-Cys{r2}: {s}')

    print('== 4. artificial exposure (kept CA with >=8 cut heavy-atom neighbours <6 A) ==')
    heavy = [a for a in atoms if not a['atom'].startswith('H')]
    for a in heavy:
        if a['seq'] not in kept or a['atom'] != 'CA': continue
        lost = sum(1 for b in heavy if b['seq'] not in kept and a['seq'] != b['seq']
                   and d2(a['p'], b['p']) < 36.0)
        if lost >= 8:
            print(f'  {a["seq"]} {a["res"]}: lost {lost} contacts to cropped-out residues')

    print('== 5. span check ==')
    from run_glp1r_colabdesign import HOT, BINDER_LEN
    hot = {int(h[1:]) for h in HOT.split(',')}
    ecd = [p for c, r, p in ca if r in hot and r <= 139]
    tm  = [p for c, r, p in ca if r in hot and r >= 140]
    d = max(dist(a, b) for a in ecd for b in tm)
    print(f'  max ECD-hotspot <-> TM-hotspot CA distance: {d:.1f} A')
    if os.path.exists(RP_PDB):
        rp = parse(RP_PDB)
        pep = {a['seq']: a['p'] for a in rp if a['chain'] == 'P' and a['atom'] == 'CA'}
        if 7 in pep and 36 in pep:
            span = dist(pep[7], pep[36])
            print(f'  GLP-1 His7->Arg36 CA span (30 aa): {span:.1f} A -> ~{span/29:.1f} A/aa in this geometry')
            print(f'  configured BINDER_LEN={BINDER_LEN}: reaches ~{span/29*(BINDER_LEN-1):.0f} A in the same extended geometry')

if __name__ == '__main__':
    main()
