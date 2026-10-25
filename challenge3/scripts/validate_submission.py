#!/usr/bin/env python3
"""
validate_submission.py -- self-check for a Challenge 3 submission CSV. CPU, stdlib-only.

Usage: python3 validate_submission.py my_submission.csv

Checks (official rules from the challenge page + competition FAQ, fetched 24 Oct 2026):
  - required columns: name, sequence, molecule_class
  - molecule_class in {protein, nanobody, scfv, fab_kappa, fab_lambda}
  - protein/single-chain: 10-250 aa, standard 20 amino acids only
  - fab_*: VH:VL joined by ':' per FAQ 5
  - unique names, unique sequences (exact duplicates are a hard filter)
  - ranking: CSV order IS the ranking (top row = most preferred); prints the order
  - novelty self-check: flags any design sharing an exact >=8-aa window with the
    natural agonist GLP-1(7-37) / GLP-1(7-36)NH2 or exendin-4. This is a REVIEW
    flag for the human, not a verdict: the official de novo / zero-shot filter is
    run by the organizers (sequence AND structural diversity, CDR diversity for
    antibodies). We cannot run their filter locally.
  - tag warning: flags sequences containing common tag motifs (Twin-Strep, His6,
    FLAG) fused in the design itself; the design will get a C-terminal Twin-Strep
    capture tag in the assay, so paratope must NOT be at the C-terminus (manual check).

This script never sends anything anywhere. Local validation only.
"""
import csv, sys, itertools

STD = set('ACDEFGHIKLMNPQRSTVWY')
CLASSES = {'protein', 'nanobody', 'scfv', 'fab_kappa', 'fab_lambda'}
# Reference agonists to stay away from (literature sequences; GLP-1 from PDB 6X18 chain P).
EXCLUDED = {
    'GLP-1(7-37)':    'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG',
    'GLP-1(7-36)NH2': 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR',
    'exendin-4':      'HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS',
}
TAG_MOTIFS = {'Twin-Strep': 'WSHPQFEK', 'His6': 'HHHHHH', 'FLAG': 'DYKDDDDK'}
K = 8

def kmers(s, k=K):
    return {s[i:i+k] for i in range(len(s)-k+1)}

EXCLUDED_KMERS = {name: kmers(s) for name, s in EXCLUDED.items()}

def fail(msg): print(f'  FAIL  {msg}')
def warn(msg): print(f'  WARN  {msg}')
def ok(msg):   print(f'  ok    {msg}')

def main(path):
    rows = list(csv.DictReader(open(path)))
    if not rows: sys.exit('empty CSV')
    cols = {c.strip().lower() for c in (rows[0].keys() or [])}
    missing = {'name', 'sequence', 'molecule_class'} - cols
    if missing: sys.exit(f'missing required columns: {missing}')
    errors = 0
    names, seqs = set(), set()
    print(f'{len(rows)} designs; CSV row order is the submission ranking:')
    for i, row in enumerate(rows, 1):
        g = {k.strip().lower(): (v or '').strip() for k, v in row.items() if k}
        name, seq, cls = g['name'], g['sequence'].upper(), g['molecule_class'].lower()
        label = f'row {i} [{name}]'
        if not name: fail(f'{label}: empty name'); errors += 1
        if name in names: fail(f'{label}: duplicate name'); errors += 1
        names.add(name)
        if seq in seqs: fail(f'{label}: duplicate sequence (hard filter)'); errors += 1
        seqs.add(seq)
        if cls not in CLASSES:
            fail(f'{label}: molecule_class {cls!r} not in {sorted(CLASSES)}'); errors += 1
        core = seq.replace(':', '')
        bad = set(core) - STD
        if bad: fail(f'{label}: non-standard residues {sorted(bad)} (Aib/modified AAs cannot be submitted as plain sequence)'); errors += 1
        if cls == 'protein' and not (10 <= len(core) <= 250):
            fail(f'{label}: protein length {len(core)} outside 10-250'); errors += 1
        if cls.startswith('fab') and ':' not in seq:
            fail(f'{label}: fab must be VH:VL joined by colon'); errors += 1
        if cls == 'nanobody' and not (100 <= len(core) <= 160):
            warn(f'{label}: nanobody length {len(core)} unusual (expect ~110-140; ANARCI annotation decides)')
        for tname, motif in TAG_MOTIFS.items():
            if motif in core: warn(f'{label}: contains {tname} motif - tag fusions need manual review')
        hits = [n for n, ks in EXCLUDED_KMERS.items() if kmers(core) & ks]
        if hits:
            warn(f'{label}: shares >=8-aa exact window with {hits} - NOVELTY REVIEW (de novo rule prohibits starting from known binders)')
        print(f'  {i:2d}. {name}  {cls}  {len(core)} aa')
    print()
    if errors: print(f'{errors} FAIL(s) - fix before anyone considers submitting'); sys.exit(1)
    ok('format checks passed; WARNs need human review; official novelty filter cannot be run locally')

if __name__ == '__main__':
    if len(sys.argv) != 2: sys.exit('usage: validate_submission.py submission.csv')
    main(sys.argv[1])
