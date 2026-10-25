#!/usr/bin/env python3
"""Empaqueta sidechains de candidatas con FASPR y evalua choques all-atom binder-receptor.
Criterio choque: dist < r_vdw_i + r_vdw_j - 0.4 A (solape vdW >0.4). Tambien cuenta crudos <2.0 A.
Sidechains MODELADOS (rotamer library), backbone fijo, receptor real fijo: cribado, no validacion."""
import subprocess, math, os, sys
# Rutas parametrizadas por el verificador (25 oct): BASE = dir de trabajo con bb_binder.pdb,
# receptor_A.pdb, mpnn/, FASPR/. Uso: python3 pack_eval.py [BASE]  (default /tmp/pose)
BASE = sys.argv[1] if len(sys.argv) > 1 else '/tmp/pose'
import numpy as np

VDW = {'N':1.55, 'C':1.70, 'O':1.52, 'S':1.80}
def elem(name): 
    c = name.strip()[0]
    return c if c in VDW else 'C'

def lee_fasta(path):
    recs, name, seq = [], None, ''
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if name: recs.append((name, seq))
            name, seq = line[1:], ''
        else: seq += line
    if name: recs.append((name, seq))
    return recs[1:]

def load_atoms(path):
    seen, out = set(), []
    for line in open(path):
        if line.startswith('ATOM') and line[16] in (' ', 'A'):
            k = (line[21], int(line[22:26]), line[12:16].strip())
            if k not in seen:
                seen.add(k)
                out.append({'chain': line[21], 'resseq': int(line[22:26]), 'resname': line[17:20].strip(),
                            'name': line[12:16].strip(), 'xyz': [float(line[30:38]), float(line[38:46]), float(line[46:54])]})
    return out

receptor = load_atoms(BASE + '/receptor_A.pdb')
rec_xyz = np.array([a['xyz'] for a in receptor]); rec_r = np.array([VDW[elem(a['name'])] for a in receptor])

os.makedirs(BASE + '/packed', exist_ok=True)
rows = []
for tag, path in (('LIBRE', BASE + '/mpnn/out/seqs/scaffold_v5.fa'), ('RESTR', BASE + '/mpnn/out_restr/seqs/scaffold_v5.fa')):
    for name, seq in lee_fasta(path):
        sid = [p.split('=')[1] for p in name.split(',') if 'sample' in p][0]
        open(BASE + '/seq_tmp.txt', 'w').write(seq + '\n')
        subprocess.run(['./faspr', '-i', BASE + '/bb_binder.pdb', '-o', fBASE + '/packed/{tag}_{sid}.pdb', '-s', BASE + '/seq_tmp.txt'],
                       cwd=BASE + '/FASPR', capture_output=True, check=True)
        binder = load_atoms(fBASE + '/packed/{tag}_{sid}.pdb')
        b_xyz = np.array([a['xyz'] for a in binder]); b_r = np.array([VDW[elem(a['name'])] for a in binder])
        # distancias binder x receptor
        d = np.sqrt(((b_xyz[:, None, :] - rec_xyz[None, :, :])**2).sum(-1))
        thr = b_r[:, None] + rec_r[None, :] - 0.4
        clash = d < thr
        ncl = int(clash.sum())
        i, j = np.unravel_index(np.argmin(d), d.shape)
        # choques por posicion binder
        pos_idx = {}
        for bi, a in enumerate(binder):
            pos_idx.setdefault(a['resseq'], []).append(bi)
        per_pos = {r: int(clash[idx].any()) for r, idx in pos_idx.items()}
        pos_clash = sorted(r for r, v in per_pos.items() if v)
        rows.append((tag, sid, ncl, float(d.min()), f"{binder[i]['name']} B{binder[i]['resseq']}-{receptor[j]['name']} A{receptor[j]['resseq']}", pos_clash))
        print(f"{tag} {sid:>9}: choques_vdW={ncl:>3}  dmin={d.min():.2f} ({binder[i]['name']} B{binder[i]['resseq']} - {receptor[j]['name']} A{receptor[j]['resseq']})  pos_con_choque={pos_clash}")
print('\nRESUMEN ordenado por choques:')
for tag, sid, ncl, dmin, par, pc in sorted(rows, key=lambda x: x[2]):
    print(f"  {tag} {sid:>9} choques={ncl:>3} dmin={dmin:.2f}")
