#!/usr/bin/env python3
"""FASPR con receptor como CONTEXTO FIJO (secuencia receptor en minusculas) + binder empaquetado.
Choques contados contra el receptor ORIGINAL (coords 6X18 exactas, no la salida FASPR).
Criterio: solape vdW >0.4 A. Control positivo: nativo 6X18 P vs R = 3 choques, dmin 2.46."""
import subprocess, os, sys
# Rutas parametrizadas por el verificador (25 oct). Uso: pack_eval_v6.py [BASE] [FASPR_DIR]
BASE = sys.argv[1] if len(sys.argv)>1 else '/tmp/pose'
FDIR = sys.argv[2] if len(sys.argv)>2 else BASE + '/FASPR'
import numpy as np

VDW = {'N':1.55, 'C':1.70, 'O':1.52, 'S':1.80}
def elem(n): return n.strip()[0] if n.strip()[0] in VDW else 'C'
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
def load(path, chain):
    seen, out = set(), []
    for line in open(path):
        if line.startswith('ATOM') and line[16] in (' ', 'A') and line[21] == chain:
            k = (int(line[22:26]), line[12:16].strip())
            if k not in seen:
                seen.add(k)
                out.append((int(line[22:26]), line[12:16].strip(), line[17:20].strip(), [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return out

receptor = load('complejo_v6.pdb', 'A')
rec_xyz = np.array([a[3] for a in receptor]); rec_r = np.array([VDW[elem(a[1])] for a in receptor])
rec_lower = open('recseq_v6_lower.txt').read().strip()
os.makedirs(BASE + '/packed_v6', exist_ok=True)
rows = []
for tag, path in (('LIBRE6', BASE + '/mpnn6/out/seqs/scaffold_v6.fa'),):
    for name, seq in lee_fasta(path):
        sid = [p.split('=')[1] for p in name.split(',') if 'sample' in p][0]
        open(BASE + '/seq_ctx.txt', 'w').write(rec_lower + seq + '\n')
        subprocess.run(['./faspr', '-i', BASE + '/complejo_v6.pdb', '-o', BASE + f'/packed_v6/{tag}_{sid}.pdb', '-s', BASE + '/seq_ctx.txt'],
                       cwd=FDIR, capture_output=True, check=True)
        binder = load(BASE + f'/packed_v6/{tag}_{sid}.pdb', 'B')
        b_xyz = np.array([a[3] for a in binder]); b_r = np.array([VDW[elem(a[1])] for a in binder])
        d = np.sqrt(((b_xyz[:, None, :] - rec_xyz[None, :, :])**2).sum(-1))
        clash = d < (b_r[:, None] + rec_r[None, :] - 0.4)
        i, j = np.unravel_index(np.argmin(d), d.shape)
        pos_clash = sorted(set(binder[k][0] for k in range(len(binder)) if clash[k].any()))
        rows.append((tag, sid, int(clash.sum()), float(d.min()), f"{binder[i][1]} B{binder[i][0]}-{receptor[j][1]} A{receptor[j][0]}", pos_clash))
        print(f"{tag} {sid:>9}: choques={int(clash.sum()):>3} dmin={d.min():.2f} ({binder[i][1]} B{binder[i][0]}-{receptor[j][1]} A{receptor[j][0]}) pos={pos_clash}")
print('\nRESUMEN (nativo=3 choques, dmin 2.46):')
for tag, sid, ncl, dmin, par, pc in sorted(rows, key=lambda x: x[2]):
    print(f"  {tag} {sid:>9} choques={ncl:>3} dmin={dmin:.2f}")
