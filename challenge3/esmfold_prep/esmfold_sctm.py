#!/usr/bin/env python3
# esmfold_sctm.py - auto-consistencia estructural AISLADA con ESMFold (facebook/esmfold_v1).
# Pliega cada secuencia SOLA y compara el CA-trace predicho con el backbone disenado (Kabsch CA-RMSD).
# QUE CONTESTA: si la secuencia, aislada, prefiere el pliegue disenado segun ESMFold.
# QUE NO CONTESTA: union a GLP1R, interfaz, agonismo, ni pliegue en presencia del receptor.
# Uso: python3 esmfold_sctm.py --fasta d.fa [...] --ref scaffold_v6.pdb --ref-n-ca 36 --out salida/
# Requiere GPU (>=16 GB recomendado, fp16) — NO corre en esta CPU (1 GB RAM).
import argparse, json, math, os
import numpy as np

def ca_trace_pdb_text(text):
    pts = []
    for line in text.splitlines():
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(pts)

def plddt_pdb_text(text):
    vals = [float(line[60:66]) for line in text.splitlines() if line.startswith('ATOM') and line[12:16].strip() == 'CA']
    return sum(vals) / len(vals) if vals else float('nan')

def ca_trace_pdb_file(path, n_ca):
    # agrupa por ID de cadena (col 22); NO depender de TER (scaffold_v*.pdb no los escribe)
    ch = {}
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            ch.setdefault(line[21], []).append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    hits = [v for v in ch.values() if len(v) == n_ca]
    assert len(hits) == 1, f'esperaba 1 cadena con {n_ca} CA en {path}, hay {len(hits)}: {[(k, len(v)) for k, v in ch.items()]}'
    return np.array(hits[0])

def kabsch_rmsd(P, Q):
    assert P.shape == Q.shape and P.shape[0] > 2
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    U = V @ D @ Wt
    P_rot = Pc @ U
    diff = P_rot - Qc
    per = np.sqrt((diff ** 2).sum(1))
    return float(np.sqrt((diff ** 2).sum() / P.shape[0])), float(per.max())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', nargs='+', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--ref-n-ca', type=int, required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--fp16', action='store_true', help='media precision (recomendado en GPU 16GB)')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    import torch
    from transformers import EsmForProteinFolding
    model = EsmForProteinFolding.from_pretrained('facebook/esmfold_v1')
    if a.fp16:
        model = model.half()
        model.esm = model.esm.float()  # LM en fp32, trunk fp16 (patron del tutorial HF)
    model = model.cuda().eval()
    ref = ca_trace_pdb_file(a.ref, a.ref_n_ca)
    for fa in a.fasta:
        name = os.path.splitext(os.path.basename(fa))[0]
        seq = ''.join(l.strip() for l in open(fa) if not l.startswith('>'))
        assert len(seq) == a.ref_n_ca, f'{name}: len {len(seq)} != ref {a.ref_n_ca}'
        with torch.no_grad():
            pdb_str = model.infer_pdb(seq)
        outp = os.path.join(a.out, f'{name}_esmfold.pdb')
        open(outp, 'w').write(pdb_str)
        pred = ca_trace_pdb_text(pdb_str)
        r, mx = kabsch_rmsd(pred, ref)
        rec = {'name': name, 'len': len(seq), 'mean_plddt': round(plddt_pdb_text(pdb_str), 2),
               'ca_rmsd_vs_ref': round(r, 3), 'max_ca_dev': round(mx, 3), 'pdb': outp}
        print(json.dumps(rec), flush=True)

if __name__ == '__main__':
    main()
