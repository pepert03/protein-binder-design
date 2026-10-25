#!/usr/bin/env python3
# ensemble_rmsd.py - RMSD entre miembros del ensemble (MODELs) de un PDB de ColabDesign
# Uso: python3 ensemble_rmsd.py <complejo.pdb>
# Salida: RMSD global (todas las cadenas) y RMSD del binder (cadena B), sobre CA,
#         tras superposicion de Kabsch sobre CA del RECEPTOR (cadena A) - asi el RMSD
#         del binder refleja diferencia de pose, no solo de conformacion interna.
# Altloc solo ' '/A; dedup atomos; requiere >=2 MODELs con mismas cadenas.
import math, sys

def leer(path):
    models, cur, seen = {}, None, set()
    for line in open(path):
        if line.startswith('MODEL'):
            cur = int(line[10:14]); models[cur] = {'A': {}, 'B': {}}; seen = set()
        elif line.startswith('ENDMDL'):
            cur = None
        elif line.startswith('ATOM') and cur is not None:
            if line[16] not in (' ', 'A'): continue
            ch = line[21]
            if ch not in ('A', 'B'): continue
            rs = int(line[22:26]); an = line[12:16].strip()
            if (ch, rs, an) in seen: continue
            seen.add((ch, rs, an))
            if an == 'CA':
                models[cur][ch][rs] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return models

def kabsch_rmsd(P, Q, fit_idx, eval_idx):
    # P,Q: listas de coords (misma ordenacion). fit_idx: residuos para superponer; eval_idx: para RMSD.
    import numpy as np
    P = np.array(P); Q = np.array(Q)
    Pf, Qf = P[fit_idx], Q[fit_idx]
    Pc, Qc = Pf - Pf.mean(0), Qf - Qf.mean(0)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1, 1, d])
    R = V @ D @ Wt
    Pa = (P - Pf.mean(0)) @ R
    Qa = Q - Qf.mean(0)
    diff = Pa[eval_idx] - Qa[eval_idx]
    return math.sqrt((diff**2).sum() / len(eval_idx))

m = leer(sys.argv[1])
ids = sorted(m)
assert len(ids) >= 2, 'se requieren >=2 MODELs'
a, b = m[ids[0]], m[ids[1]]
comunes_a = sorted(set(a['A']) & set(a['A'].keys()) & set(b['A']))
comunes_b = sorted(set(a['B']) & set(b['B']))
comunes = [('A', r) for r in comunes_a] + [('B', r) for r in comunes_b]
P = [a[ch][r] for ch, r in comunes]; Q = [b[ch][r] for ch, r in comunes]
fit = list(range(len(comunes_a)))                       # CA del receptor
todos = list(range(len(comunes)))
bind = list(range(len(comunes_a), len(comunes)))
print(f'MODELs {ids[0]} vs {ids[1]}: receptor CA={len(comunes_a)}, binder CA={len(comunes_b)}')
print(f'RMSD global (fit receptor): {kabsch_rmsd(P, Q, fit, todos):.2f} A')
print(f'RMSD binder (fit receptor): {kabsch_rmsd(P, Q, fit, bind):.2f} A')
