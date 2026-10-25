#!/usr/bin/env python3
"""Mapa CB inequívoco con la fórmula EXACTA de ProteinMPNN (protein_mpnn_utils.py l.968-970,
commit del repo clonado; coeficientes publicados por dauparas).
Validación externa: aplicada a residuos reales del receptor 6X18, |CA-CB| debe ~1.53 A."""
import math

def read_atoms(path, altok=(' ', 'A')):
    atoms = []
    for line in open(path):
        if not line.startswith('ATOM'):
            continue
        if line[16] not in altok:
            continue
        atoms.append({'name': line[12:16].strip(), 'chain': line[21],
                      'resseq': int(line[22:26]), 'resname': line[17:20].strip(),
                      'x': float(line[30:38]), 'y': float(line[38:46]), 'z': float(line[46:54])})
    # dedupe (mismo atom name/residuo queda el primero = altloc ' ' o 'A')
    seen, out = set(), []
    for a in atoms:
        k = (a['chain'], a['resseq'], a['name'])
        if k not in seen:
            seen.add(k); out.append(a)
    return out

def cb_mpnn(N, CA, C):
    b = [CA[i]-N[i] for i in range(3)]
    c = [C[i]-CA[i] for i in range(3)]
    a = [b[1]*c[2]-b[2]*c[1], b[2]*c[0]-b[0]*c[2], b[0]*c[1]-b[1]*c[0]]
    return [CA[i] - 0.58273431*a[i] + 0.56802827*b[i] - 0.54067466*c[i] for i in range(3)]

atoms = read_atoms(__import__("sys").argv[1] if len(__import__("sys").argv)>1 else "scaffold_v6.pdb")
receptor = [a for a in atoms if a['chain'] == 'A']
binder = [a for a in atoms if a['chain'] == 'B']

# 1) VALIDACION EXTERNA: formula sobre el receptor real (6X18)
dists = []
byres = {}
for a in receptor:
    byres.setdefault(a['resseq'], {})[a['name']] = [a['x'], a['y'], a['z']]
for resseq, d in byres.items():
    if all(k in d for k in ('N', 'CA', 'C')) and d.get('resname') != 'GLY':
        cb = cb_mpnn(d['N'], d['CA'], d['C'])
        dists.append(math.dist(cb, d['CA']))
print(f"VALIDACION 6X18: n={len(dists)} residuos con N/CA/C | |CA-CB| min={min(dists):.3f} med={sorted(dists)[len(dists)//2]:.3f} max={max(dists):.3f} A (esperado ~1.53)")

# 2) MAPA CB binder: distancia minima CB->receptor por posicion
byres_b = {}
for a in binder:
    byres_b.setdefault(a['resseq'], {})[a['name']] = [a['x'], a['y'], a['z']]
print("\nMAPA CB (formula MPNN validada), distancia minima CB->atomo receptor:")
flags = []
for resseq in sorted(byres_b):
    d = byres_b[resseq]
    cb = cb_mpnn(d['N'], d['CA'], d['C'])
    dmin, atom = min(((math.dist(cb, [a['x'], a['y'], a['z']]), f"{a['name']} A{a['resseq']}") for a in receptor))
    tag = '  <<< <2.0' if dmin < 2.0 else ''
    if dmin < 2.0: flags.append((resseq, dmin, atom))
    print(f"  B{resseq:>2}: {dmin:5.2f} A  (vs {atom}){tag}")
print(f"\nCONFLICTOS CB<2.0 A (formula validada): {[(r, round(d,2)) for r, d, _ in flags]}")

# 3) contraste con la formula simplificada mia anterior (v=-(u+w) normalizado, 1.53 A)
def cb_simple(N, CA, C):
    u = [(N[i]-CA[i])/math.dist(N, CA) for i in range(3)]
    w = [(C[i]-CA[i])/math.dist(C, CA) for i in range(3)]
    s = [u[i]+w[i] for i in range(3)]; n = math.sqrt(sum(v*v for v in s))
    return [CA[i] - 1.53*s[i]/n for i in range(3)]
flags2 = []
for resseq in sorted(byres_b):
    d = byres_b[resseq]
    cb = cb_simple(d['N'], d['CA'], d['C'])
    dmin = min(math.dist(cb, [a['x'], a['y'], a['z']]) for a in receptor)
    if dmin < 2.0: flags2.append((resseq, round(dmin, 2)))
print(f"CONFLICTOS con formula simplificada anterior: {flags2}")
