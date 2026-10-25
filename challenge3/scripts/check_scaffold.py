#!/usr/bin/env python3
# check_v3.py - Validacion geometrica profunda del scaffold (cadena B):
#   angulos N-CA-C / CA-C-N / C-N-CA (refs 111/116/121 +-3), torsion omega (planaridad,
#   ref |180| +-10), phi/psi (helicoidal ~ -57/-47 +-20), QUIRALIDAD en CA (signo del
#   diedro impropio comparado con una helice real de 6X18_R) y colisiones INTERNAS
#   (pares no enlazados de la cadena B a <2.0 A, excluyendo vecinos de enlace).
# Uso: python3 check_v3.py <scaffold.pdb> [referencia_real.pdb]
import math, sys

def cargar(path):
    A, B = {}, {}
    for line in open(path):
        if line.startswith('ATOM') and line[16] in (' ', 'A'):
            ch = line[21]; rs = int(line[22:26]); an = line[12:16].strip()
            d = A if ch == 'A' else (B if ch == 'B' else None)
            if d is not None: d.setdefault(rs, {})[an] = tuple(float(line[30+8*i:38+8*i]) for i in range(3))
    return A, B

def sub(a, b): return tuple(x-y for x, y in zip(a, b))
def dot(a, b): return sum(x*y for x, y in zip(a, b))
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a):
    n = math.sqrt(dot(a, a)); return (a[0]/n, a[1]/n, a[2]/n)

def angulo(p1, p2, p3):
    v1, v2 = norm(sub(p1, p2)), norm(sub(p3, p2))
    return math.degrees(math.acos(max(-1, min(1, dot(v1, v2)))))

def diedro(p1, p2, p3, p4):
    b0 = norm(sub(p1, p2)); b1 = norm(sub(p3, p2)); b2 = norm(sub(p4, p3))
    v = sub(b0, tuple(dot(b0, b1)*x for x in b1)); w = sub(b2, tuple(dot(b2, b1)*x for x in b1))
    x, y = dot(v, w), dot(cross(b1, v), w)
    return math.degrees(math.atan2(y, x))

A, B = cargar(sys.argv[1])
resid = sorted(r for r in B if all(k in B[r] for k in ('N', 'CA', 'C', 'O')))
mal = 0
print('== Angulos y torsiones ==')
for i, r in enumerate(resid):
    a1 = angulo(B[r]['N'], B[r]['CA'], B[r]['C'])
    flags = []
    if abs(a1 - 111.0) > 3: flags.append(f'N-CA-C {a1:.1f}'); mal += 1
    if i < len(resid) - 1:
        r2 = resid[i+1]
        a2 = angulo(B[r]['CA'], B[r]['C'], B[r2]['N'])
        a3 = angulo(B[r]['C'], B[r2]['N'], B[r2]['CA'])
        om = diedro(B[r]['CA'], B[r]['C'], B[r2]['N'], B[r2]['CA'])
        if abs(a2 - 116.0) > 3: flags.append(f'CA-C-N {a2:.1f}'); mal += 1
        if abs(a3 - 121.0) > 3: flags.append(f'C-N-CA {a3:.1f}'); mal += 1
        if abs(abs(om) - 180.0) > 10: flags.append(f'omega {om:.1f}'); mal += 1
    if 0 < i < len(resid) - 1:
        r0, r2 = resid[i-1], resid[i+1]
        phi = diedro(B[r0]['C'], B[r]['N'], B[r]['CA'], B[r]['C'])
        psi = diedro(B[r]['N'], B[r]['CA'], B[r]['C'], B[r2]['N'])
        if abs(phi + 57.8) > 20 or abs(psi + 47.0) > 20: flags.append(f'phi/psi {phi:.0f}/{psi:.0f}'); mal += 1
    if flags: print(f'  B{r}: {", ".join(flags)}')
print(f'Fuera de tolerancia: {mal}')

print('== Quiralidad (diedro impropio C-N-CA-CB esperado ~ -33.9 en L) ==')
def cb_l(n, ca, c):
    v1, v2 = norm(sub(n, ca)), norm(sub(c, ca))
    bis = norm(tuple(-(a+b) for a, b in zip(v1, v2)))
    perp = norm(cross(v1, v2))
    # dos soluciones tetraedricas; elegir la L: probar ambas y quedarnos con diedro impropio <0
    for s in (1, -1):
        cb = tuple(ca[i] + 1.53*(0.816*bis[i] + s*0.577*perp[i]) for i in range(3))
        imp = diedro(c, n, ca, cb)
        if imp < 0: return cb, imp
    return cb, imp
signos = set()
for r in resid:
    cb, imp = cb_l(B[r]['N'], B[r]['CA'], B[r]['C'])
    signos.add('L' if imp < 0 else 'D')
print('  clases de quiralidad en scaffold:', signos)
if len(sys.argv) > 2:
    BR = {}
    for line in open(sys.argv[2]):
        if line.startswith('ATOM') and line[16] in (' ', 'A'):
            rs = int(line[22:26]); an = line[12:16].strip()
            BR.setdefault(rs, {})[an] = tuple(float(line[30+8*i:38+8*i]) for i in range(3))
    rref = sorted(r for r in BR if all(k in BR[r] for k in ('N', 'CA', 'C', 'CB')))[10]
    imp_ref = diedro(BR[rref]['C'], BR[rref]['N'], BR[rref]['CA'], BR[rref]['CB'])
    print(f'  referencia real (res {rref}): impropio {imp_ref:.1f} ({"L" if imp_ref < 0 else "D"})')

print('== Colisiones internas cadena B (no enlazadas, <2.0 A) ==')
ENL = {('N', 'CA'), ('CA', 'C'), ('C', 'O')}
col = 0
for i, r1 in enumerate(resid):
    for j, r2 in enumerate(resid):
        if j < i: continue
        for a1, p1 in B[r1].items():
            for a2, p2 in B[r2].items():
                if i == j and ({a1, a2} == {a1} or frozenset((a1, a2)) in {frozenset(e) for e in ENL}): continue
                if j == i+1 and ((a1 == 'C' and a2 in ('N', 'CA')) or (a1 == 'O' and a2 == 'N')): continue
                if math.dist(p1, p2) < 2.0:
                    print(f'  B{r1}-{a1} vs B{r2}-{a2}: {math.dist(p1, p2):.2f}'); col += 1
print(f'Colisiones internas: {col}')
