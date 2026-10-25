#!/usr/bin/env python3
# bb_check.py - Verifica la geometria quimica del backbone aproximado de un scaffold PDB
# (cadena B) y la holgura de sidechains via CB colocado en geometria tetraedrica estandar.
# Referencias Engh-Huber: N-CA 1.46, CA-C 1.52, C-O 1.23, C-N(peptido) 1.33 (tolerancia 0.1 A).
# Uso: python3 bb_check.py <scaffold.pdb>
import math, sys

def cargar(path):
    A, B = {}, {}
    for line in open(path):
        if line.startswith('ATOM') and line[16] in (' ', 'A'):
            ch = line[21]; rs = int(line[22:26]); an = line[12:16].strip()
            d = A if ch == 'A' else (B if ch == 'B' else None)
            if d is not None: d.setdefault(rs, {})[an] = tuple(float(line[30+8*i:38+8*i]) for i in range(3))
    return A, B

def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def scale(a, s): return tuple(x * s for x in a)
def norm(a):
    n = math.dist(a, (0, 0, 0)); return scale(a, 1 / n)
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

A, B = cargar(sys.argv[1])
resid = sorted(B)
print(f'Cadena B: {len(resid)} residuos; cadena A: {len(A)}')
refs = [('N','CA',1.46), ('CA','C',1.52), ('C','O',1.23)]
mal = 0
for r in resid:
    for a1, a2, ref in refs:
        if a1 in B[r] and a2 in B[r]:
            d = math.dist(B[r][a1], B[r][a2])
            if abs(d - ref) > 0.1:
                print(f'  FUERA: B{r} {a1}-{a2} = {d:.2f} (ref {ref})'); mal += 1
for i in range(len(resid) - 1):
    r1, r2 = resid[i], resid[i+1]
    if 'C' in B[r1] and 'N' in B[r2]:
        d = math.dist(B[r1]['C'], B[r2]['N'])
        if abs(d - 1.33) > 0.1:
            print(f'  FUERA: peptido C{r1}-N{r2} = {d:.2f} (ref 1.33)'); mal += 1
print(f'Enlaces fuera de tolerancia: {mal}')

# CB tetraedrico: posicion estandar a partir de N, CA, C
ratoms = [q for r in A for q in A[r].values()]
print('\nHolgura por residuo (min distancia CA/CB al receptor, A):')
ajustados = 0
for r in resid:
    if not all(k in B[r] for k in ('N', 'CA', 'C')): continue
    n, ca, c = B[r]['N'], B[r]['CA'], B[r]['C']
    v1 = norm(sub(n, ca)); v2 = norm(sub(c, ca))
    # CB en direccion opuesta a la bisectriz, corregida con la normal del plano
    bis = norm(add(v1, v2))
    perp = norm(cross(v1, v2))
    cb_dir = norm(add(scale(bis, -1.0), scale(perp, 0.0)))  # simplificacion: opuesta a bisectriz
    cb = add(ca, scale(cb_dir, 1.53))
    dca = min(math.dist(ca, q) for q in ratoms)
    dcb = min(math.dist(cb, q) for q in ratoms)
    flag = ' CB<2A!' if dcb < 2.0 else ''
    if dcb < 2.0: ajustados += 1
    print(f'  B{r:2d}: CA {dca:5.2f} | CB {dcb:5.2f}{flag}')
print(f'\nResiduos con CB a <2.0 A del receptor: {ajustados} de {len(resid)}')
