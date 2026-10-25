#!/usr/bin/env python3
# orient.py - Analisis de orientacion extremo-a-extremo de un binder sobre GLP1R (crop ch3)
# Uso: python3 orient.py <complejo.pdb>
# Entrada: PDB con 1+ MODELs, cadena A = receptor (crop ch3, numeracion resseq original),
#          cadena B = binder (36 aa). Altloc: solo ' '/A; atomos duplicados descartados.
# Salida: por MODEL, lista de residuos del binder con contacto <=6 A a la cadena A,
#         clasificados por zona del receptor (ECD <=139, ECL loops, TM) segun UniProt P43220.
# Interpretacion: con el modo de union nativo como referencia (two-domain), el extremo N
# del peptido deberia asociarse al bolsillo de activacion (U367-391) y el C al ECD.
import math, sys

# resseq del crop -> numeracion UniProt P43220 (offset constante +28 dentro de cada fragmento)
FRAGS = [(1,23,29),(33,101,61),(108,136,136),(148,223,176),(238,300,266),(321,376,349)]
RES2UP = {r: u0 + (r - a) for a, b, u0 in FRAGS for r in range(a, b+1)}
ECL = set(range(202,228)) | set(range(291,306)) | set(range(371,384))  # ECL1, ECL2, ECL3 (UniProt)

def zona(u):
    if u <= 139: return 'ECD'
    if u in ECL: return 'ECL'
    if u >= 140: return 'TM'
    return '?'

models, cur = {}, None
seen = set()
for line in open(sys.argv[1]):
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
        models[cur][ch].setdefault(rs, []).append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

for mid in sorted(models):
    m = models[mid]
    zones = {}
    for rb in m['B']:
        best = (1e9, None)
        for rt in m['A']:
            for x, y, z in m['B'][rb]:
                for tx, ty, tz in m['A'][rt]:
                    d = math.sqrt((x-tx)**2 + (y-ty)**2 + (z-tz)**2)
                    if d < best[0]: best = (d, rt)
        if best[1] is not None and best[0] <= 6.0:
            u = RES2UP.get(best[1], 0)
            zones[rb] = (zona(u), round(best[0], 1), u)
    print('MODEL', mid)
    print('  binder_N(1-12):', [(r,) + zones[r] for r in sorted(zones) if r <= 12])
    print('  binder_C(27-36):', [(r,) + zones[r] for r in sorted(zones) if r >= 27])
