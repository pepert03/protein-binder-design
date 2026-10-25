#!/usr/bin/env python3
# orient_gate.py - Gate de orientacion extremo-a-extremo (criterio de priorizacion de ESTA
# estrategia, no requisito universal). Uso: python3 orient_gate.py <complejo.pdb>
# Criterios (por MODEL; PASS requiere que TODOS los MODELs pasen):
#   G1: >=3 residuos de binder 1-12 con contacto <=6 A al bolsillo de activacion (U367-391)
#   G2: >=3 residuos de binder 27-36 con contacto <=6 A al ECD (U<=139)
#   G3: 0 residuos de binder 27-36 con contacto <=3 A al core TM (U140-404, excl. ECL)
#   G4 (esterico): 0 pares de ATOMOS intercadena <2 A (misma semantica que gate_pose.py)
# Zonas UniProt P43220. Altloc solo ' '/A; dedup. Numeracion resseq original + offset +28.
import math, sys

FRAGS = [(1,23,29),(33,101,61),(108,136,136),(148,223,176),(238,300,266),(321,376,349)]
RES2UP = {r: u0 + (r - a) for a, b, u0 in FRAGS for r in range(a, b+1)}
ECL = set(range(202,228)) | set(range(291,306)) | set(range(371,384))
BOLSILLO = set(range(367,392))  # bolsillo de activacion nativo (TM6/ECL3)

def cargar(path):
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
            models[cur][ch].setdefault(rs, []).append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return models

def dmin(p1, p2):
    m = 1e9
    for x, y, z in p1:
        for tx, ty, tz in p2:
            d = math.sqrt((x-tx)**2 + (y-ty)**2 + (z-tz)**2)
            if d < m: m = d
    return m

todo_ok = True
for mid, m in sorted(cargar(sys.argv[1]).items()):
    choques = 0
    for rb in m['B']:
        for rt in m['A']:
            for p1 in m['B'][rb]:
                for p2 in m['A'][rt]:
                    if math.dist(p1, p2) < 2.0: choques += 1
    n_bols = sum(1 for rb in m['B'] if rb <= 12 and any(dmin(m['B'][rb], m['A'][rt]) <= 6.0 for rt in m['A'] if RES2UP.get(rt, 0) in BOLSILLO))
    c_ecd = sum(1 for rb in m['B'] if rb >= 27 and any(dmin(m['B'][rb], m['A'][rt]) <= 6.0 for rt in m['A'] if 0 < RES2UP.get(rt, 999) <= 139))
    c_tm_duro = [rb for rb in m['B'] if rb >= 27 and any(dmin(m['B'][rb], m['A'][rt]) <= 3.0 for rt in m['A'] if RES2UP.get(rt, 0) >= 140 and RES2UP.get(rt, 0) not in ECL)]
    g1, g2, g3, g4 = n_bols >= 3, c_ecd >= 3, len(c_tm_duro) == 0, choques == 0
    ok = g1 and g2 and g3 and g4
    todo_ok = todo_ok and ok
    print(f'MODEL {mid}: G1_N-bolsillo={n_bols}({"OK" if g1 else "FALLO"}) G2_C-ECD={c_ecd}({"OK" if g2 else "FALLO"}) G3_C-en-TM={c_tm_duro}({"OK" if g3 else "FALLO"}) G4_choques={choques}({"OK" if g4 else "FALLO"}) -> {"PASS" if ok else "FAIL"}')
print('GATE GLOBAL:', 'PASS' if todo_ok else 'FAIL')
