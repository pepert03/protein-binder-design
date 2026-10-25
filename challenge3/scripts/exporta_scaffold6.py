#!/usr/bin/env python3
# exporta_scaffold6.py - Scaffold v6: pose de muestreo_v6.py (packing poly-ALA en bucle).
# Mismo exportador verificado que v5 (geometria NeRF calibrada + Kabsch + crop FRAGS).
import math, os
import numpy as np
exec(open(os.environ.get('MUESTREO_V6', '/tmp/pose/muestreo_v6.py')).read().split("# ---- muestreo global")[0])  # build_helix, kabsch, trace_ca, pose_atoms, gate_geo, pack_choques, byres

PARAMS = ([109.981, 134.129, 111.779], [0.2222, 0.7393, -0.6357], 2.8347)
anclaje, e, roll = np.array(PARAMS[0]), np.array(PARAMS[1]), PARAMS[2]
pts = pose_atoms(anclaje, e, roll)
ch, n_bols, d_ecd = gate_geo(pts)
ncl = pack_choques(pts)
print(f"pose v6: bb_choques={ch} n_bols={n_bols} d_ecd={d_ecd:.2f} packing_polyALA={ncl}")

FRAGS = [(1,23,29),(33,101,61),(108,136,136),(148,223,176),(238,300,266),(321,376,349)]
UP2CROP = {u0 + (r - a): r for a, b, u0 in FRAGS for r in range(a, b+1)}
out = []; serial = 0
for line in open('/tmp/ch3pack/challenge3/target/6X18_RP.pdb'):
    if line.startswith('ATOM') and line[21] == 'R':
        u_res = int(line[22:26])
        if u_res not in UP2CROP: continue
        serial += 1
        out.append(f'ATOM  {serial:5d} {line[12:16]}{line[16]}{line[17:20]} A{UP2CROP[u_res]:4d}{line[26:54]} 1.00  0.00          {line[76:78]}')
ele = {'N':'N','CA':'C','C':'C','O':'O'}
for i in range(36):
    for j, an in enumerate(('N','CA','C','O')):
        p = pts[i*4+j]; serial += 1
        out.append('ATOM  %5d %-4s%s%-3s %s%4d%s   %8.3f%8.3f%8.3f  1.00  0.00          %2s' % (serial, an, ' ', 'ALA', 'B', i+1, ' ', p[0], p[1], p[2], ele[an]))
open(os.environ.get('V6_OUT', '/tmp/pose/scaffold_v6.pdb'),'w').write('MODEL        1\n' + '\n'.join(out) + '\nENDMDL\nEND\n')
print('scaffold_v6.pdb escrito:', len(out), 'atomos')
