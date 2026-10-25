#!/usr/bin/env python3
# exporta_scaffold2.py - Igual que exporta_scaffold.py pero el MODELO DE BACKBONE del
# muestreo es identico al exportado (N,CA,C,O) y el desempate premia contactos N-bolsillo.
# Determinista (seed fija). Exporta scaffold_v2.pdb listo para orient_gate.py.
import math, random, os, sys
import numpy as np
random.seed(23); np.random.seed(23)
# carga R, ratoms, c_bol, c_ecd, ecd_atoms, e0 desde muestreo_helice.py (mismo directorio);
# argv[1] = 6X18_RP.pdb (lo usa muestreo_helice), argv[2] = salida (defecto scaffold_v2.pdb)
_mh = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'muestreo_helice.py')
exec(open(_mh).read().split('best = ')[0])

FRAGS = [(1,23,29),(33,101,61),(108,136,136),(148,223,176),(238,300,266),(321,376,349)]
UP2CROP = {u0 + (r - a): r for a, b, u0 in FRAGS for r in range(a, b+1)}
bol_atoms = np.array([q for r in range(367,392) if r in R for q in R[r].values()])

def helice4(anclaje, e, roll, n=36):
    u = np.cross(e, [1,0,0])
    if np.linalg.norm(u) < 1e-6: u = np.cross(e, [0,1,0])
    u /= np.linalg.norm(u); v = np.cross(e, u)
    pts = []
    for i in range(n):
        t = i*1.5; ang = math.radians(100*i) + roll
        ca = anclaje + e*t + 2.3*(math.cos(ang)*u + math.sin(ang)*v)
        pts += [ca - 1.32*e, ca, ca + 1.52*e, ca + 1.52*e + 0.8*v]
    return np.array(pts)

def score4(pts):
    d = np.sqrt(((pts[:,None,:] - ratoms[None,:,:])**2).sum(-1))
    ch = int((d < 2.0).sum())
    nterm = pts[:12*4]
    dn = np.sqrt(((nterm[:,None,:] - bol_atoms[None,:,:])**2).sum(-1))
    n_bols = int(sum(1 for i in range(12) if dn[i*4:(i+1)*4].min() <= 6.0))
    cterm = pts[-4:]
    d_ecd = np.sqrt(((cterm[:,None,:] - ecd_atoms[None,:,:])**2).sum(-1)).min()
    return (ch, -n_bols, d_ecd)

best = ((1e9, 0, 1e9), None)
for k in range(2500):
    anclaje = c_bol + np.random.uniform(-4, 4, 3)
    tilt = math.radians(random.uniform(0, 40))
    ar = np.random.normal(size=3)
    e = e0*math.cos(tilt) + ar/np.linalg.norm(ar)*math.sin(tilt); e /= np.linalg.norm(e)
    roll = random.uniform(0, 2*math.pi)
    s = score4(helice4(anclaje, e, roll))
    if s < best[0]: best = (s, (anclaje.copy(), e.copy(), roll))
a0, e0b, r0 = best[1]
for k in range(2500):
    anclaje = a0 + np.random.uniform(-1.5, 1.5, 3)
    tilt = math.radians(random.uniform(-8, 8))
    ar = np.random.normal(size=3)
    e = e0b*math.cos(tilt) + ar/np.linalg.norm(ar)*math.sin(tilt); e /= np.linalg.norm(e)
    roll = r0 + random.uniform(-0.35, 0.35)
    s = score4(helice4(anclaje, e, roll))
    if s < best[0]: best = (s, (anclaje.copy(), e.copy(), roll)); a0, e0b, r0 = best[1]
(ch, nnb, d_ecd), (anclaje, e, roll) = best
print(f'pose final: choques={ch}, N-bolsillo={-nnb}, C-ECD={d_ecd:.1f} A')
print('PARAMS:', repr([list(anclaje), list(e), roll]))

u = np.cross(e, [1,0,0])
if np.linalg.norm(u) < 1e-6: u = np.cross(e, [0,1,0])
u /= np.linalg.norm(u); v = np.cross(e, u)
out = []; serial = 0; skipped = 0
for line in open('/tmp/ch3pack/challenge3/target/6X18_RP.pdb'):
    if line.startswith('ATOM') and line[21] == 'R':
        u_res = int(line[22:26])
        if u_res not in UP2CROP: skipped += 1; continue
        serial += 1
        out.append(f'ATOM  {serial:5d} {line[12:16]}{line[16]}{line[17:20]} A{UP2CROP[u_res]:4d}{line[26:54]} 1.00  0.00          {line[76:78]}')
ele = {'N':'N','CA':'C','C':'C','O':'O'}
for i in range(36):
    t = i*1.5; ang = math.radians(100*i) + roll
    ca = anclaje + e*t + 2.3*(math.cos(ang)*u + math.sin(ang)*v)
    atoms = {'N': ca - 1.32*e, 'CA': ca, 'C': ca + 1.52*e, 'O': ca + 1.52*e + 0.8*v}
    for an, p in atoms.items():
        serial += 1
        out.append('ATOM  %5d %-4s%s%-3s %s%4d%s   %8.3f%8.3f%8.3f  1.00  0.00          %2s' % (serial, an, ' ', 'ALA', 'B', i+1, ' ', p[0], p[1], p[2], ele[an]))
open(sys.argv[2] if len(sys.argv) > 2 else 'scaffold_v2.pdb','w').write('MODEL        1\n' + '\n'.join(out) + '\nENDMDL\nEND\n')
print('scaffold escrito:', sys.argv[2] if len(sys.argv) > 2 else 'scaffold_v2.pdb', '|', len(out), 'atomos; receptor saltados fuera del crop:', skipped)
