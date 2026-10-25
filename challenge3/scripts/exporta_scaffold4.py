#!/usr/bin/env python3
# exporta_scaffold3.py - Scaffold v3: hélice ideal con GEOMETRIA REAL de backbone
# (NeRF: N-CA 1.46, CA-C 1.52, C-N 1.33, C-O 1.23; phi=-57.8 psi=-47.0 omega=180),
# ajustada por Kabsch (sobre CA) a la pose muestreada target-only. numpy local; la
# verificación del resultado es pura Python (verifica_scaffold_v3.py).
import math, random, os, sys
import numpy as np
random.seed(23); np.random.seed(23)
# argv[1] = 6X18_RP.pdb (lo usa muestreo_helice), argv[2] = salida (defecto scaffold_v4.pdb)
_mh = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'muestreo_helice.py')
exec(open(_mh).read().split('best = ')[0])  # R, ratoms, c_bol, c_ecd, ecd_atoms, e0

# ---- helice ideal con geometria real (NeRF simplificado) ----
B_N_CA, B_CA_C, B_C_N, B_C_O = 1.46, 1.52, 1.33, 1.23
A_NCA_C = math.radians(111.0); A_CA_C_N = math.radians(116.0); A_C_N_CA = math.radians(121.0)
PHI, PSI, OMEGA = math.radians(-57.8-180), math.radians(-47.0-180), math.radians(0.0)  # convencion calibrada vs diedro estandar y helice real 6X18

def nerf(a, b, c, l, theta, chi):
    # NeRF estandar, calibrado contra helice real 6X18 (phi -57.8 / psi -47 / omega +-180, CA-CA 3.80)
    bc = (c - b) / np.linalg.norm(c - b); ab = a - b
    n = np.cross(ab, bc); n /= np.linalg.norm(n); m = np.cross(n, bc)
    return c + (-l*math.cos(theta)*bc + l*math.sin(theta)*math.cos(chi)*m + l*math.sin(theta)*math.sin(chi)*n)

def build_helix(n=36):
    N = np.array([0.0, 1.46, 0.0]); CA = np.array([0.0, 0.0, 0.0]); C = np.array([1.52, 0.0, 0.0])
    out = []
    for i in range(n):
        O = C + 1.23 * (C - CA) / np.linalg.norm(C - CA)  # O provisional; se corrige abajo
        out.append({'N': N.copy(), 'CA': CA.copy(), 'C': C.copy(), 'O': O})
        if i < n - 1:
            Nn = nerf(N, CA, C, B_C_N, A_CA_C_N, PSI)
            CAn = nerf(CA, C, Nn, B_N_CA, A_C_N_CA, OMEGA)
            Cn = nerf(C, Nn, CAn, B_CA_C, A_NCA_C, PHI)
            N, CA, C = Nn, CAn, Cn
    # O correcto: bisectriz CA-C-N(+1)
    for i, r in enumerate(out):
        if i < n - 1:
            d1 = (r['CA'] - r['C'])/np.linalg.norm(r['CA'] - r['C']); d2 = (out[i+1]['N'] - r['C'])/np.linalg.norm(out[i+1]['N'] - r['C'])
        else:
            d1 = (r['CA'] - r['C'])/np.linalg.norm(r['CA'] - r['C']); d2 = d1
        r['O'] = r['C'] - B_C_O * (d1 + d2) / np.linalg.norm(d1 + d2)  # O opuesto a (CA + N_sig)
    return out

HELIX = build_helix()
CA_CANON = np.array([r['CA'] for r in HELIX])
ALL_CANON = np.array([r[a] for r in HELIX for a in ('N', 'CA', 'C', 'O')])

def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt

def trace_ca(anclaje, e, roll, n=36):
    u = np.cross(e, [1,0,0])
    if np.linalg.norm(u) < 1e-6: u = np.cross(e, [0,1,0])
    u /= np.linalg.norm(u); v = np.cross(e, u)
    return np.array([anclaje + e*(i*1.5) + 2.3*(math.cos(math.radians(100*i)+roll)*u + math.sin(math.radians(100*i)+roll)*v) for i in range(n)])

def pose_atoms(anclaje, e, roll):
    tr = trace_ca(anclaje, e, roll)
    Rk = kabsch(CA_CANON, tr)
    return (ALL_CANON - CA_CANON.mean(0)) @ Rk + tr.mean(0)

bol_atoms = np.array([q for r in range(367,392) if r in R for q in R[r].values()])
def score(pts):
    d = np.sqrt(((pts[:,None,:] - ratoms[None,:,:])**2).sum(-1))
    ch = int((d < 2.0).sum())
    dn = np.sqrt(((pts[:48][:,None,:] - bol_atoms[None,:,:])**2).sum(-1))
    n_bols = int(sum(1 for i in range(12) if dn[i*4:(i+1)*4].min() <= 6.0))
    d_ecd = np.sqrt(((pts[-4:][:,None,:] - ecd_atoms[None,:,:])**2).sum(-1)).min()
    return (ch, -n_bols, d_ecd)

best = ((1e9, 0, 1e9), None)
for k in range(2200):
    anclaje = c_bol + np.random.uniform(-4, 4, 3)
    tilt = math.radians(random.uniform(0, 40))
    ar = np.random.normal(size=3)
    e = e0*math.cos(tilt) + ar/np.linalg.norm(ar)*math.sin(tilt); e /= np.linalg.norm(e)
    roll = random.uniform(0, 2*math.pi)
    s = score(pose_atoms(anclaje, e, roll))
    s = (s[0], -s[1], s[2]) if (-s[1] >= 3 and s[2] <= 6.0) else (1e6 + s[0], -s[1], s[2])  # restriccion dura
    if s < best[0]: best = (s, (anclaje.copy(), e.copy(), roll))
if best[1] is None: raise SystemExit('ninguna pose cumple n_bols>=3 y C-ECD<=6')
a0, e0b, r0 = best[1]
for k in range(2200):
    anclaje = a0 + np.random.uniform(-1.5, 1.5, 3)
    tilt = math.radians(random.uniform(-8, 8))
    ar = np.random.normal(size=3)
    e = e0b*math.cos(tilt) + ar/np.linalg.norm(ar)*math.sin(tilt); e /= np.linalg.norm(e)
    roll = r0 + random.uniform(-0.35, 0.35)
    s = score(pose_atoms(anclaje, e, roll))
    s = (s[0], -s[1], s[2]) if (-s[1] >= 3 and s[2] <= 6.0) else (1e6 + s[0], -s[1], s[2])
    if s < best[0]: best = (s, (anclaje.copy(), e.copy(), roll)); a0, e0b, r0 = best[1]
(ch, nnb, d_ecd), (anclaje, e, roll) = best
print(f'pose v4: choques={ch}, N-bolsillo={-nnb}, C-ECD={d_ecd:.1f} A')
print('PARAMS_V4:', repr([[float(x) for x in anclaje], [float(x) for x in e], float(roll)]))

# exportar con backbone real
FRAGS = [(1,23,29),(33,101,61),(108,136,136),(148,223,176),(238,300,266),(321,376,349)]
UP2CROP = {u0 + (r - a): r for a, b, u0 in FRAGS for r in range(a, b+1)}
tr = trace_ca(anclaje, e, roll)
Rk = kabsch(CA_CANON, tr)
placed = (ALL_CANON - CA_CANON.mean(0)) @ Rk + tr.mean(0)
out = []; serial = 0
for line in open(sys.argv[1] if len(sys.argv) > 1 else '6X18_RP.pdb'):
    if line.startswith('ATOM') and line[21] == 'R':
        u_res = int(line[22:26])
        if u_res not in UP2CROP: continue
        serial += 1
        out.append(f'ATOM  {serial:5d} {line[12:16]}{line[16]}{line[17:20]} A{UP2CROP[u_res]:4d}{line[26:54]} 1.00  0.00          {line[76:78]}')
ele = {'N':'N','CA':'C','C':'C','O':'O'}
for i in range(36):
    for j, an in enumerate(('N','CA','C','O')):
        p = placed[i*4+j]; serial += 1
        out.append('ATOM  %5d %-4s%s%-3s %s%4d%s   %8.3f%8.3f%8.3f  1.00  0.00          %2s' % (serial, an, ' ', 'ALA', 'B', i+1, ' ', p[0], p[1], p[2], ele[an]))
open(sys.argv[2] if len(sys.argv) > 2 else 'scaffold_v4.pdb','w').write('MODEL        1\n' + '\n'.join(out) + '\nENDMDL\nEND\n')
print('scaffold_v4.pdb escrito:', len(out), 'atomos')
