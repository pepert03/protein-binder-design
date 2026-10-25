#!/usr/bin/env python3
# muestreo_helice.py - Prueba de existencia: busca por muestreo Monte Carlo una pose de
# helice ideal (36 aa) LIBRE DE CHOQUES de backbone usando SOLO coordenadas del receptor.
# DOF: punto de anclaje N en el bolsillo (U367-391), direccion del eje (cono alrededor del
# eje bolsillo->ECD), giro. Objetivo: 0 choques <2.0 A; desempate: extremo C cerca del ECD.
# Si existe, el scaffolding por extremos es viable sin tocar coordenadas del ligando.
# LIMITE: pseudo-backbone aproximado (CA + N/C estimados) sin sidechains del binder; los
# conteos NO son comparables con el gate all-atom (gate_pose.py) de los disenos.
# Uso: python3 muestreo_helice.py <6X18_RP.pdb>
import math, random, sys
import numpy as np

random.seed(7); np.random.seed(7)

def cargar_R(path):
    R = {}
    for line in open(path):
        if line.startswith('ATOM') and line[21] == 'R' and line[17:20].strip() != 'HOH':
            rs = int(line[22:26]); an = line[12:16].strip()
            R.setdefault(rs, {})[an] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return R

R = cargar_R(sys.argv[1] if len(sys.argv) > 1 else '6X18_RP.pdb')
ratoms = np.array([q for r in R for q in R[r].values()])

def centroid(residues):
    pts = np.array([R[r]['CA'] for r in residues if r in R and 'CA' in R[r]])
    return pts.mean(0)

c_bol = centroid(range(367,392)); c_ecd = centroid(range(120,139))
ecd_atoms = np.array([q for r in range(29,140) if r in R for q in R[r].values()])
e0 = c_ecd - c_bol; e0 = e0/np.linalg.norm(e0)

def helice(anclaje, e, roll, n=36):
    u = np.cross(e, [1,0,0])
    if np.linalg.norm(u) < 1e-6: u = np.cross(e, [0,1,0])
    u = u/np.linalg.norm(u); v = np.cross(e, u)
    pts = []
    for i in range(n):
        t = i*1.5; ang = math.radians(100*i) + roll
        ca = anclaje + e*t + 2.3*(math.cos(ang)*u + math.sin(ang)*v)
        pts += [ca, ca - 1.3*e, ca + 1.3*e]
    return np.array(pts)

def score(pts):
    d = np.sqrt(((pts[:,None,:] - ratoms[None,:,:])**2).sum(-1))
    choques = int((d < 2.0).sum())
    cterm = pts[-3:]  # backbone del ultimo residuo
    d_ecd = np.sqrt(((cterm[:,None,:] - ecd_atoms[None,:,:])**2).sum(-1)).min()
    return choques, d_ecd

best = (1e9, 1e9, None)
for k in range(4000):
    anclaje = c_bol + np.random.uniform(-4, 4, 3)
    tilt = math.radians(random.uniform(0, 40))
    axis_rand = np.random.normal(size=3)
    e = e0*math.cos(tilt) + axis_rand/np.linalg.norm(axis_rand)*math.sin(tilt)
    e = e/np.linalg.norm(e)
    roll = random.uniform(0, 2*math.pi)
    pts = helice(anclaje, e, roll)
    ch, d_ecd = score(pts)
    if (ch, d_ecd) < best[:2]: best = (ch, d_ecd, (anclaje, e, roll))
    if k % 2000 == 1999: print(f'  {k+1} muestras: mejor choques={best[0]}, C-ECD={best[1]:.1f} A')

print(f'\nMEJOR: choques={best[0]}, distancia extremo C al ECD={best[1]:.1f} A')
if best[0] == 0:
    print('VIABLE: existe pose de helice ideal sin choques con anclaje en el bolsillo (solo coordenadas del target).')
else:
    print('NO demostrado con este muestreo: aumentar muestras o permitir kink (dos segmentos).')
