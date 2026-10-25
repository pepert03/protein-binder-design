#!/usr/bin/env python3
# verifica_scaffold_v2.py - Verificador PURO PYTHON (sin numpy) para Umbrel.
# Reconstruye la pose de scaffold_v2 desde 10 numeros (anclaje, eje, giro) y comprueba:
#   0 choques backbone(<2.0 A) vs receptor 6X18, N-term en bolsillo (U367-391), C-term en ECD.
# Uso: python3 verifica_scaffold_v2.py [6X18_RP.pdb]
import math, sys

ANCLAJE = (110.96145509357878, 131.20303604138286, 117.79731715151486)
EJE = (0.12458546295067466, 0.7641837728928252, -0.6328519761118336)
ROLL = 4.7813817536481205
NRES = 36

def cargar_R(path):
    R = {}
    for line in open(path):
        if line.startswith('ATOM') and line[21] == 'R' and line[17:20].strip() != 'HOH':
            rs = int(line[22:26])
            R.setdefault(rs, []).append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return R

def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a):
    n = math.sqrt(sum(x*x for x in a)); return (a[0]/n, a[1]/n, a[2]/n)
def add(a, b, s=1.0):
    return (a[0]+b[0]*s, a[1]+b[1]*s, a[2]+b[2]*s)

u = norm(cross(EJE, (1,0,0)))
v = cross(EJE, u)
helix = []
for i in range(NRES):
    t = i*1.5; ang = math.radians(100*i) + ROLL
    ca = add(add(ANCLAJE, EJE, t), add(u, v, 0), 0)
    off = (2.3*math.cos(ang)*u[0] + 2.3*math.sin(ang)*v[0],
           2.3*math.cos(ang)*u[1] + 2.3*math.sin(ang)*v[1],
           2.3*math.cos(ang)*u[2] + 2.3*math.sin(ang)*v[2])
    ca = add(add(ANCLAJE, EJE, t), off)
    helix.append([add(ca, EJE, -1.32), ca, add(ca, EJE, 1.52), add(add(ca, EJE, 1.52), v, 0.8)])

R = cargar_R(sys.argv[1] if len(sys.argv) > 1 else '6X18_RP.pdb')
ratoms = [q for r in R for q in R[r]]
bol = [q for r in range(367, 392) if r in R for q in R[r]]
ecd = [q for r in range(29, 140) if r in R for q in R[r]]

choques = sum(1 for p in (p for res in helix for p in res) for q in ratoms if math.dist(p, q) < 2.0)
n_bols = sum(1 for i in range(12) if min(math.dist(p, q) for p in helix[i] for q in bol) <= 6.0)
d_ecd = min(math.dist(p, q) for p in helix[-1] for q in ecd)
print(f'choques={choques} (esperado 0) | N-bolsillo={n_bols} (esperado 12) | C-ECD={d_ecd:.1f} A (esperado 2.5)')
print('REPRODUCIDO' if choques == 0 and n_bols == 12 and abs(d_ecd - 2.5) < 0.5 else 'NO REPRODUCIDO')
