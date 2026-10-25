#!/usr/bin/env python3
"""Construccion tetraedrica canonica de CB: validada contra CB reales de 6X18.
CB forma angulo tetraedrico (109.47 grad) con N-CA y C-CA; quiralidad decidida por datos reales."""
import math

def read_atoms(path, altok=(' ', 'A')):
    seen, out = set(), []
    for line in open(path):
        if line.startswith('ATOM') and line[16] in altok:
            k = (line[21], int(line[22:26]), line[12:16].strip())
            if k not in seen:
                seen.add(k)
                out.append((line[21], int(line[22:26]), line[12:16].strip(),
                            [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return out

def sub(a, b): return [a[i]-b[i] for i in range(3)]
def norm(v): n = math.sqrt(sum(x*x for x in v)); return [x/n for x in v]
def dot(a, b): return sum(a[i]*b[i] for i in range(3))
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def cb_tetra(N, CA, C, signo):
    u, w = norm(sub(N, CA)), norm(sub(C, CA))
    n = norm(cross(u, w))
    uw = dot(u, w)
    a = -1.0/(3.0*(1.0+uw))
    c = signo*math.sqrt(max(0.0, 1.0 - 2*a*a*(1+uw)))
    d = [a*(u[i]+w[i]) + c*n[i] for i in range(3)]
    return [CA[i] + 1.53*d[i] for i in range(3)]

def cb_mpnn(N, CA, C):
    b = sub(CA, N); c = sub(C, CA); a = cross(b, c)
    return [CA[i]-0.58273431*a[i]+0.56802827*b[i]-0.54067466*c[i] for i in range(3)]

atoms = read_atoms('scaffold_v5.pdb')
byres = {}
for ch, rs, nm, xyz in atoms:
    byres.setdefault((ch, rs), {})[nm] = xyz

# 1) quiralidad desde datos reales + error de cada construccion
err_t_plus, err_t_minus, err_m = [], [], []
for (ch, rs), d in byres.items():
    if ch == 'A' and all(k in d for k in ('N', 'CA', 'C', 'CB')):
        err_t_plus.append(math.dist(cb_tetra(d['N'], d['CA'], d['C'], +1), d['CB']))
        err_t_minus.append(math.dist(cb_tetra(d['N'], d['CA'], d['C'], -1), d['CB']))
        err_m.append(math.dist(cb_mpnn(d['N'], d['CA'], d['C']), d['CB']))
def stat(v, name):
    s = sorted(v)
    print(f"{name}: n={len(v)} med={s[len(v)//2]:.3f} max={max(v):.3f}")
stat(err_t_plus, 'tetra signo +1 vs CB real')
stat(err_t_minus, 'tetra signo -1 vs CB real')
stat(err_m, 'MPNN           vs CB real')
gana = 'signo +1' if sorted(err_t_plus)[len(err_t_plus)//2] < sorted(err_t_minus)[len(err_t_minus)//2] else 'signo -1'
print('quiralidad ganadora segun datos reales:', gana)

# 2) divergencia por posicion entre construcciones validadas, cadena B
sgn = +1 if gana == 'signo +1' else -1
receptor = [xyz for (ch, rs), d in byres.items() if ch == 'A' for nm, xyz in d.items()]
print('\npos  |d(CB_mpnn,CB_tetra)|  dmin_receptor_mpnn  dmin_receptor_tetra')
for rs in sorted(r for (c, r) in byres if c == 'B'):
    d = byres[('B', rs)]
    cm = cb_mpnn(d['N'], d['CA'], d['C']); ct = cb_tetra(d['N'], d['CA'], d['C'], sgn)
    div = math.dist(cm, ct)
    dm = min(math.dist(cm, r) for r in receptor); dt = min(math.dist(ct, r) for r in receptor)
    flag = '  <<< <2.0' if min(dm, dt) < 2.0 else ''
    print(f"B{rs:>2}     {div:5.3f} A            {dm:5.2f}          {dt:5.2f}{flag}")
