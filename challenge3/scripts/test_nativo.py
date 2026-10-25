import math, sys
import numpy as np

def cargar_R(path):
    R = {}
    for line in open(path):
        if line.startswith('ATOM') and line[21] == 'R' and line[17:20].strip() != 'HOH':
            rs = int(line[22:26]); an = line[12:16].strip()
            R.setdefault(rs, {})[an] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return R
R = cargar_R(sys.argv[1] if len(sys.argv) > 1 else '6X18_RP.pdb')
ratoms = np.array([q for r in R for q in R[r].values()])

B_N_CA, B_CA_C, B_C_N, B_C_O = 1.46, 1.52, 1.33, 1.23
A_NCA_C = math.radians(111.0); A_CA_C_N = math.radians(116.0); A_C_N_CA = math.radians(121.0)
PHI, PSI, OMEGA = math.radians(-57.8-180), math.radians(-47.0-180), math.radians(0.0)

def nerf(a, b, c, l, theta, chi):
    bc = (c - b) / np.linalg.norm(c - b); ab = a - b
    n = np.cross(ab, bc); n /= np.linalg.norm(n); m = np.cross(n, bc)
    return c + (-l*math.cos(theta)*bc + l*math.sin(theta)*math.cos(chi)*m + l*math.sin(theta)*math.sin(chi)*n)

def build_helix(n=36):
    N = np.array([0.0, 1.46, 0.0]); CA = np.array([0.0, 0.0, 0.0]); C = np.array([1.52, 0.0, 0.0])
    out = []
    for i in range(n):
        out.append({'N': N.copy(), 'CA': CA.copy(), 'C': C.copy()})
        if i < n - 1:
            Nn = nerf(N, CA, C, B_C_N, A_CA_C_N, PSI)
            CAn = nerf(CA, C, Nn, B_N_CA, A_C_N_CA, OMEGA)
            Cn = nerf(C, Nn, CAn, B_CA_C, A_NCA_C, PHI)
            N, CA, C = Nn, CAn, Cn
    for i, r in enumerate(out):
        d1 = (r['CA'] - r['C']) / np.linalg.norm(r['CA'] - r['C'])
        d2 = ((out[i+1]['N'] - r['C']) / np.linalg.norm(out[i+1]['N'] - r['C'])) if i < n-1 else d1
        r['O'] = r['C'] + B_C_O * (d1 + d2) / np.linalg.norm(d1 + d2)
    return out

HELIX = build_helix(30)
CA_CANON = np.array([r['CA'] for r in HELIX])
ALL_CANON = np.array([r[a] for r in HELIX for a in ('N', 'CA', 'C', 'O')])

def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt

P = {}
for line in open(sys.argv[1] if len(sys.argv) > 1 else '6X18_RP.pdb'):
    if line.startswith('ATOM') and line[21] == 'P' and line[17:20].strip() != 'HOH':
        rs = int(line[22:26]); an = line[12:16].strip()
        P.setdefault(rs, {})[an] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
nat_ca = np.array([P[r]['CA'] for r in sorted(P) if 'CA' in P[r]])
print('CA nativos:', len(nat_ca))
Rk = kabsch(CA_CANON, nat_ca)
placed = (ALL_CANON - CA_CANON.mean(0)) @ Rk + nat_ca.mean(0)
fitted = (CA_CANON - CA_CANON.mean(0)) @ Rk + nat_ca.mean(0)
rmsd = float(np.sqrt(((fitted - nat_ca)**2).sum() / len(nat_ca)))
d = np.sqrt(((placed[:, None, :] - ratoms[None, :, :]) ** 2).sum(-1))
print(f'RMSD ajuste: {rmsd:.2f} A')
print(f'choques backbone real sobre traza nativa: {int((d < 2.0).sum())}; min {d.min():.2f} A')
