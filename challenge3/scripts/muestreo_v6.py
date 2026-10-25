#!/usr/bin/env python3
# muestreo_v6.py - Busqueda Monte Carlo v6: gate geometrico v5 + PACKING poly-ALA EN EL BUCLE.
# Proxy de secuencia DECLARADO Y FIJO: poly-ALA (36 x A) = sidechain minima; la pose pasa solo si
# incluso sidechains minimas caben. NO es optimizacion de secuencia. Receptor: 6X18 cadena R real,
# rigido, contexto fijo en FASPR (minusculas). Choques: solape vdW >0.4 A vs coords receptor originales.
# Control unico conocido: nativo 6X18 P-R = 3 choques, dmin 2.46 A. Umbral <=5 EXPLORATORIO (no validado).
# Determinista: seeds fijas. Registra: muestras, descartes por etapa, histograma, objetivo.
import math, random, subprocess, os, sys
# Rutas por entorno (RP_PDB/V6_DIR/FASPR_DIR) - parche del verificador 25 oct
import numpy as np

SEED = 20261025
random.seed(SEED); np.random.seed(SEED)
VDW = {'N':1.55, 'C':1.70, 'O':1.52, 'S':1.80}
def elem(n): return n.strip()[0] if n.strip()[0] in VDW else 'C'

# ---- receptor: 6X18_R, solo residuos backbone-completo (para FASPR y choques) ----
RP = os.environ.get('RP_PDB', '/tmp/ch3pack/challenge3/target/6X18_RP.pdb')
byres = {}
for line in open(RP):
    if line.startswith('ATOM') and line[21] == 'R' and line[16] in (' ', 'A'):
        byres.setdefault(int(line[22:26]), {})[line[12:16].strip()] = (line[17:20].strip(), [float(line[30:38]), float(line[38:46]), float(line[46:54])])
completos = sorted(r for r, s in byres.items() if {'N','CA','C','O'} <= set(s))
TRES2UNO = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
rec_seq = ''.join(TRES2UNO[byres[r]['CA'][0]] for r in completos)
rec_atoms = [(r, nm, xyz) for r in completos for nm, (rn, xyz) in byres[r].items()]
rec_xyz = np.array([a[2] for a in rec_atoms]); rec_r = np.array([VDW[elem(a[1])] for a in rec_atoms])
print(f"receptor: {len(completos)} residuos completos de {len(byres)}; {len(rec_atoms)} atomos")

# ---- geometria de helice ideal verificada (misma que v5) ----
B_N_CA, B_CA_C, B_C_N, B_C_O = 1.46, 1.52, 1.33, 1.23
A_NCA_C = math.radians(111.0); A_CA_C_N = math.radians(116.0); A_C_N_CA = math.radians(121.0)
PHI, PSI, OMEGA = math.radians(-57.8-180), math.radians(-47.0-180), math.radians(0.0)
def nerf(a, b, c, l, theta, chi):
    bc = (c-b)/np.linalg.norm(c-b); ab = a-b
    n = np.cross(ab, bc); n /= np.linalg.norm(n); m = np.cross(n, bc)
    return c + (-l*math.cos(theta)*bc + l*math.sin(theta)*math.cos(chi)*m + l*math.sin(theta)*math.sin(chi)*n)
def build_helix(n=38):
    N = np.array([0.0,1.46,0.0]); CA = np.array([0.0,0.0,0.0]); C = np.array([1.52,0.0,0.0])
    out = []
    for i in range(n):
        out.append({'N':N.copy(),'CA':CA.copy(),'C':C.copy(),'O':None})
        if i < n-1:
            Nn = nerf(N,CA,C,B_C_N,A_CA_C_N,PSI); CAn = nerf(CA,C,Nn,B_N_CA,A_C_N_CA,OMEGA); Cn = nerf(C,Nn,CAn,B_CA_C,A_NCA_C,PHI)
            N,CA,C = Nn,CAn,Cn
    for i,r in enumerate(out):
        if i < n-1:
            d1 = (r['CA']-r['C'])/np.linalg.norm(r['CA']-r['C']); d2 = (out[i+1]['N']-r['C'])/np.linalg.norm(out[i+1]['N']-r['C'])
        else:
            d1 = (r['CA']-r['C'])/np.linalg.norm(r['CA']-r['C']); d2 = d1
        r['O'] = r['C'] - B_C_O*(d1+d2)/np.linalg.norm(d1+d2)
    return out
HELIX = build_helix()[1:37]
CA_CANON = np.array([r['CA'] for r in HELIX]); ALL_CANON = np.array([r[a] for r in HELIX for a in ('N','CA','C','O')])
def kabsch(P,Q):
    Pc,Qc = P-P.mean(0), Q-Q.mean(0)
    V,S,Wt = np.linalg.svd(Pc.T@Qc); D = np.diag([1,1,np.sign(np.linalg.det(V@Wt))])
    return V@D@Wt
def trace_ca(anclaje,e,roll,n=36):
    u = np.cross(e,[1,0,0])
    if np.linalg.norm(u)<1e-6: u = np.cross(e,[0,1,0])
    u/=np.linalg.norm(u); v = np.cross(e,u)
    return np.array([anclaje+e*(i*1.5)+2.3*(math.cos(math.radians(100*i)+roll)*u+math.sin(math.radians(100*i)+roll)*v) for i in range(n)])
def pose_atoms(anclaje,e,roll):
    tr = trace_ca(anclaje,e,roll); Rk = kabsch(CA_CANON,tr)
    return (ALL_CANON-CA_CANON.mean(0))@Rk + tr.mean(0)

# ---- zonas (mismo criterio que v5) ----
R = {r: {nm: xyz for nm,(rn,xyz) in byres[r].items()} for r in completos}
def centroid(residues):
    return np.array([R[r]['CA'] for r in residues if r in R]).mean(0)
c_bol = centroid(range(367,392)); c_ecd = centroid(range(120,139))
ecd_atoms = np.array([q for r in range(29,140) if r in R for q in R[r].values()])
bol_atoms = np.array([q for r in range(367,392) if r in R for q in R[r].values()])
e0 = c_ecd-c_bol; e0 = e0/np.linalg.norm(e0)

def gate_geo(pts):
    d = np.sqrt(((pts[:,None,:]-rec_xyz[None,:,:])**2).sum(-1))
    ch = int((d<2.0).sum())
    dn = np.sqrt(((pts[:48][:,None,:]-bol_atoms[None,:,:])**2).sum(-1))
    n_bols = int(sum(1 for i in range(12) if dn[i*4:(i+1)*4].min() <= 6.0))
    d_ecd = float(np.sqrt(((pts[-4:][:,None,:]-ecd_atoms[None,:,:])**2).sum(-1)).min())
    return ch, n_bols, d_ecd

# ---- packing poly-ALA en contexto ----
VDIR = os.environ.get('V6_DIR', '/tmp/pose/v6'); FASPR_DIR = os.environ.get('FASPR_DIR', '/tmp/pose/FASPR')
os.makedirs(VDIR, exist_ok=True)
def fmt(serial, name, resname, chain, resseq, x, y, z):
    an = (' '+name).ljust(4) if len(name) < 4 else name
    return f"ATOM  {serial:>5} {an}{resname:>4} {chain}{resseq:>4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {name[0]:>2}"
def escribe_complejo(pts):
    L = []; s = 1
    for r in completos:
        for nm,(rn,xyz) in byres[r].items():
            L.append(fmt(s, nm, rn, 'A', r, *xyz)); s += 1
    for i in range(36):
        for j, nm in enumerate(('N','CA','C','O')):
            x,y,z = pts[i*4+j]
            L.append(fmt(s, nm, 'ALA', 'B', i+1, x, y, z)); s += 1
    open(VDIR + '/tmp.pdb','w').write('\n'.join(L)+'\nEND\n')
def pack_choques(pts):
    escribe_complejo(pts)
    open(VDIR + '/seq.txt','w').write(rec_seq.lower() + 'A'*36 + '\n')
    r = subprocess.run(['./faspr','-i',VDIR + '/tmp.pdb','-o',VDIR + '/tmp_out.pdb','-s',VDIR + '/seq.txt'],
                       cwd=FASPR_DIR, capture_output=True)
    if not os.path.exists(VDIR + '/tmp_out.pdb'): return None
    seen_b, bxyz, br = set(), [], []
    for line in open(VDIR + '/tmp_out.pdb'):
        if line.startswith('ATOM') and line[21] == 'B':
            k = (int(line[22:26]), line[12:16].strip())
            if k not in seen_b:
                seen_b.add(k); bxyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])]); br.append(VDW[elem(line[12:16])])
    b_xyz = np.array(bxyz); b_r = np.array(br)
    d = np.sqrt(((b_xyz[:,None,:]-rec_xyz[None,:,:])**2).sum(-1))
    return int((d < (b_r[:,None]+rec_r[None,:]-0.4)).sum())

# ---- muestreo global + refino, packing en supervivientes geometricos ----
REG = {'muestras':0, 'descarte_choque_bb':0, 'descarte_nbols':0, 'descarte_ecd':0, 'supervivientes_geo':0, 'empaquetadas':0}
hist = []
def evalua(anclaje, e, roll, empacotar=True):
    pts = pose_atoms(anclaje, e, roll)
    ch, n_bols, d_ecd = gate_geo(pts)
    REG['muestras'] += 1
    if ch > 0: REG['descarte_choque_bb'] += 1; return None
    if n_bols < 3: REG['descarte_nbols'] += 1; return None
    if d_ecd > 6.0: REG['descarte_ecd'] += 1; return None
    REG['supervivientes_geo'] += 1
    if not empacotar: return (999, d_ecd, (anclaje.copy(), e.copy(), roll))
    ncl = pack_choques(pts)
    REG['empaquetadas'] += 1
    if ncl is None: return None
    hist.append(ncl)
    return (ncl, d_ecd, (anclaje.copy(), e.copy(), roll))


stage = sys.argv[1] if len(sys.argv) > 1 else 'global'
best = (10**9, 1e9, None)
if stage == 'global':
    # colecta laxa: guarda poses con 0 choques bb y n_bols>=2, d_ecd<=10; ranking geometrico; packing despues
    cands = []
    N_GLOBAL = 4400
    for k in range(N_GLOBAL):
        anclaje = c_bol + np.random.uniform(-4,4,3)
        tilt = math.radians(random.uniform(0,40)); ar = np.random.normal(size=3)
        e = e0*math.cos(tilt) + ar/np.linalg.norm(ar)*math.sin(tilt); e /= np.linalg.norm(e)
        roll = random.uniform(0, 2*math.pi)
        pts = pose_atoms(anclaje, e, roll)
        ch, n_bols, d_ecd = gate_geo(pts)
        REG['muestras'] += 1
        if ch > 0: REG['descarte_choque_bb'] += 1; continue
        if n_bols < 2 or d_ecd > 10.0: REG['descarte_nbols'] += 1; continue
        cands.append(((-n_bols, d_ecd), anclaje.copy(), e.copy(), roll))
    cands.sort(key=lambda c: c[0])
    np.save(VDIR + '/cands.npy', np.array([[*c[1], *c[2], c[3]] for c in cands])) if cands else None
    with open(VDIR + '/cands_meta.txt','w') as f:
        f.write(f"{len(cands)} candidatos geometricos laxos de {N_GLOBAL}\n")
    print(f"GLOBAL laxo: {len(cands)} candidatos de {N_GLOBAL} (0 choques bb, n_bols>=2, d_ecd<=10)")
    print(f"registro: {REG}")
elif stage == 'pack':
    # packing poly-ALA de los N mejores geometricos
    NPACK = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    V = np.load(VDIR + '/cands.npy')
    res = []
    for i, v in enumerate(V[:NPACK]):
        pts = pose_atoms(v[:3], v[3:6], float(v[6]))
        ncl = pack_choques(pts)
        _, n_bols, d_ecd = gate_geo(pts)
        if ncl is not None:
            res.append((ncl, d_ecd, -n_bols, i))
            hist.append(ncl)
        if i % 50 == 49: print(f"  pack {i+1}/{min(NPACK,len(V))}", flush=True)
    res.sort()
    np.save(VDIR + '/pack_rank.npy', np.array([r[:3] for r in res]))
    with open(VDIR + '/pack_rank_idx.txt','w') as f:
        for r in res: f.write(f"{r[3]} {r[0]} {r[1]:.2f} {r[2]}\n")
    h = np.array(hist)
    print(f"PACK: n={len(res)} | choques min={h.min()} p25={np.percentile(h,25):.0f} med={np.median(h):.0f} max={h.max()} | <=5: {(h<=5).sum()}")
    print(f"mejores 5 (choques, d_ecd, -n_bols, idx): {res[:5]}")

elif stage == 'refine_v5':
    # refinamiento local ALREDEDOR de la pose v5 (mi scaffold de novo previo, target-only)
    a0 = np.array([113.337, 138.316, 109.455]); e0b = np.array([-0.0613, 0.8296, -0.5550]); r0 = 1.8819
    best = (10**9, 1e9, None)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 800
    for k in range(N):
        anclaje = a0 + np.random.uniform(-3,3,3)
        tilt = math.radians(random.uniform(-20,20)); ar = np.random.normal(size=3)
        e = e0b*math.cos(tilt)+ar/np.linalg.norm(ar)*math.sin(tilt); e/=np.linalg.norm(e)
        roll = r0 + random.uniform(-0.7,0.7)
        pts = pose_atoms(anclaje, e, roll)
        ch, n_bols, d_ecd = gate_geo(pts)
        REG['muestras'] += 1
        if ch > 0: REG['descarte_choque_bb'] += 1; continue
        if n_bols < 3 or d_ecd > 6.0: REG['descarte_nbols'] += 1; continue
        REG['supervivientes_geo'] += 1
        ncl = pack_choques(pts)
        REG['empaquetadas'] += 1
        if ncl is None: continue
        hist.append(ncl)
        if (ncl, d_ecd) < (best[0], best[1]):
            best = (ncl, d_ecd, (anclaje.copy(), e.copy(), roll))
        if k % 100 == 99: print(f"  refino_v5 {k+1}/{N}: mejor={best[0]}, superv={REG['supervivientes_geo']}", flush=True)
    print(f"REFINO_V5: muestras={N} superv={REG['supervivientes_geo']} mejor packing={best[0]}, C-ECD={best[1]:.2f}")
    if best[2] is not None:
        print(f"anclaje={np.round(best[2][0],3).tolist()} eje={np.round(best[2][1],4).tolist()} roll={best[2][2]:.4f}")
        np.save(VDIR + '/best_final.npy', np.concatenate([best[2][0], best[2][1], [best[2][2]]]))
        with open(VDIR + '/best_final.txt','w') as f: f.write(f"{best[0]} {best[1]}")
    print(f"registro: {REG}")
    if len(hist):
        h = np.array(hist); print(f"choques refino_v5: min={h.min()} med={np.median(h):.0f} max={h.max()} | <=5: {(h<=5).sum()}/{len(h)}")

elif stage == 'refine_wide':
    a0 = np.array([111.93, 135.843, 111.966]); e0b = np.array([0.1114, 0.7591, -0.6413]); r0 = 1.8939
    best = (10**9, 1e9, None)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
    for k in range(N):
        anclaje = a0 + np.random.uniform(-5,5,3)
        tilt = math.radians(random.uniform(-30,30)); ar = np.random.normal(size=3)
        e = e0b*math.cos(tilt)+ar/np.linalg.norm(ar)*math.sin(tilt); e/=np.linalg.norm(e)
        roll = r0 + random.uniform(-1.0,1.0)
        pts = pose_atoms(anclaje, e, roll)
        ch, n_bols, d_ecd = gate_geo(pts)
        REG['muestras'] += 1
        if ch > 0: REG['descarte_choque_bb'] += 1; continue
        if n_bols < 2 or d_ecd > 8.0: REG['descarte_nbols'] += 1; continue
        REG['supervivientes_geo'] += 1
        ncl = pack_choques(pts)
        REG['empaquetadas'] += 1
        if ncl is None: continue
        hist.append(ncl)
        if (ncl, d_ecd) < (best[0], best[1]):
            best = (ncl, d_ecd, (anclaje.copy(), e.copy(), roll))
        if k % 250 == 249: print(f"  wide {k+1}/{N}: mejor={best[0]}, superv={REG['supervivientes_geo']}", flush=True)
    print(f"WIDE: muestras={N} superv={REG['supervivientes_geo']} mejor packing={best[0]}, C-ECD={best[1]:.2f}")
    if best[2] is not None:
        print(f"anclaje={np.round(best[2][0],3).tolist()} eje={np.round(best[2][1],4).tolist()} roll={best[2][2]:.4f}")
        np.save(VDIR + '/best_wide.npy', np.concatenate([best[2][0], best[2][1], [best[2][2]]]))
        with open(VDIR + '/best_wide.txt','w') as f: f.write(f"{best[0]} {best[1]}")
    print(f"registro: {REG}")
    if len(hist):
        h = np.array(hist); print(f"choques wide: min={h.min()} p25={np.percentile(h,25):.0f} med={np.median(h):.0f} max={h.max()} | <=5: {(h<=5).sum()}/{len(h)}")

elif stage == 'refine2':
    v = np.load(VDIR + '/best_wide.npy')
    a0, e0b, r0 = v[:3], v[3:6], float(v[6])
    best = (10**9, 1e9, None)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 1200
    for k in range(N):
        anclaje = a0 + np.random.uniform(-1.5,1.5,3)
        tilt = math.radians(random.uniform(-10,10)); ar = np.random.normal(size=3)
        e = e0b*math.cos(tilt)+ar/np.linalg.norm(ar)*math.sin(tilt); e/=np.linalg.norm(e)
        roll = r0 + random.uniform(-0.35,0.35)
        pts = pose_atoms(anclaje, e, roll)
        ch, n_bols, d_ecd = gate_geo(pts)
        REG['muestras'] += 1
        if ch > 0: REG['descarte_choque_bb'] += 1; continue
        if n_bols < 2 or d_ecd > 8.0: REG['descarte_nbols'] += 1; continue
        REG['supervivientes_geo'] += 1
        ncl = pack_choques(pts)
        REG['empaquetadas'] += 1
        if ncl is None: continue
        hist.append(ncl)
        if (ncl, d_ecd) < (best[0], best[1]):
            best = (ncl, d_ecd, (anclaje.copy(), e.copy(), roll))
        if k % 300 == 299: print(f"  refine2 {k+1}/{N}: mejor={best[0]}, superv={REG['supervivientes_geo']}", flush=True)
    print(f"REFINE2: muestras={N} superv={REG['supervivientes_geo']} mejor packing={best[0]}, C-ECD={best[1]:.2f}")
    if best[2] is not None:
        print(f"anclaje={np.round(best[2][0],3).tolist()} eje={np.round(best[2][1],4).tolist()} roll={best[2][2]:.4f}")
        np.save(VDIR + '/best_final.npy', np.concatenate([best[2][0], best[2][1], [best[2][2]]]))
        with open(VDIR + '/best_final.txt','w') as f: f.write(f"{best[0]} {best[1]}")
    print(f"registro: {REG}")
    if len(hist):
        h = np.array(hist); print(f"choques refine2: min={h.min()} med={np.median(h):.0f} max={h.max()} | <=5: {(h<=5).sum()}/{len(h)}")
elif stage == 'refine':
    V = np.load(VDIR + '/cands.npy')
    best_idx = int(open(VDIR + '/pack_rank_idx.txt').readline().split()[0])
    v = V[best_idx]
    a0, e0b, r0 = v[:3], v[3:6], float(v[6])
    bn = int(open(VDIR + '/pack_rank_idx.txt').readline().split()[1])
    best = (bn, 1e9, None)
    for k in range(1200):
        anclaje = a0 + np.random.uniform(-1.5,1.5,3)
        tilt = math.radians(random.uniform(-8,8)); ar = np.random.normal(size=3)
        e = e0b*math.cos(tilt)+ar/np.linalg.norm(ar)*math.sin(tilt); e/=np.linalg.norm(e)
        roll = r0 + random.uniform(-0.35,0.35)
        pts = pose_atoms(anclaje, e, roll)
        ch, n_bols, d_ecd = gate_geo(pts)
        REG['muestras'] += 1
        if ch > 0: REG['descarte_choque_bb'] += 1; continue
        ncl = pack_choques(pts)
        REG['empaquetadas'] += 1
        if ncl is None: continue
        hist.append(ncl)
        if (ncl, d_ecd) < (best[0], best[1]):
            best = (ncl, d_ecd, (anclaje.copy(), e.copy(), roll)); a0, e0b, r0 = best[2]
        if k % 200 == 199: print(f"  refino {k+1}/1200: mejor={best[0]}, C-ECD={best[1]:.2f}", flush=True)
    print(f"REFINO: mejor packing={best[0]}, C-ECD={best[1]:.2f}")
    print(f"anclaje={np.round(best[2][0],3).tolist()} eje={np.round(best[2][1],4).tolist()} roll={best[2][2]:.4f}")
    np.save(VDIR + '/best_final.npy', np.concatenate([best[2][0], best[2][1], [best[2][2]]]))
    with open(VDIR + '/best_final.txt','w') as f: f.write(f"{best[0]} {best[1]}")
    print(f"registro: {REG}")
    h = np.array(hist); print(f"choques refino: min={h.min()} med={np.median(h):.0f} max={h.max()} | <=5: {(h<=5).sum()}/{len(h)}")
