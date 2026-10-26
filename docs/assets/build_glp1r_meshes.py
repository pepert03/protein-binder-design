#!/usr/bin/env python3
"""Build the GLP1R (challenge 3) meshes for docs/index.html and merge them into
docs/assets/structures.js.

Receptor: 6X18 chain R (cryo-EM, experimental), drawn as a CA-trace tube (no
fake cartoon, no surface). Binders: 6X18 chain P (native GLP-1, experimental
pose) and the three de novo candidates gA_3 / gB_2 / gB_1 in their Boltz 2
predicted poses (seed 1), each rigidly aligned onto 6X18 chain R through the
166 core-TM CAs (UniProt P43220 FT TRANSMEM segments). Everything shown is a
CA trace. Candidate poses are model predictions and are NOT reproducible
between seeds (7-12 A CA RMSD); nothing here is experimentally validated.

Stdlib + numpy. Run from the repo root: python3 docs/assets/build_glp1r_meshes.py
"""
import json, base64, struct, re
import numpy as np

ROOT = '.'
EM = ROOT + '/challenge3/target/6X18_RP.pdb'
CANDS = {
    'gA_3': ('0.93', ROOT + '/challenge3/results/zeroshot_409_v2/filter_outputs/out/gA_3/boltz_results_f409_gA_3/predictions/f409_gA_3/f409_gA_3_model_0.cif'),
    'gB_2': ('0.92', ROOT + '/challenge3/results/zeroshot_409_v2/filter_outputs/out/gB_2/boltz_results_f409_gB_2/predictions/f409_gB_2/f409_gB_2_model_0.cif'),
    'gB_1': ('0.83', ROOT + '/challenge3/results/zeroshot_409_v2/filter_outputs/out/gB_1/boltz_results_f409_gB_1/predictions/f409_gB_1/f409_gB_1_model_0.cif'),
}
TM = []
for a, b in [(140,164),(176,201),(228,251),(266,290),(306,328),(349,370),(384,404)]:
    TM += list(range(a, b+1))

def pdb_ca(path, chain):
    d = {}
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA' and line[21] == chain:
            rn = int(line[22:26])
            if chain == 'R' and not (24 <= rn <= 463):
                continue  # drop waters numbered 501-520
            d[rn] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return d

def cif_ca(path, chain, off=0):
    d = {}
    for line in open(path):
        if not line.startswith('ATOM'):
            continue
        f = line.split()
        if f[3] == 'CA' and f[9] == chain:
            d[int(f[6]) + off] = np.array([float(f[10]), float(f[11]), float(f[12])])
    return d

def kabsch(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Pc.T @ Qc; U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, Q.mean(0) - (R @ P.mean(0))

def tube(points, radius, sides=8, break_dist=5.0):
    """CA-trace tube. Breaks faces across gaps > break_dist. Returns (V (n,3), N (n,3), F (m,3))."""
    n = len(points)
    V = np.zeros((n * sides, 3)); N = np.zeros((n * sides, 3))
    for i, p in enumerate(points):
        p_prev = points[i - 1] if i > 0 else p
        p_next = points[i + 1] if i < n - 1 else p
        t = p_next - p_prev
        t = t / (np.linalg.norm(t) + 1e-9)
        ref = np.array([0.0, 0.0, 1.0]) if abs(t[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        u = np.cross(t, ref); u /= np.linalg.norm(u)
        v = np.cross(t, u)
        for k in range(sides):
            ang = 2 * np.pi * k / sides
            nrm = np.cos(ang) * u + np.sin(ang) * v
            V[i * sides + k] = p + radius * nrm
            N[i * sides + k] = nrm
    F = []
    for i in range(n - 1):
        if np.linalg.norm(points[i + 1] - points[i]) > break_dist:
            continue  # chain gap: no false connection
        for k in range(sides):
            a = i * sides + k; b = i * sides + (k + 1) % sides
            c = (i + 1) * sides + k; d = (i + 1) * sides + (k + 1) % sides
            F.append((a, c, b)); F.append((b, c, d))
    return V, N, np.array(F, dtype=np.uint16)

def pack_mesh(V, N, F):
    lo = V.min(0)
    span = float((V.max(0) - lo).max())
    scale = span / 32000.0
    q = np.round((V - lo) / scale).astype('<i2')
    nn = np.round(N * 120).clip(-128, 127).astype('<i1')
    blob = q.tobytes() + nn.tobytes()
    pad = (-len(blob)) % 4
    blob += b'\0' * pad
    blob += np.ascontiguousarray(F.astype('<u2')).tobytes()
    return {'b64': base64.b64encode(blob).decode(), 'nv': int(len(V)), 'nf': int(len(F)),
            'lo': [float(x) for x in lo], 'span': span, 'has_col': False, 'i16': True}

def main():
    emR = pdb_ca(EM, 'R')
    emP = pdb_ca(EM, 'P')
    rec_keys = sorted(emR)
    # recentre: the renderer orbits and measures radius from the origin,
    # so the receptor CA mean must sit at (0,0,0) in the encoded frame
    C = np.array([emR[k] for k in rec_keys]).mean(0)
    rec_mesh = pack_mesh(*tube([emR[k] - C for k in rec_keys], radius=1.3))
    designs = []
    pep_keys = sorted(emP)
    designs.append({'name': 'GLP-1 (native)', 'tag': 'experimental pose in 6X18 (cryo-EM)',
                    'mesh': pack_mesh(*tube([emP[k] - C for k in pep_keys], radius=1.7))})
    emR_arr = {k: emR[k] for k in rec_keys}
    for name, (iptm, cif) in CANDS.items():
        pR = cif_ca(cif, 'A', off=23)
        pB = cif_ca(cif, 'B')
        common = sorted(set(TM) & set(pR) & set(emR_arr))
        assert len(common) == 166, (name, len(common))
        R, t = kabsch(np.array([pR[k] for k in common]), np.array([emR_arr[k] for k in common]))
        bk = sorted(pB)
        B = (R @ np.array([pB[k] for k in bk]).T).T + t - C
        designs.append({'name': name,
                        'tag': f'Boltz 2 prediction (seed 1), ipTM {iptm}, aligned on the 6X18 TM core; not reproducible between seeds; not validated',
                        'mesh': pack_mesh(*tube(list(B), radius=1.7))})
    # merge into structures.js
    t = open(ROOT + '/docs/assets/structures.js').read()
    m = re.search(r'export const DATA = (\{.*\});?\s*$', t, re.S)
    data = json.loads(m.group(1))
    data['glp1r'] = {'receptor': rec_mesh, 'designs': designs}
    with open(ROOT + '/docs/assets/structures.js', 'w') as f:
        f.write('export const DATA = ' + json.dumps(data) + ';\n')
    print('glp1r receptor: nv', rec_mesh['nv'], 'nf', rec_mesh['nf'])
    for d in designs:
        print('design', d['name'], 'nv', d['mesh']['nv'], 'nf', d['mesh']['nf'])

if __name__ == '__main__':
    main()
