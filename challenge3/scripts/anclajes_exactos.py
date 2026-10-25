#!/usr/bin/env python3
# Metrica EXACTA reconciliada con verificador: CA del PRIMER residuo del binder vs U367-391;
# CA del ULTIMO residuo vs U120-138. Control nativo con la MISMA definicion.
# Post-relajado: mapeo por orden de cadena/fragmento (pdbfixer renumero).
import math
def load(path):
    ch, order = {}, []
    for line in open(path):
        if line.startswith('ATOM') and line[16] in (' ', 'A'):
            c = line[21]
            if c not in ch: ch[c] = []; order.append(c)
            ch[c].append({'res': int(line[22:26]), 'name': line[12:16].strip(),
                'xyz': [float(line[30:38]), float(line[38:46]), float(line[46:54])]})
    return ch, order
FRAGS = [(1,23,29),(33,101,61),(108,136,136),(148,223,176),(238,300,266),(321,376,349)]
CROP2U = {}
for a, b, u0 in FRAGS:
    for r in range(a, b+1): CROP2U[r] = u0 + (r-a)
def zonas_crop():
    bol, ecd = set(), set()
    for cr, u in CROP2U.items():
        if 367 <= u <= 391: bol.add(cr)
        if 120 <= u <= 138: ecd.add(cr)
    return bol, ecd
def dmin(A, B): return min(math.dist(a, b) for a in A for b in B)

# --- nativo (frame U directo)
ch, _ = load('nativo_RP.pdb')
caP = [a for a in ch['P'] if a['name'] == 'CA']
caR = [a for a in ch['R'] if a['name'] == 'CA']
bol = [a['xyz'] for a in caR if 367 <= a['res'] <= 391]
ecd = [a['xyz'] for a in caR if 120 <= a['res'] <= 138]
print(f"NATIVO: B1->bolsillo={dmin([caP[0]['xyz']], bol):.2f} | ultimo->ECD={dmin([caP[-1]['xyz']], ecd):.2f}")
# nativo relajado: peptido = cadena de 30 CA
ch, order = load('nativo_RP_rel.pdb')
pep = next(v for v in ch.values() if len([a for a in v if a['name']=='CA']) == 30)
rec = [a for c, v in ch.items() for a in v if v is not pep]
caP2 = [a for a in pep if a['name'] == 'CA']
caR2 = [a for a in rec if a['name'] == 'CA']  # renum: usar orden -> resseq nativo por indice
# reconstruir zonas por indice: orden de CA receptor post == orden nativo pre
caR_pre = [a for a in ch and []] # no
bol2 = [caR2[i]['xyz'] for i, a in enumerate(caR) if 367 <= a['res'] <= 391]
ecd2 = [caR2[i]['xyz'] for i, a in enumerate(caR) if 120 <= a['res'] <= 138]
print(f"NATIVO REL: B1->bolsillo={dmin([caP2[0]['xyz']], bol2):.2f} | ultimo->ECD={dmin([caP2[-1]['xyz']], ecd2):.2f}")

# --- candidata pre/post (crop): zonas por indice del frame pre
bolc, ecdc = zonas_crop()
chp, _ = load('cand2_complex.pdb')
caR_pre = [a for a in chp['A'] if a['name'] == 'CA']
caB_pre = [a for a in chp['B'] if a['name'] == 'CA']
bolp = [a['xyz'] for a in caR_pre if a['res'] in bolc]
ecdp = [a['xyz'] for a in caR_pre if a['res'] in ecdc]
print(f"CAND PRE: B1->bolsillo={dmin([caB_pre[0]['xyz']], bolp):.2f} | B36->ECD={dmin([caB_pre[-1]['xyz']], ecdp):.2f}")
chq, _ = load('cand2_complex_rel.pdb')
caR_post = [a for c in 'ABCDEF' for a in chq[c] if a['name'] == 'CA']
caB_post = [a for a in chq['G'] if a['name'] == 'CA']
bolq = [caR_post[i]['xyz'] for i, a in enumerate(caR_pre) if a['res'] in bolc]
ecdq = [caR_post[i]['xyz'] for i, a in enumerate(caR_pre) if a['res'] in ecdc]
print(f"CAND REL: B1->bolsillo={dmin([caB_post[0]['xyz']], bolq):.2f} | B36->ECD={dmin([caB_post[-1]['xyz']], ecdq):.2f}")
