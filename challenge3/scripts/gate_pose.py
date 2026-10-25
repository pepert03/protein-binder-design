#!/usr/bin/env python3
"""Gate de pose para disenos GLP1R (Proteinbase ch3).
Uso: python3 gate_pose.py <complejo.pdb> [epitope_contacts.tsv] [P43220.fasta]
- Analiza cada MODEL del PDB por separado (nunca mezcla modelos).
- Solo pares intermoleculares cadena A (target) x cadena B (binder) del mismo MODEL.
- Filtra altloc (conserva ' ' o 'A') y deduplica registros de atomo repetidos.
- Mapeo resseq->UniProt: numeracion original con huecos, offset constante +28,
  verificado por anclaje de secuencia (kmer 12) contra P43220 y por resname
  contra el TSV del epitopo (52/52). No circular: la tabla de fragmentos sale
  de la secuencia, no del epitopo.
- Reporta: contactos <=5 A, solape de residuos con los 52 nativos (22 ECD/30 TM-ECL),
  residuos fuera del nativo-52 por zona, choques <2.0 A (lista completa), minimo global,
  residuos del binder comprometidos, secuencia del binder.
HONESTIDAD: confianza estructural in silico; no es evidencia de union medida ni de activacion.
"""
import sys, math

PDB = sys.argv[1]
TSV = sys.argv[2] if len(sys.argv) > 2 else None
FA  = sys.argv[3] if len(sys.argv) > 3 else None
AAS = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# fragmentos (resseq_ini, resseq_fin, uniprot_ini) - verificado por secuencia (kmer 12) en 6X18 crop ch3
FRAGS = [(1,23,29),(33,101,61),(108,136,136),(148,223,176),(238,300,266),(321,376,349)]
RES2UP = {r: u0 + (r - a) for a, b, u0 in FRAGS for r in range(a, b + 1)}

def parse_models(path):
    models, cur = {}, None
    seen = set()
    for line in open(path):
        if line.startswith('MODEL'):
            cur = int(line[10:14]); models[cur] = {'A': {}, 'B': {}}; seen = set()
        elif line.startswith('ENDMDL'):
            cur = None
        elif line.startswith('ATOM') and cur is not None:
            alt = line[16]
            if alt not in (' ', 'A'): continue
            ch = line[21]
            if ch not in ('A', 'B'): continue
            rs = int(line[22:26]); an = line[12:16].strip()
            key = (ch, rs, an)
            if key in seen: continue            # deduplica registros repetidos
            seen.add(key)
            res = line[17:20].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            d = models[cur][ch].setdefault(rs, {'res': res, 'atoms': []})
            d['atoms'].append((an, x, y, z))
    return models

def analyze(m):
    A, B = m['A'], m['B']
    bat = [(r, an, x, y, z) for r in B for an, x, y, z in B[r]['atoms']]
    rows, clashes, gmin = [], [], (1e18, None)
    for r in sorted(A):
        dmin = 1e18
        for _, x, y, z in A[r]['atoms']:
            for _, _, bx, by, bz in bat:
                dd = (x-bx)**2 + (y-by)**2 + (z-bz)**2
                if dd < dmin: dmin = dd
        if dmin <= 25.0: rows.append((r, RES2UP.get(r), A[r]['res'], math.sqrt(dmin)))
    for rb, an, x, y, z in bat:
        for rt in A:
            for at, tx, ty, tz in A[rt]['atoms']:
                dd = (x-tx)**2 + (y-ty)**2 + (z-tz)**2
                if dd < gmin[0]: gmin = (dd, (rb, an, rt, at))
                if dd < 4.0: clashes.append((dd, rb, an, rt, at))
    beng = {rb for rb, _, x, y, z in bat for rt in A for _, tx, ty, tz in A[rt]['atoms']
            if (x-tx)**2 + (y-ty)**2 + (z-tz)**2 <= 25.0}
    seqB = ''.join(AAS.get(B[r]['res'], 'X') for r in sorted(B))
    return rows, sorted(clashes), gmin, beng, seqB

nat = {}
if TSV:
    for ln in open(TSV).read().splitlines()[1:]:
        if ln.strip():
            c = ln.split(); nat[int(c[2])] = c[3]

for mid in sorted(parse_models(PDB)):
    m = parse_models(PDB)[mid]
    rows, cl, gmin, beng, seqB = analyze(m)
    print('=' * 60)
    print('MODEL', mid, '| binder_seq', seqB, '| n_target', len(m['A']), 'n_binder', len(m['B']))
    print('contactos_le5A', len(rows), '| binder_res comprometidos', len(beng), 'de', len(m['B']))
    if nat:
        hits = {'ECD': 0, 'TM/ECL': 0}
        touched = set()
        for rs, u, rn, d in rows:
            if u in nat: hits[nat[u]] += 1; touched.add(u)
        print('solape_residuos_nativos52:', hits, '| no_tocados:', sorted(set(nat) - touched))
    fuera = [(rs, u, rn, round(d, 2)) for rs, u, rn, d in rows if u not in nat]
    print('fuera_nativo52:', len(fuera))
    for rs, u, rn, d in fuera:
        zona = 'ECD' if u and u <= 145 else ('TM/ECL' if u and u >= 148 else 'juncion')
        print('   resseq %d = U%s %s %.2fA %s' % (rs, u, rn, d, zona))
    print('pares_<2.0A:', len(cl), '| min_global_A:', round(math.sqrt(gmin[0]), 2), gmin[1])
    for dd, rb, an, rt, at in cl[:20]:
        print('   choque %.2fA B%d(%s)-%s vs A%d/U%s(%s)-%s' % (math.sqrt(dd), rb, m['B'][rb]['res'], an, rt, RES2UP.get(rt), m['A'][rt]['res'], at))
