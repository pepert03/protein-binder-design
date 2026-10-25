#!/usr/bin/env python3
# viabilidad_scaffold.py - Chequeo geometrico de viabilidad de scaffolding por extremos
# usando SOLO coordenadas del target 6X18 (R = receptor, P = peptido nativo como referencia
# de corredor; NO se copia su secuencia). REGLAS OFICIALES (de novo/zero-shot): el peptido
# nativo sirve SOLO para calibrar geometria/filtros, NUNCA como base de diseno.
# Uso: python3 viabilidad_scaffold.py <6X18_RP.pdb> Responde: (a) que luz tiene el corredor entre el
# ECD y el bolsillo de activacion, (b) cuantos residuos hacen falta para cruzarlo en helice,
# (c) si una varilla helicoidal en el marco nativo tendria choques de backbone.
import math, sys

def cargar(path):
    A, B = {}, {}
    for line in open(path):
        if line.startswith('ATOM') and line[17:20].strip() != 'HOH':
            ch = line[21]; rs = int(line[22:26]); an = line[12:16].strip()
            d = A if ch == 'R' else (B if ch == 'P' else None)
            if d is not None:
                d.setdefault(rs, {})[an] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return A, B

R, P = cargar(sys.argv[1] if len(sys.argv) > 1 else '6X18_RP.pdb')
pep_res = sorted(P)
print(f'Peptido nativo: resseq {pep_res[0]}-{pep_res[-1]} ({len(pep_res)} res)')

def dmin(p, q):
    return min(math.dist(a, b) for a in p for b in q)

# (a) Luz del corredor: min distancia de cada CA del peptido al receptor (atomos pesados)
print('\nLuz del corredor (CA peptido -> receptor, min A):')
apretados = 0
for r in pep_res:
    if 'CA' not in P[r]: continue
    ca = [P[r]['CA']]
    m = min(dmin(ca, list(R[t].values())) for t in R)
    flag = ' <2A!' if m < 2.0 else ''
    if m < 2.0: apretados += 1
    if r <= 20 or m < 2.0: print(f'  P{r}: {m:.2f}{flag}')
print(f'  CA con luz <2.0 A: {apretados} de {len(pep_res)}')

# (b) Distancia ECD(U120-138) <-> bolsillo(U367-391) en el RECEPTOR (superficie a superficie)
ecd = [r for r in R if 120 <= r <= 138]; bol = [r for r in R if 367 <= r <= 391]
d_eb = min(dmin(list(R[a].values()), list(R[b].values())) for a in ecd for b in bol)
print(f'\nDistancia minima superficie ECD(120-138) <-> bolsillo(367-391): {d_eb:.1f} A')

# Span del extremo del peptido nativo: CA de P7 (N, insertado) a CA del ultimo residuo
caN, caC = P[pep_res[0]]['CA'], P[pep_res[-1]]['CA']
span = math.dist(caN, caC)
print(f'Span CA N-term(P{pep_res[0]}) -> C-term(P{pep_res[-1]}): {span:.1f} A en {len(pep_res)} res')
# Helice ideal: 1.5 A/residuo de avance axial
import math as m2
need = span / 1.5
print(f'Una helice ideal que cubra ese span necesita ~{m2.ceil(need)} res (tenemos {len(pep_res)}: holgado={len(pep_res) >= need})')

# (c) Choques de backbone del marco nativo: N/CA/C/O del peptido vs receptor <2.0 A
bb = ('N', 'CA', 'C', 'O')
nch = sum(1 for r in pep_res for an in bb if an in P[r] for t in R for q in R[t].values() if math.dist(P[r][an], q) < 2.0)
print(f'Choques backbone(<2.0 A) del marco nativo vs receptor: {nch}')
