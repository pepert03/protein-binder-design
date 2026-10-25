#!/usr/bin/env python3
# audita_relax.py - Auditoria con mapeo EXPLICITO secuencia->residuo->atomo entre pre y post.
# Uso: python3 audita_relax.py <pre.pdb> <post.pdb> <cadena_ligando_pre>
# Mapeo: concatena cadenas post en orden; alinea por SECUENCIA de resnames (debe ser identica
# a la de la cadena pre correspondiente, en orden); empareja atomos por (indice_residuo, nombre).
# Solo atomos pesados comunes por nombre; OXT/HTERMs (anadidos por capping) excluidos y contados.
import sys, math
import numpy as np
VDW = {'N':1.55, 'C':1.70, 'O':1.52, 'S':1.80}
def elem(n): return n.strip()[0] if n.strip()[0] in VDW else 'C'

def load(path):
    chains, order = {}, []
    for line in open(path):
        if line.startswith('ATOM') and line[16] in (' ', 'A'):
            c = line[21]
            if c not in chains: chains[c] = []; order.append(c)
            chains[c].append({'res': int(line[22:26]), 'resname': line[17:20].strip(),
                'name': line[12:16].strip(), 'elem': line[76:78].strip() or elem(line[12:16]),
                'xyz': [float(line[30:38]), float(line[38:46]), float(line[46:54])]})
    return chains, order

def residues(atoms):
    out, cur, key = [], [], None
    for a in atoms:
        k = a['res']
        if k != key: cur = []; out.append((a['resname'], cur)); key = k
        cur.append(a)
    return out

def mapea(pre_atoms, post_chains_atoms, etiqueta):
    pre_res = residues(pre_atoms)
    post_atoms = [a for atoms in post_chains_atoms for a in atoms]
    post_res = residues(post_atoms)
    assert len(pre_res) == len(post_res), f"{etiqueta}: {len(pre_res)} res pre != {len(post_res)} post"
    mism = [i for i, ((rn1, _), (rn2, _)) in enumerate(zip(pre_res, post_res)) if rn1 != rn2]
    assert not mism, f"{etiqueta}: resnames difieren en {mism[:5]}"
    mapa = []  # (i_res, {name: (pre_atom, post_atom)})
    n_oxt = 0
    for i, ((rn, r1), (_, r2)) in enumerate(zip(pre_res, post_res)):
        n_oxt += sum(1 for a in r2 if a['name'] == 'OXT')
        d1 = {a['name']: a for a in r1}; d2 = {a['name']: a for a in r2}
        comunes = {n: (d1[n], d2[n]) for n in d1 if n in d2 and d1[n]['elem'] in VDW and n != 'OXT'}
        mapa.append((i, rn, comunes))
    pares = [(p, q) for _, _, c in mapa for p, q in c.values()]
    print(f"{etiqueta}: {len(mapa)} residuos alineados por secuencia, {len(pares)} atomos pesados comunes, OXT post excluidos: {n_oxt}")
    return mapa, pares

def metricas(mapa_r, mapa_b):
    R = [p for _, _, c in mapa_r for p, q in c.items() for p in (p,)][:0]  # no usado
    rpost = [q for _, _, c in mapa_r for _, q in c.values()]
    bpost = [q for _, _, c in mapa_b for _, q in c.values()]
    rpre = [p for _, _, c in mapa_r for p, _ in c.values()]
    bpre = [p for _, _, c in mapa_b for p, _ in c.values()]
    def choq(A, B):
        ax = np.array([a['xyz'] for a in A]); ar = np.array([VDW[a['elem']] for a in A])
        bx = np.array([b['xyz'] for b in B]); br = np.array([VDW[b['elem']] for b in B])
        d = np.sqrt(((bx[:,None,:]-ax[None,:,:])**2).sum(-1))
        return int((d < (br[:,None]+ar[None,:]-0.4)).sum()), float(d.min()), int((d < 5.0).sum())
    return choq(rpre, bpre), choq(rpost, bpost)

def rmsd_ca(mapa):
    v = [math.dist(p['xyz'], q['xyz']) for _, _, c in mapa for n, (p, q) in c.items() if n == 'CA']
    return len(v), math.sqrt(sum(x*x for x in v)/len(v)), max(v)

pre_path, post_path, lig = sys.argv[1], sys.argv[2], sys.argv[3]
pre, order_pre = load(pre_path); post, order_post = load(post_path)
lig_pre = pre[lig]
rec_pre = [a for c in pre for a in pre[c] if c != lig] if len(pre) > 2 else pre[[c for c in pre if c != lig][0]]
# ligando post: cadena cuya secuencia de resnames == la del ligando pre (match exacto)
seq_pre = ''.join(r for r, _ in residues(lig_pre))
hits = []
for c in order_post:
    s = ''.join(r for r, _ in residues(post[c]))
    if s == seq_pre: hits.append([c])
if not hits:
    # ligando partido en fragmentos: concatenacion de cadenas contiguas que reproduzca seq_pre
    acc, usado = '', []
    for c in order_post:
        s = ''.join(r for r, _ in residues(post[c]))
        if seq_pre.startswith(acc + s) or (acc + s) == (acc + s)[:len(seq_pre)] and seq_pre.startswith(acc+s):
            acc += s; usado.append(c)
            if acc == seq_pre: break
    assert acc == seq_pre, 'ligando no localizado en post por secuencia'
    hits = [usado]
lig_post_chains = hits[0]
rec_post_chains = [c for c in order_post if c not in lig_post_chains]
mapa_b, _ = mapea(lig_pre, [post[c] for c in lig_post_chains], f'LIGANDO {lig}->{"+".join(lig_post_chains)}')
mapa_r, _ = mapea(rec_pre, [post[c] for c in rec_post_chains], f'RECEPTOR->{"+".join(rec_post_chains)}')
(c0, d0, n0), (c1, d1, n1) = metricas(mapa_r, mapa_b)
nb, rb, mb = rmsd_ca(mapa_b); nr, rr, mr = rmsd_ca(mapa_r)
print(f"choques {c0}->{c1} | dmin {d0:.2f}->{d1:.2f} | contactos<5A {n0}->{n1}")
print(f"ligando CA: n={nb} RMSD={rb:.3f} max={mb:.3f} | receptor CA: n={nr} RMSD={rr:.3f} max={mr:.3f}")
