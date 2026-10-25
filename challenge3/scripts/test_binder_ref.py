#!/usr/bin/env python3
# test_binder_ref.py - REGRESION: referencia binder para evaluaciones (ESMFold/RMSD/packing).
# Falla EXPLICITAMENTE si: no hay cadena 'B', no tiene 36 CA, o el backbone no es identico
# entre scaffold_v6 y cand2_complex. PROHIBIDO fallback por ventana/RMSD sobre el receptor
# (ese bug selecciono un segmento del receptor, resseq 148-183, como si fuera el binder).
# Uso: python3 test_binder_ref.py [raiz_repo]   (exit 0 = OK, exit 1 = FALLO)
import sys, math

BB = {'N', 'CA', 'C', 'O'}

def parse(path):
    out = []
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')):
            out.append((line[21], int(line[22:26]), line[12:16].strip(),
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return out

def fail(msg):
    print(f'FALLO REGRESION binder-ref: {msg}')
    sys.exit(1)

def ca_por_cadena(atoms):
    ch = {}
    for c, r, n, x in atoms:
        if n == 'CA':
            ch.setdefault(c, []).append((r, x))
    return ch

def main(root):
    scaf = parse(f'{root}/challenge3/scaffolds/scaffold_v6.pdb')
    cand = parse(f'{root}/challenge3/results/relax1/cand2_complex.pdb')

    for tag, atoms in [('scaffold_v6', scaf), ('cand2_complex', cand)]:
        ch = ca_por_cadena(atoms)
        if 'B' not in ch:
            fail(f'{tag}: no existe cadena B (cadenas: {sorted(ch)})')
        if len(ch['B']) != 36:
            fail(f'{tag}: cadena B tiene {len(ch["B"])} CA, esperados 36')
        hits = [c for c, v in ch.items() if len(v) == 36]
        if hits != ['B']:
            fail(f'{tag}: cadenas con 36 CA = {hits}, esperado solo [B]; NO elegir por RMSD/ventana')

    sb = [(r, n, x) for c, r, n, x in scaf if c == 'B' and n in BB]
    cb = [(r, n, x) for c, r, n, x in cand if c == 'B' and n in BB]
    if len(sb) != len(cb):
        fail(f'backbone binder: {len(sb)} vs {len(cb)} atomos')
    dmax = max(math.dist(a[2], b[2]) for a, b in zip(sb, cb) if (a[0], a[1]) == (b[0], b[1]))
    if dmax > 1e-3:
        fail(f'backbone binder scaffold vs complejo difiere {dmax:.4f} A (>1e-3)')

    sa = [(r, n, x) for c, r, n, x in scaf if c == 'A' and n in BB]
    ca_ = {(r, n): x for c, r, n, x in cand if c == 'A' and n in BB}
    miss = [(r, n) for r, n, x in sa if (r, n) not in ca_]
    if miss:
        fail(f'backbone receptor: {len(miss)} atomos del scaffold ausentes en el complejo')
    dmax = max(math.dist(x, ca_[(r, n)]) for r, n, x in sa)
    if dmax > 1e-3:
        fail(f'backbone receptor scaffold vs complejo difiere {dmax:.4f} A (>1e-3)')

    print('OK regresion binder-ref: cadena B unica con 36 CA; backbone binder y receptor identicos scaffold<->complejo (sin fallback por RMSD)')

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
