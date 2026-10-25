#!/usr/bin/env python3
# orient_gate.py v2 - Gate de orientacion extremo-a-extremo (GLP1R ch3), ROBUSTO.
# - Soporta PDB con o sin registros MODEL (sin MODEL = un solo modelo, con aviso).
# - ERRORES explicitos (exit 2) si: no hay atomos, falta cadena A o B, alguna linea ATOM
#   de cadena A/B tiene resseq no parseable (tipico de columnas desalineadas), o faltan CA.
# - Criterios (por modelo; PASS requiere que TODOS pasen):
#   G1: >=3 res binder 1-12 con contacto <=6 A al bolsillo (U367-391)
#   G2: >=3 res binder 27-36 con contacto <=6 A al ECD (U<=139)
#   G3: 0 res binder 27-36 a <=3 A del core TM (U140-404 excl. ECL)
#   G4: 0 pares de ATOMOS intercadena <2.0 A (semantica de gate_pose.py)
# - COMPLETITUD: exige el receptor crop ch3 COMPLETO (316 resseqs del mapeo verificado
#   P43220) y binder 36 res (1..36); datos incompletos o mapeo equivocado = exit 2.
# - NOTA DE ALCANCE: el gate mide orientacion/contactos/choques all-atom del PDB dado;
#   NO certifica quimica del backbone (usar bb_check.py + revision de angulos/diedros).
# Uso: python3 orient_gate.py <complejo.pdb>   (exit 0=PASS, 1=FAIL, 2=ERROR de entrada)
import math, sys

FRAGS = [(1,23,29),(33,101,61),(108,136,136),(148,223,176),(238,300,266),(321,376,349)]
RES2UP = {r: u0 + (r - a) for a, b, u0 in FRAGS for r in range(a, b+1)}
ECL = set(range(202,228)) | set(range(291,306)) | set(range(371,384))
BOLSILLO = set(range(367,392))

def fallo(msg):
    print('ERROR ENTRADA:', msg); sys.exit(2)

def cargar(path):
    models, cur, seen = {}, None, set()
    aviso_single = False
    n_atom = 0
    try:
        fh = open(path)
    except OSError as e:
        fallo(f'no se puede abrir {path}: {e}')
    for line in fh:
        if line.startswith('MODEL'):
            cur = int(line[10:14]); models[cur] = {'A': {}, 'B': {}, '_seen': set()}; seen = models[cur]['_seen']
        elif line.startswith('ENDMDL'):
            cur = None
        elif line.startswith('ATOM'):
            n_atom += 1
            ch = line[21] if len(line) > 21 else '?'
            if ch not in ('A', 'B'):
                if len(line) > 22 and line[22] in ('A', 'B'):
                    fallo(f'columnas desalineadas (cadena en col 23): {line.rstrip()!r}')
                continue
            if cur is None:
                cur = 1; models[cur] = {'A': {}, 'B': {}, '_seen': set()}; seen = models[cur]['_seen']; aviso_single = True
            if len(line) < 54: fallo(f'linea ATOM corta ({len(line)} cols): {line.rstrip()!r}')
            if line[16] not in (' ', 'A'): continue
            try:
                rs = int(line[22:26])
            except ValueError:
                fallo(f'resseq no parseable (columnas desalineadas?): {line.rstrip()!r}')
            an = line[12:16].strip()
            if (ch, rs, an) in seen: continue
            seen.add((ch, rs, an))
            models[cur][ch].setdefault(rs, []).append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if n_atom == 0: fallo('el PDB no contiene registros ATOM')
    if not models: fallo('no se encontro ningun modelo con cadenas A/B')
    for mid, m in models.items():
        if not m['A']: fallo(f'MODEL {mid}: cadena A (receptor) vacia')
        if not m['B']: fallo(f'MODEL {mid}: cadena B (binder) vacia')
        for ch in ('A', 'B'):
            if not any(t[0] == ch and t[2] == 'CA' for t in m['_seen']):
                fallo(f'MODEL {mid}: cadena {ch} sin atomos CA')
        del m['_seen']
    if aviso_single: print('AVISO: PDB sin registros MODEL; tratado como un unico modelo.')
    esperado_a = set(RES2UP)
    esperado_b = set(range(1, 37))
    for mid, m in models.items():
        ra, rb = set(m['A']), set(m['B'])
        if ra != esperado_a:
            fallo(f'MODEL {mid}: receptor incompleto o mal mapeado ({len(ra)} res; faltan {len(esperado_a - ra)}, inesperados {len(ra - esperado_a)})')
        if rb != esperado_b:
            fallo(f'MODEL {mid}: binder incompleto o mal numerado ({len(rb)} res; esperados 1..36)')
    return models, aviso_single

def dmin(p1, p2):
    return min(math.dist(a, b) for a in p1 for b in p2)

def main():
    if len(sys.argv) < 2: fallo('falta argumento: <complejo.pdb>')
    models, aviso_single = cargar(sys.argv[1])
    todo_ok = True
    for mid, m in sorted(models.items()):
        choques = 0
        for rb in m['B']:
            for rt in m['A']:
                for p1 in m['B'][rb]:
                    for p2 in m['A'][rt]:
                        if math.dist(p1, p2) < 2.0: choques += 1
        n_bols = sum(1 for rb in m['B'] if rb <= 12 and any(dmin(m['B'][rb], m['A'][rt]) <= 6.0 for rt in m['A'] if RES2UP.get(rt, 0) in BOLSILLO))
        c_ecd = sum(1 for rb in m['B'] if rb >= 27 and any(dmin(m['B'][rb], m['A'][rt]) <= 6.0 for rt in m['A'] if 0 < RES2UP.get(rt, 999) <= 139))
        c_tm_duro = [rb for rb in m['B'] if rb >= 27 and any(dmin(m['B'][rb], m['A'][rt]) <= 3.0 for rt in m['A'] if RES2UP.get(rt, 0) >= 140 and RES2UP.get(rt, 0) not in ECL)]
        g1, g2, g3, g4 = n_bols >= 3, c_ecd >= 3, len(c_tm_duro) == 0, choques == 0
        ok = g1 and g2 and g3 and g4
        todo_ok = todo_ok and ok
        print(f'MODEL {mid}: G1_N-bolsillo={n_bols}({"OK" if g1 else "FALLO"}) G2_C-ECD={c_ecd}({"OK" if g2 else "FALLO"}) G3_C-en-TM={c_tm_duro}({"OK" if g3 else "FALLO"}) G4_choques={choques}({"OK" if g4 else "FALLO"}) -> {"PASS" if ok else "FAIL"}')
    print('GATE GLOBAL:', 'PASS' if todo_ok else 'FAIL')
    sys.exit(0 if todo_ok else 1)

main()
