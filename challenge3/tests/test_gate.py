#!/usr/bin/env python3
# test_gate.py - Regresiones de orient_gate.py (puro Python, sin numpy).
# Ejecutar desde cualquier cwd: python3 challenge3/tests/test_gate.py
import subprocess, sys, os, tempfile

BASE = os.path.dirname(os.path.abspath(__file__))
GATE = os.path.join(BASE, '..', 'scripts', 'orient_gate.py')
G3_1 = os.path.join(BASE, '..', 'results', 'traj1', 'g3_1.pdb')
SCAF2 = os.path.join(BASE, '..', 'scaffolds', 'scaffold_v2_DISCARDED.pdb')
SCAF4 = os.path.join(BASE, '..', 'scaffolds', 'scaffold_v4.pdb')
SCAF5 = os.path.join(BASE, '..', 'scaffolds', 'scaffold_v5.pdb')
TMP = tempfile.mkdtemp(prefix='gate_tests_')

def run(path):
    p = subprocess.run(['python3', GATE, path], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr

def mkpdb(lines, name):
    path = os.path.join(TMP, name)
    open(path, 'w').write('\n'.join(lines) + '\n')
    return path

def atom(serial, an, ch, rs, x, y, z, shift=0):
    base = 'ATOM  %5d %-4s%s%-3s %s%4d%s   %8.3f%8.3f%8.3f  1.00  0.00          %2s' % (serial, an, ' ', 'ALA', ch, rs, ' ', x, y, z, 'C')
    if shift: base = base[:21] + ' ' + base[21:]
    return base

res = []
# T1: PDB desalineado -> exit 2 con mensaje claro
p = mkpdb(['MODEL        1', atom(1, 'CA', 'A', 1, 0, 0, 0, shift=1), atom(2, 'CA', 'B', 1, 5, 5, 5, shift=1), 'ENDMDL', 'END'], 't1.pdb')
rc, out = run(p); res.append(('T1 desalineado -> ERROR claro', rc == 2 and ('desalinead' in out or 'resseq' in out)))
# T2: sin MODEL -> aviso + ahora ERROR de completitud (receptor sintetico de 1 res)
p = mkpdb([atom(1, 'CA', 'A', 1, 0, 0, 0), atom(2, 'CA', 'B', 1, 3, 3, 3), 'END'], 't2.pdb')
rc, out = run(p); res.append(('T2 sin MODEL + receptor incompleto -> aviso y ERROR', rc == 2 and 'AVISO' in out and 'receptor incompleto' in out))
# T3: falta cadena B -> exit 2
p = mkpdb(['MODEL        1', atom(1, 'CA', 'A', 1, 0, 0, 0), 'ENDMDL'], 't3.pdb')
rc, out = run(p); res.append(('T3 sin cadena B -> ERROR claro', rc == 2 and 'cadena B (binder) vacia' in out))
# T4: archivo sin ATOM -> exit 2
p = mkpdb(['HEADER    VACIO'], 't4.pdb')
rc, out = run(p); res.append(('T4 sin ATOM -> ERROR claro', rc == 2 and 'no contiene registros ATOM' in out))
# T5: g3_1 real (receptor completo) -> FAIL (exit 1)
rc, out = run(G3_1); res.append(('T5 g3_1 -> FAIL exit 1', rc == 1 and 'GATE GLOBAL: FAIL' in out))
# T6: scaffold_v2 (geometria de pose valida para el gate, quimica descartada) -> PASS (exit 0)
rc, out = run(SCAF2); res.append(('T6 scaffold_v2 -> PASS exit 0', rc == 0 and 'GATE GLOBAL: PASS' in out))
# T7: receptor incompleto (scaffold_v2 menos un fragmento de cadena A) -> exit 2
lines = [l.rstrip('\n') for l in open(SCAF2) if not (l.startswith('ATOM') and l[21] == 'A' and 33 <= int(l[22:26]) <= 101)]
p = mkpdb(lines, 't7.pdb')
rc, out = run(p); res.append(('T7 receptor incompleto -> ERROR claro', rc == 2 and 'receptor incompleto' in out))
# T8: binder incompleto (scaffold_v2 menos B36) -> exit 2
lines = [l.rstrip('\n') for l in open(SCAF2) if not (l.startswith('ATOM') and l[21] == 'B' and int(l[22:26]) == 36)]
p = mkpdb(lines, 't8.pdb')
rc, out = run(p); res.append(('T8 binder incompleto -> ERROR claro', rc == 2 and 'binder incompleto' in out))

# T9: scaffold_v4 real -> PASS (exit 0)
rc, out = run(SCAF4); res.append(('T9 scaffold_v4 -> PASS exit 0', rc == 0 and 'GATE GLOBAL: PASS' in out))

# T10: scaffold_v5 (input vigente) -> PASS (exit 0)
rc, out = run(SCAF5); res.append(('T10 scaffold_v5 -> PASS exit 0', rc == 0 and 'GATE GLOBAL: PASS' in out))

for nombre, ok in res: print(('OK ' if ok else 'FALLO '), nombre)
print('TODOS PASS' if all(ok for _, ok in res) else 'HAY FALLOS')
sys.exit(0 if all(ok for _, ok in res) else 1)
