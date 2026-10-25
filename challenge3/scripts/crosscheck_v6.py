# crosscheck_v4.py - Verificacion INDEPENDIENTE con biopython de scaffold_v6.pdb:
# (a) phi/psi/omega/angulos/enlaces medidos con Bio.PDB.vectors (codigo distinto al mio)
# (b) RMSD de mi helice ideal canonica vs helice REAL de 6X18 (U233-244) sobre N/CA/C/O
import math
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.vectors import calc_dihedral, calc_angle
from Bio.PDB.Superimposer import Superimposer as SVDSuperimposer

p = PDBParser(QUIET=True)
st = p.get_structure('s', (__import__("sys").argv[1] if len(__import__("sys").argv)>1 else "scaffold_v6.pdb"))
B = st[0]['B']
resid = [r for r in B if r.id[0] == ' ']
ph, ps, om, angs = [], [], [], []
for i in range(1, len(resid)-1):
    r0, r1, r2 = resid[i-1], resid[i], resid[i+1]
    ph.append(math.degrees(calc_dihedral(r0['C'].get_vector(), r1['N'].get_vector(), r1['CA'].get_vector(), r1['C'].get_vector())))
    ps.append(math.degrees(calc_dihedral(r1['N'].get_vector(), r1['CA'].get_vector(), r1['C'].get_vector(), r2['N'].get_vector())))
    om.append(math.degrees(calc_dihedral(r1['CA'].get_vector(), r1['C'].get_vector(), r2['N'].get_vector(), r2['CA'].get_vector())))
    angs.append(math.degrees(calc_angle(r1['N'].get_vector(), r1['CA'].get_vector(), r1['C'].get_vector())))
print(f'biopython: phi {min(ph):.1f}..{max(ph):.1f} | psi {min(ps):.1f}..{max(ps):.1f} | omega {min(om):.1f}..{max(om):.1f} | N-CA-C {min(angs):.1f}..{max(angs):.1f}')
def bond(r, a, b): return float(np.linalg.norm(r[a].coord - r[b].coord))
print('N-CA', f"{min(bond(r,'N','CA') for r in resid):.3f}..{max(bond(r,'N','CA') for r in resid):.3f}",
      '| CA-C', f"{min(bond(r,'CA','C') for r in resid):.3f}..{max(bond(r,'CA','C') for r in resid):.3f}",
      '| C-O', f"{min(bond(r,'C','O') for r in resid):.3f}..{max(bond(r,'C','O') for r in resid):.3f}")
pep = [float(np.linalg.norm(resid[i]['C'].coord - resid[i+1]['N'].coord)) for i in range(len(resid)-1)]
print('C-N pept', f'{min(pep):.3f}..{max(pep):.3f}')

# (b) ideal (res 1-12 del scaffold) vs real U233-244
ideal_atoms, real_atoms = [], []
ref = p.get_structure('r', __import__("os").environ.get("RP_PDB", "/tmp/ch3pack/challenge3/target/6X18_RP.pdb"))[0]['R']
rmap = {r.id[1]: r for r in ref}
for j, ures in enumerate(range(233, 245)):
    rb = resid[j]
    for an in ('N', 'CA', 'C', 'O'):
        ideal_atoms.append(rb[an].coord); real_atoms.append(rmap[ures][an].coord)
def kabsch_rmsd(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    R = V @ D @ Wt
    return float(np.sqrt((((Pc @ R) - Qc) ** 2).sum() / len(P)))
print(f'RMSD ideal(res1-12) vs real U233-244 (N/CA/C/O, 48 atomos): {kabsch_rmsd(np.array(ideal_atoms), np.array(real_atoms)):.3f} A')
