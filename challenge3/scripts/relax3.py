#!/usr/bin/env python3
# relax3.py - PROTOCOLO A REAL: receptor fijado por MASA 0 (sin restriccion armonica).
# Verificacion empirica exigida: receptor CA RMSD debe ser EXACTAMENTE 0.
# amber14 vacio, L-BFGS tol=10 maxIter=500. Sistemas: C1 nativo(P libre), C2 poly-ALA(B), C3 LIBRE6_2(B).
import math
import numpy as np
from openmm.app import PDBFile, ForceField, Simulation
from openmm import LangevinMiddleIntegrator, Platform, LocalEnergyMinimizer
from openmm.unit import nanometer, kilojoule_per_mole, dalton
from pdbfixer import PDBFixer

VDW = {'N':1.55, 'C':1.70, 'O':1.52, 'S':1.80}
def elem(n): return n.strip()[0] if n.strip()[0] in VDW else 'C'
def cadenas(path):
    ch = {}
    for line in open(path):
        if line.startswith('ATOM') and line[16] in (' ', 'A'):
            e = line[76:78].strip() or elem(line[12:16])
            ch.setdefault(line[21], []).append((e, line[12:16].strip(), [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return ch
def ligand_chain(ch, n_ca):
    hits = [c for c, v in ch.items() if sum(1 for a in v if a[1] == 'CA') == n_ca]
    assert len(hits) == 1, hits
    return hits[0]
def choques(path, n_ca):
    ch = cadenas(path); lc = ligand_chain(ch, n_ca)
    rec = [a for c, v in ch.items() if c != lc for a in v if a[0] in VDW and a[1] != 'OXT']
    lig = [a for a in ch[lc] if a[0] in VDW and a[1] != 'OXT']
    rx = np.array([a[2] for a in rec]); rr = np.array([VDW[a[0]] for a in rec])
    bx = np.array([a[2] for a in lig]); br = np.array([VDW[a[0]] for a in lig])
    d = np.sqrt(((bx[:,None,:]-rx[None,:,:])**2).sum(-1))
    return int((d < (br[:,None]+rr[None,:]-0.4)).sum()), float(d.min()), int((d < 5.0).sum())
def rmsd_ca(pre, post, n_ca):
    chp, chq = cadenas(pre), cadenas(post)
    lp, lq = ligand_chain(chp, n_ca), ligand_chain(chq, n_ca)
    ca_p = [a[2] for a in chp[lp] if a[1] == 'CA']; ca_q = [a[2] for a in chq[lq] if a[1] == 'CA']
    v = [math.dist(ca_p[i], ca_q[i]) for i in range(n_ca)]
    rl, ml = math.sqrt(sum(x*x for x in v)/len(v)), max(v)
    ca_rp = [a[2] for c, v in chp.items() if c != lp for a in v if a[1] == 'CA']
    ca_rq = [a[2] for c, v in chq.items() if c != lq for a in v if a[1] == 'CA']
    v = [math.dist(ca_rp[i], ca_rq[i]) for i in range(min(len(ca_rp), len(ca_rq)))]
    return rl, ml, math.sqrt(sum(x*x for x in v)/len(v)), max(v)
def relaja_masa0(src, dst, n_ca):
    fx = PDBFixer(filename=src)
    fx.findMissingResidues(); fx.findMissingAtoms(); fx.addMissingAtoms(); fx.addMissingHydrogens(7.0)
    ff = ForceField('amber14-all.xml')
    system = ff.createSystem(fx.topology)
    chain_ca = {}
    for chn in fx.topology.chains():
        chain_ca[chn.id] = sum(1 for r in chn.residues() for a in r.atoms() if a.name == 'CA')
    lig_id = next(cid for cid, n in chain_ca.items() if n == n_ca)
    nfix = 0
    for i, atom in enumerate(fx.topology.atoms()):
        if atom.residue.chain.id != lig_id:
            system.setParticleMass(i, 0.0); nfix += 1
    sim = Simulation(fx.topology, system, LangevinMiddleIntegrator(0, 0.01, 0.0), Platform.getPlatformByName('CPU'))
    sim.context.setPositions(fx.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    LocalEnergyMinimizer.minimize(sim.context, 10.0, 500)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    with open(dst, 'w') as fh:
        PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(), fh)
    return e0.value_in_unit(kilojoule_per_mole), e1.value_in_unit(kilojoule_per_mole), nfix
for tag, src, nca in (('C1 NATIVO', 'nativo_RP.pdb', 30),
                      ('C2 poly-ALA', 'ala6_complex.pdb', 36),
                      ('C3 LIBRE6_2', 'cand2_complex.pdb', 36)):
    dst = src.replace('.pdb', '_relA.pdb')
    c0, d0, n0 = choques(src, nca)
    e0, e1, nfix = relaja_masa0(src, dst, nca)
    c1, d1, n1 = choques(dst, nca)
    rl, ml, rr, mr = rmsd_ca(src, dst, nca)
    print(f"{tag}: choques {c0}->{c1} | dmin {d0:.2f}->{d1:.2f} | contactos {n0}->{n1} | lig CA RMSD {rl:.3f} (max {ml:.3f}) | RECEPTOR CA RMSD {rr:.6f} (max {mr:.6f}) | E {e0:.0f}->{e1:.0f} | masa0 {nfix}", flush=True)
