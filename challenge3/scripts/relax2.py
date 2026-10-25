#!/usr/bin/env python3
# relax2.py - Receptor FIJADO con restriccion k=10000 kJ/mol/nm2 a TODOS sus atomos (RMSD se verifica ~0);
# ligando libre. amber14 vacio, L-BFGS tol=10 maxIter=500. Mismo protocolo en los 3 sistemas.
# Sistemas: C1 nativo (P libre), C2 poly-ALA empaquetado (B libre), C3 LIBRE6_2 (B libre).
import math
import numpy as np
from openmm.app import PDBFile
from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, LocalEnergyMinimizer
from openmm.unit import nanometer, kilojoule_per_mole
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
    assert len(hits) == 1, f"cadenas con {n_ca} CA: {hits}"
    return hits[0]

def choques(path, n_ca):
    ch = cadenas(path)
    lc = ligand_chain(ch, n_ca)
    rec = [a for c, v in ch.items() if c != lc for a in v if a[0] in VDW and a[1] != 'OXT']
    lig = [a for a in ch[lc] if a[0] in VDW and a[1] != 'OXT']
    rx = np.array([a[2] for a in rec]); rr = np.array([VDW[a[0]] for a in rec])
    bx = np.array([a[2] for a in lig]); br = np.array([VDW[a[0]] for a in lig])
    d = np.sqrt(((bx[:,None,:]-rx[None,:,:])**2).sum(-1))
    return int((d < (br[:,None]+rr[None,:]-0.4)).sum()), float(d.min()), int((d < 5.0).sum())

def rmsd_ca(pre, post, n_ca):
    chp, chq = cadenas(pre), cadenas(post)
    lp = ligand_chain(chp, n_ca); lq = ligand_chain(chq, n_ca)
    ca_p = [a[2] for a in chp[lp] if a[1] == 'CA']; ca_q = [a[2] for a in chq[lq] if a[1] == 'CA']
    v = [math.dist(ca_p[i], ca_q[i]) for i in range(n_ca)]
    rl = math.sqrt(sum(x*x for x in v)/len(v)); ml = max(v)
    ca_rp = [a[2] for c, v in chp.items() if c != lp for a in v if a[1] == 'CA']
    ca_rq = [a[2] for c, v in chq.items() if c != lq for a in v if a[1] == 'CA']
    n = min(len(ca_rp), len(ca_rq))
    v = [math.dist(ca_rp[i], ca_rq[i]) for i in range(n)]
    return rl, ml, math.sqrt(sum(x*x for x in v)/n), max(v)

def relaja_fijo(src, dst, n_ca):
    fx = PDBFixer(filename=src)
    fx.findMissingResidues(); fx.findMissingAtoms(); fx.addMissingAtoms(); fx.addMissingHydrogens(7.0)
    from openmm.app import ForceField, Simulation
    ff = ForceField('amber14-all.xml')
    system = ff.createSystem(fx.topology)
    f = CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    f.addGlobalParameter('k', 10000.0)
    f.addPerParticleParameter('x0'); f.addPerParticleParameter('y0'); f.addPerParticleParameter('z0')
    # identificar ligando por n CA por cadena en la topologia
    chain_ca = {}
    for ch in fx.topology.chains():
        chain_ca[ch.id] = sum(1 for r in ch.residues() for a in r.atoms() if a.name == 'CA')
    lig_id = next(cid for cid, n in chain_ca.items() if n == n_ca)
    nfix = 0
    for i, atom in enumerate(fx.topology.atoms()):
        if atom.residue.chain.id != lig_id:
            p = fx.positions[i].value_in_unit(nanometer)
            f.addParticle(i, [p[0], p[1], p[2]]); nfix += 1
    system.addForce(f)
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
    dst = src.replace('.pdb', '_rel2.pdb')
    c0, d0, n0 = choques(src, nca)
    e0, e1, nfix = relaja_fijo(src, dst, nca)
    c1, d1, n1 = choques(dst, nca)
    rl, ml, rr, mr = rmsd_ca(src, dst, nca)
    print(f"{tag}: choques {c0}->{c1} | dmin {d0:.2f}->{d1:.2f} | contactos {n0}->{n1} | lig CA RMSD {rl:.3f} (max {ml:.3f}) | receptor CA RMSD {rr:.4f} (max {mr:.4f}) | E {e0:.0f}->{e1:.0f} | atomos fijados {nfix}", flush=True)
