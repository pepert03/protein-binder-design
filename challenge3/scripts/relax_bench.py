#!/usr/bin/env python3
# relax_bench.py - EXPERIMENTO PEQUENO pre-declarado: la relajacion local (misma para ambos)
# absorbe los choques del protocolo rigido? Protocolo: OpenMM 8.6.1, amber14, vacio,
# Modeller addHydrogens, restriccion armonica CA k=100 kJ/mol/nm^2, LocalEnergyMinimizer
# tol=10 kJ/mol/nm max 500 iter. Control: nativo 6X18 R-P. Lectura: choques vdW>0.4 y dmin
# antes/despues con el MISMO criterio que el cribado. Comparable: mismo protocolo ambos sistemas.
import sys
import numpy as np
from openmm.app import PDBFile, Modeller, ForceField, Simulation
from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform
from openmm.unit import kilojoule_per_mole, nanometer

VDW = {'N':1.55, 'C':1.70, 'O':1.52, 'S':1.80}
def elem(n): return n.strip()[0] if n.strip()[0] in VDW else 'C'
def choques(pdb_path, ca, cb):
    seen, chains = set(), {}
    for line in open(pdb_path):
        if (line.startswith('ATOM') or line.startswith('HETATM')) and line[16] in (' ', 'A'):
            k = (line[21], line[22:26].strip(), line[12:16].strip())
            if k in seen: continue
            seen.add(k)
            e = line[76:78].strip() or elem(line[12:16])
            if e == 'H': continue
            chains.setdefault(line[21], []).append((e, [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    # robusto: cadena con mas atomos = receptor; la otra = ligando (ca/cb solo orientativos)
    ids = sorted(chains, key=lambda c: -len(chains[c]))
    A, B = chains[ids[0]], chains[ids[1]]
    ax = np.array([a[1] for a in A]); ar = np.array([VDW.get(a[0], 1.70) for a in A])
    bx = np.array([b[1] for b in B]); br = np.array([VDW.get(b[0], 1.70) for b in B])
    d = np.sqrt(((ax[:,None,:]-bx[None,:,:])**2).sum(-1))
    return int((d < (ar[:,None]+br[None,:]-0.4)).sum()), float(d.min())

def relaja(src, dst):
    from pdbfixer import PDBFixer
    fx = PDBFixer(filename=src)
    fx.findMissingResidues(); fx.findMissingAtoms(); fx.addMissingAtoms(); fx.addMissingHydrogens(7.0)
    mod = Modeller(fx.topology, fx.positions)
    ff = ForceField('amber14-all.xml')
    system = ff.createSystem(mod.topology)
    f = CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    f.addGlobalParameter('k', 100.0)
    f.addPerParticleParameter('x0'); f.addPerParticleParameter('y0'); f.addPerParticleParameter('z0')
    pos0 = mod.positions
    for i, atom in enumerate(mod.topology.atoms()):
        if atom.name == 'CA':
            p = pos0[i].value_in_unit(nanometer)
            f.addParticle(i, [p[0], p[1], p[2]])
    system.addForce(f)
    integ = LangevinMiddleIntegrator(0, 0.01, 0.0)
    sim = Simulation(mod.topology, system, integ, Platform.getPlatformByName('CPU'))
    sim.context.setPositions(pos0)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    from openmm import LocalEnergyMinimizer
    LocalEnergyMinimizer.minimize(sim.context, 10.0, 500)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    st = sim.context.getState(getPositions=True)
    with open(dst, 'w') as fh:
        from openmm.app import PDBFile as PF
        PF.writeFile(sim.topology, st.getPositions(), fh)
    return e0.value_in_unit(kilojoule_per_mole), e1.value_in_unit(kilojoule_per_mole)

for tag, src, ca, cb in (('NATIVO 6X18 R-P', 'nativo_RP.pdb', 'R', 'P'),
                         ('CANDIDATA LIBRE6_2', 'cand2_complex.pdb', 'A', 'B')):
    c0, d0 = choques(src, ca, cb)
    e0, e1 = relaja(src, src.replace('.pdb', '_rel.pdb'))
    c1, d1 = choques(src.replace('.pdb', '_rel.pdb'), ca, cb)
    print(f"{tag}: choques {c0} -> {c1} | dmin {d0:.2f} -> {d1:.2f} | E {e0:.0f} -> {e1:.0f}", flush=True)
