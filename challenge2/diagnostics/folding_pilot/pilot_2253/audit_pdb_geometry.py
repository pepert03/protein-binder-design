"""Offline geometry audit. No affinity/pH claims. Explicit chain roles required."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from scipy.spatial.distance import cdist

def audit(path, binder_chain, target_chains, expected_sequence=None):
    rows=[]
    for model in PDBParser(QUIET=False).get_structure('audit',path):
        if binder_chain in target_chains: raise ValueError('Binder also listed as target')
        if any(c not in model for c in [binder_chain]+target_chains): raise ValueError('Missing explicit chain')
        binder=model[binder_chain]; targets=[model[c] for c in target_chains]
        residues=[r for r in binder if r.id[0]==' ']
        sequence=''.join(seq1(r.resname, custom_map={'UNK':'X'}) for r in residues)
        if expected_sequence is not None and sequence != expected_sequence: raise ValueError('Exact binder sequence mismatch')
        ba=[a for r in residues for a in r if a.element not in ('H','D')]
        ta=[a for c in targets for r in c if r.id[0]==' ' for a in r if a.element not in ('H','D')]
        if not ba or not ta: raise ValueError('No binder/target heavy atoms')
        d=cdist([a.coord for a in ba],[a.coord for a in ta]); contacts={}
        for cutoff in (4,5,8):
            indices=np.where(d<cutoff)
            contacts[str(cutoff)]={'heavy_atom_pairs':len(indices[0]),'binder_residues':len({ba[i].parent.id for i in indices[0]}),'target_residues':len({(ta[j].parent.parent.id,ta[j].parent.id) for j in indices[1]})}
        gaps=[]
        for c in [binder]+targets:
            rs=[r for r in c if r.id[0]==' ' and 'CA' in r]
            for a,b in zip(rs,rs[1:]):
                dist=float(np.linalg.norm(a['CA'].coord-b['CA'].coord))
                if b.id[1]-a.id[1]!=1 or dist>4.5:
                    gaps.append({'chain':c.id,'from_residue':a.id[1],'to_residue':b.id[1],'CA_distance_A':dist})
        ca=np.array([r['CA'].coord for r in residues if 'CA' in r]); rg=float(np.sqrt(np.mean(np.sum((ca-ca.mean(0))**2,axis=1)))) if len(ca) else None
        rows.append({'model':model.id+1,'sequence':sequence,'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),'binder_length':len(sequence),'target_residues':sum(1 for c in targets for r in c if r.id[0]==' '),'min_heavy_distance_A':float(d.min()),'pairs_under2A':int((d<2).sum()),'contacts_A':contacts,'binder_CA_radius_gyration_A':rg,'chain_discontinuities':gaps})
    return {'file_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'binder_chain':binder_chain,'target_chains':target_chains,'method':'Every MODEL; heavy atoms exclude H/D; strict distances; gaps flagged, not declared covalent errors; Rg is descriptive, not acceptance test','models':rows}

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('pdb');p.add_argument('--binder',required=True);p.add_argument('--targets',required=True);p.add_argument('--sequence');p.add_argument('--out',required=True);a=p.parse_args()
    result=audit(a.pdb,a.binder,a.targets.split(','),a.sequence)
    Path(a.out).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
