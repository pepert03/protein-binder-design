"""Model-isolated His-cation proximity. A heuristic, not H-bond/pH prediction."""
from Bio.PDB import PDBParser
from scipy.spatial.distance import cdist
import numpy as np
from io import StringIO
RING={'ND1','NE2'}
CATION={'LYS':{'NZ'},'ARG':{'NH1','NH2','NE'},'HIS':RING}

def model_counts(pdb_text,binder_chain=None,target_chains=None,cutoff_A=5.0):
    if cutoff_A<=0:raise ValueError('Positive cutoff required')
    records=[]
    for model in PDBParser(QUIET=True).get_structure('his',StringIO(pdb_text)):
        chains=list(model)
        if len(chains)<2:raise ValueError('Target and binder needed')
        # Legacy AfDesign saves target as A and binder as B; override for other layouts.
        bid=binder_chain if binder_chain is not None else chains[-1].id
        tids=target_chains if target_chains is not None else [c.id for c in chains if c.id!=bid]
        if bid in tids or not tids or bid not in model or any(c not in model for c in tids):raise ValueError('Invalid explicit chain roles')
        target=[a for cid in tids for r in model[cid] if r.resname in CATION for a in r if a.name in CATION[r.resname]]
        xyz=np.array([a.coord for a in target]); hits=[]
        for r in model[bid]:
            if r.resname!='HIS':continue
            ring=[a.coord for a in r if a.name in RING]
            if len(ring)!=2:raise ValueError('Incomplete histidine ring')
            if len(target):
                d=cdist(np.array(ring),xyz);i,j=np.unravel_index(d.argmin(),d.shape)
                if d[i,j]<cutoff_A:
                    a=target[j];hits.append({'binder_residue':r.id[1],'min_distance_A':float(d[i,j]),'target_chain':a.parent.parent.id,'target_residue':a.parent.id[1],'target_atom':a.name})
        records.append({'model':model.id+1,'binder_chain':bid,'target_chains':tids,'his_count_within_cutoff':len(hits),'cutoff_A':cutoff_A,'contacts':hits})
    if not records:raise ValueError('No models')
    return records

def conservative_count(pdb_text,blen):
    """Legacy integer summary: minimum across models, never pooled atoms."""
    models=model_counts(pdb_text)
    structure=PDBParser(QUIET=True).get_structure('length',StringIO(pdb_text))
    for model,row in zip(structure,models):
        if len(list(model[row['binder_chain']]))!=blen:raise ValueError('Binder length mismatch')
    return min(row['his_count_within_cutoff'] for row in models)
