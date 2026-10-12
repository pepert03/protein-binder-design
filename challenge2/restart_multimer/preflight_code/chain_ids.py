"""Preserve original target protomers in ColabDesign binder multimer inputs.
ColabDesign1.1.1 binder prep collapses all target chains into one asym_id.
Call AFTER prep_inputs and BEFORE model.run/design/predict. Not GPU-tested.
"""
import numpy as np

def set_target_protomer_ids(model, target_chains=('A','B','C')):
    t=int(model._target_len);b=int(model._binder_len)
    chains=np.asarray(model._pdb['idx']['chain'])[:t]
    if len(chains)!=t: raise ValueError('target chain mapping length mismatch')
    if set(chains)!=set(target_chains): raise ValueError('unexpected target chains')
    asym=np.empty(t+b,dtype=np.int32)
    sym=np.empty(t+b,dtype=np.int32)
    entity=np.empty(t+b,dtype=np.int32)
    aa=np.asarray(model._pdb['batch']['aatype'])[:t]
    # Crops differ per protomer and may not be treated as identical entities.
    # Full trimers with exactly identical sequences share entity, distinct sym.
    sequence_entities={};entity_counts={}
    for k,ch in enumerate(target_chains):
        positions=np.flatnonzero(chains==ch)
        seq=tuple(aa[positions].tolist())
        eid=sequence_entities.setdefault(seq,len(sequence_entities))
        sid=entity_counts.get(eid,0);entity_counts[eid]=sid+1
        asym[positions]=k;sym[positions]=sid;entity[positions]=eid
    asym[t:]=len(target_chains);sym[t:]=0;entity[t:]=len(sequence_entities)
    for key,value in [('asym_id',asym),('sym_id',sym),('entity_id',entity)]:
        model._inputs[key]=value
    return {key:np.asarray(model._inputs[key]).tolist() for key in ['asym_id','sym_id','entity_id']}
