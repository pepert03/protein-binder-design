import numpy as np
from types import SimpleNamespace
from chain_ids import set_target_protomer_ids
for sequences,expected_entities in [([[1,2],[1,2],[1,2]],[0,0,0,0,0,0,1]),([[1,2],[1,3],[1,4]],[0,0,1,1,2,2,3])]:
 m=SimpleNamespace(_target_len=6,_binder_len=1,_pdb={'idx':{'chain':np.array(list('AABBCC'))},'batch':{'aatype':np.array(sum(sequences,[]))}},_inputs={})
 r=set_target_protomer_ids(m)
 assert r['asym_id']==[0,0,1,1,2,2,3]
 assert r['entity_id']==expected_entities
 assert r['sym_id']==([0,0,1,1,2,2,0] if sequences[0]==sequences[1] else [0]*7)
print('two CPU chain-mapping cases passed')
