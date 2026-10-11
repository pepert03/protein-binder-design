import sys,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from pae_scopes import pae_scopes
p=np.array([[0,0,2,2],[0,0,18,18],[4,20,6,6],[4,20,6,6.]])
x=pae_scopes(p,2,[0])
assert x['binder_to_hotspots']['mean_A']==3
assert x['binder_to_entire_target']['mean_A']==11
assert x['binder_intra']['mean_A']==6
assert x['binder_to_entire_target']['pair_count']==4
for a,t,h in [(p,0,[]),(p,2,[2]),(p,2,[0,0]),(np.full((3,3),np.nan),2,[])]:
 try:pae_scopes(a,t,h)
 except ValueError:pass
 else:raise AssertionError('Negative case failed')
print('PASS scopes, asymmetry symmetrization, units, shape and index guards')
