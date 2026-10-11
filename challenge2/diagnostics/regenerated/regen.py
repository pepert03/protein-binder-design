import jax
for _name in ['tree_map','tree_flatten','tree_unflatten','tree_leaves','tree_structure']:
 if not hasattr(jax,_name):setattr(jax,_name,getattr(jax.tree_util,_name))

import jax
if not hasattr(jax,'tree_map'):jax.tree_map=jax.tree_util.tree_map
import jax.numpy as jnp
_orig=jnp.clip
def _clip(a,*args,a_min=None,a_max=None,**kw):
 if a_min is not None:kw['min']=a_min
 if a_max is not None:kw['max']=a_max
 return _orig(a,*args,**kw)
jnp.clip=_clip
import numpy as np,json,inspect,pathlib
from colabdesign import mk_afdesign_model,clear_mem
from colabdesign.af.loss import get_pae_loss
from colabdesign.mpnn import mk_mpnn_model
print('GPU',jax.devices(),flush=True)
print('PAE_SOURCE',inspect.getsource(get_pae_loss),flush=True)
print('MPNN_API',inspect.signature(mk_mpnn_model),flush=True)
SEQ='CKDVIIIKHPYPVTPDPRWHQTVQRYRERAEREELLGDGVRNGCLCCLLLLLLLILLLKELTKVP'
L=[l for l in open('1tnf.pdb') if l.startswith('ATOM') and l[21] in 'AB']
hot={('A',r) for r in [65,67,113,115,141,143,146]}
ca=[(l[21],int(l[22:26]),np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])) for l in L if l[12:16]==' CA ']
hp=[p for c,r,p in ca if (c,r) in hot]
keep={(c,r) for c,r,p in ca if min(np.linalg.norm(p-q) for q in hp)<20}
k2=set(keep)
for ch in 'AB':
 rs=sorted(r for c,r in keep if c==ch)
 for a,b in zip(rs,rs[1:]):
  if 1<b-a<=6:k2|={(ch,x) for x in range(a,b)}
 rs=sorted(r for c,r in k2 if c==ch);frag=[];cur=[]
 for r in rs:
  if cur and r-cur[-1]>1:frag.append(cur);cur=[]
  cur.append(r)
 if cur:frag.append(cur)
 for f in frag:
  if len(f)<4:k2-={(ch,x) for x in f}
open('tnf_crop.pdb','w').write(''.join(l for l in L if (l[21],int(l[22:26])) in k2))
clear_mem();af=mk_afdesign_model(protocol='binder',use_multimer=False,data_dir='.')
af.prep_inputs(pdb_filename='tnf_crop.pdb',chain='A,B',binder_len=len(SEQ),rm_target_seq=False)
af.predict(seq=SEQ,num_models=2,num_recycles=3,seed=17001,verbose=False)
af.save_pdb('regenerated_candidate.pdb')
lg=af.aux['log'];json.dump({'sequence':SEQ,'log':{k:float(v) for k,v in lg.items() if np.ndim(v)==0},'pae_source':inspect.getsource(get_pae_loss),'mpnn_constructor':str(inspect.signature(mk_mpnn_model))},open('regenerated_candidate.json','w'),indent=2)
print('REGENERATED',lg,flush=True)
# Parse exact chain identities and save before any repair.
from Bio.PDB import PDBParser
s=PDBParser(QUIET=True).get_structure('x','regenerated_candidate.pdb')
print('CHAINS',[(c.id,len([r for r in c if 'CA' in r])) for c in s[0]],flush=True)
