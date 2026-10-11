import jax
for n in ['tree_map','tree_flatten','tree_unflatten','tree_leaves','tree_structure']:
 if not hasattr(jax,n):setattr(jax,n,getattr(jax.tree_util,n))
import jax.numpy as jnp
_orig=jnp.clip
def _clip(a,*args,a_min=None,a_max=None,**kw):
 if a_min is not None:kw['min']=a_min
 if a_max is not None:kw['max']=a_max
 return _orig(a,*args,**kw)
jnp.clip=_clip

import time,json,pathlib,numpy as np,inspect
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser
from colabdesign import mk_afdesign_model,clear_mem
from colabdesign.af.loss import get_pae_loss
import colabdesign.af.prep as prep_module,textwrap
source=textwrap.dedent(inspect.getsource(prep_module._af_prep._prep_model))
assert 'subbatch_size = None' in source and 'subbatch_size = 4' in source
patched=source.replace('subbatch_size = None','subbatch_size = 1').replace('subbatch_size = 4','subbatch_size = 1')
namespace=dict(prep_module.__dict__);exec(patched,namespace);prep_module._af_prep._prep_model=namespace['_prep_model']
pathlib.Path('prep_patch_before.txt').write_text(source);pathlib.Path('prep_patch_after.txt').write_text(patched)
print('GRADIENT_SUBBATCH_PATCH',patched,flush=True)
out=pathlib.Path('output');out.mkdir(exist_ok=True)
clear_mem();af=mk_afdesign_model(protocol='binder',use_multimer=False,data_dir='.',best_metric='i_pae',model_names=['model_1_ptm'],use_remat=True,use_bfloat16=True)
af._cfg.model.embeddings_and_evoformer.evoformer.outer_product_mean.chunk_size=32
af.prep_inputs(pdb_filename='input/full1tnf.pdb',chain='A,B,C',binder_len=50,hotspot='A65,A67,A113,A115,A141,A143,A146',rm_target_seq=False)
af.set_weights(plddt=.1,pae=.1,i_pae=.5,con=.5,i_con=1.,set_defaults=True)
af.restart(seed=30101,mode='gumbel',rm_aa='C')
print('FRESH_CONFIG',json.dumps(dict(seed=30101,binder_len=50,target_len=af._target_len,weights=af.opt['weights'],hotspot=list(map(int,af.opt['hotspot'])))),flush=True)
print('MEMORY_CONFIG',af._cfg.model.global_config,flush=True)
t0=time.time();af.run(num_models=1,num_recycles=0,backprop=True)
print('FULL_FIRST_GRADIENT_OK',time.time()-t0,flush=True)
af.restart(seed=30101,mode='gumbel',rm_aa='C')
af.design_3stage(100,30,10,num_models=1,num_recycles=0,ramp_recycles=False,verbose=1)
seq=af.get_seqs()[0];af.save_pdb(str(out/'trajectory01_optimized.pdb'))
lg={k:float(v) for k,v in af.aux['log'].items() if np.ndim(v)==0}
json.dump(dict(sequence=seq,log=lg,elapsed_s=time.time()-t0),open(out/'trajectory01_optimized.json','w'),indent=2)
print('FRESH_TRAJECTORY',seq,lg,flush=True)
# Discrete independent prediction; preserves optimized coordinates separately.
clear_mem();exec(source,namespace);prep_module._af_prep._prep_model=namespace['_prep_model']
af=mk_afdesign_model(protocol='binder',use_multimer=False,data_dir='.')
af.prep_inputs(pdb_filename='input/full1tnf.pdb',chain='A,B,C',binder_len=50,rm_target_seq=False)
af.predict(seq=seq,num_models=2,num_recycles=3,verbose=False,seed=30102)
af.save_pdb(str(out/'trajectory01_repredicted.pdb'))
lg={k:float(v) for k,v in af.aux['log'].items() if np.ndim(v)==0}
audit=[]
for f in out.glob('*.pdb'):
 for m in PDBParser(QUIET=True).get_structure(f.stem,f):
  ch=list(m);ba=[a for a in ch[-1].get_atoms() if a.element!='H'];ta=[a for c in ch[:-1] for a in c.get_atoms() if a.element!='H'];d=cdist([a.coord for a in ba],[a.coord for a in ta]);audit.append(dict(pdb=f.name,model=m.id,target_residues=sum(len(list(c)) for c in ch[:-1]),binder_residues=len(list(ch[-1])),min_heavy_distance_A=float(d.min()),pairs_under2A=int((d<2).sum()),binder_residues_under4A=len({ba[i].parent.id[1] for i,j in zip(*np.where(d<4))})))
json.dump(dict(sequence=seq,repredict_log=lg,pae_source=inspect.getsource(get_pae_loss),clash_audit=audit),open(out/'trajectory01_audit.json','w'),indent=2)
print('FRESH_REPREDICT',lg,'CLASH',audit,flush=True)
