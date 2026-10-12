"""Bounded AF2-multimer memory/chain-ID preflight, not a binder result."""
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
import pathlib,json,time,inspect,textwrap,numpy as np,hashlib
from colabdesign import mk_afdesign_model,clear_mem
import colabdesign.af.prep as prep_module
from folding_loss_callback import make_binder_intra_callback
from chain_ids import set_target_protomer_ids
source=textwrap.dedent(inspect.getsource(prep_module._af_prep._prep_model))
assert 'subbatch_size = None' in source and 'subbatch_size = 4' in source
patched=source.replace('subbatch_size = None','subbatch_size = 1').replace('subbatch_size = 4','subbatch_size = 1');ns=dict(prep_module.__dict__);exec(patched,ns);prep_module._af_prep._prep_model=ns['_prep_model']
pathlib.Path('prep_original.txt').write_text(source);pathlib.Path('prep_patched.txt').write_text(patched)
def dump(name,obj):pathlib.Path(name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
manifest=json.loads(pathlib.Path('crop_manifest.json').read_text());t=sum(manifest['per_chain'].values());assert t==277
# Confirm what the model actually receives, not only mutated host input arrays.
def ids_post(inputs,aux):
 aux['chain_ids_seen']={k:inputs[k] for k in ['asym_id','entity_id','sym_id']}
print('PREFLIGHT_START',jax.devices(),flush=True)
af=mk_afdesign_model(protocol='binder',use_multimer=True,data_dir='.',model_names=['model_1_multimer_v3'],best_metric='loss',use_remat=True,use_bfloat16=True,loss_callback=make_binder_intra_callback(t),post_callback=ids_post)
assert af._model_names==['model_1_multimer_v3'],af._model_names
af._cfg.model.embeddings_and_evoformer.evoformer.outer_product_mean.chunk_size=32
af.prep_inputs(pdb_filename='neighbour24.pdb',chain='A,B,C',binder_len=60,hotspot='A65,A67,A113,A115,A141,A143,A146',rm_target_seq=False)
assert af._target_len==277
mapping=set_target_protomer_ids(af);dump('chain_ids_host.json',mapping)
af.set_weights(pae=0,plddt=.1,binder_intra_pae=.4,i_pae=.1,con=1.,i_con=1.,set_defaults=True)
af.opt['con'].update({'num':2,'cutoff':14.,'seqsep':9,'binary':False});af._opt['con'].update(af.opt['con'])
af.opt['i_con'].update({'num':2,'cutoff':20.,'binary':False});af._opt['i_con'].update(af.opt['i_con'])
af.restart(seed=42101,mode='gumbel',rm_aa='C')
for recycles in [0,1]:
 st=time.time();af.run(num_models=1,num_recycles=recycles,backprop=True)
 seen={k:np.asarray(v).tolist() for k,v in af.aux['chain_ids_seen'].items()};assert seen==mapping
 logs={k:float(v) for k,v in af.aux['log'].items() if np.ndim(v)==0}
 dump('preflight_recycles'+str(recycles)+'.json',{'elapsed_s':time.time()-st,'jax':jax.__version__,'devices':str(jax.devices()),'models':af._model_names,'target_len':277,'binder_len':60,'chain_ids_seen':seen,'aggregate_log':logs,'status':'memory/chain feature preflight only, no optimized or validated binder'})
 print('PREFLIGHT_GRADIENT_PASS',recycles,time.time()-st,logs,flush=True)
print('PREFLIGHT_COMPLETE',flush=True)
