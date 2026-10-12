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
print('COHORT_START',jax.devices(),flush=True)
af=mk_afdesign_model(protocol='binder',use_multimer=True,data_dir='.',model_names=['model_1_multimer_v3'],best_metric='plddt',use_remat=True,use_bfloat16=True,loss_callback=make_binder_intra_callback(t),post_callback=ids_post)
assert af._model_names==['model_1_multimer_v3'],af._model_names
af._cfg.model.embeddings_and_evoformer.evoformer.outer_product_mean.chunk_size=32
af.prep_inputs(pdb_filename='neighbour24.pdb',chain='A,B,C',binder_len=60,hotspot='A65,A67,A113,A115,A141,A143,A146',rm_target_seq=False)
assert af._target_len==277
mapping=set_target_protomer_ids(af);dump('chain_ids_host.json',mapping)
af.opt['weights']['binder_intra_pae']=0.;af._opt['weights']['binder_intra_pae']=0.
af.set_weights(pae=0,plddt=.1,binder_intra_pae=0.,i_pae=.1,con=1.,i_con=1.,set_defaults=True)
af.opt['con'].update({'num':2,'cutoff':14.,'seqsep':9,'binary':False});af._opt['con'].update(af.opt['con'])
af.opt['i_con'].update({'num':2,'cutoff':20.,'binary':False});af._opt['i_con'].update(af.opt['i_con'])
af.restart(seed=42101,mode='gumbel',rm_aa='C')

from pae_scopes import pae_scopes
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
root=pathlib.Path('cohort_rescue0');root.mkdir(exist_ok=True)
started=time.time();dump('cohort_rescue0_config.json',{'seed_list':[42111],'actual_weights':af.opt['weights'],'fix':'explicit zero intra weight rescue after .4 cohort rejection; exact PDB sequence export','soft_sequence_caveat':'PDB pseudo argmax sequence is not independently validated','selection':'best binder pLDDT rather than cross-temperature loss; verify interface separately'})
def save(folder,label):
 folder.mkdir(exist_ok=True,parents=True)
 af.save_pdb(str(folder/(label+'.pdb')),get_best=False)
 seq=''.join('ARNDCQEGHILKMFPSTWYV'[int(a)] for a in np.asarray(af.aux['all']['aatype'])[0,-60:])
 for model in PDBParser(QUIET=True).get_structure(label,folder/(label+'.pdb')):
  assert ''.join(seq1(r.resname) for r in list(model.get_residues())[-60:])==seq
 raw=np.asarray(af.aux['all']['pae']);np.savez_compressed(folder/(label+'_pae.npz'),pae_A=raw)
 logs={k:float(v) for k,v in af.aux['log'].items() if np.ndim(v)==0}
 dump(folder/(label+'.json'),{'sequence':seq,'aggregate_log':logs,'per_model':[pae_scopes(x,277,af.opt['hotspot']) for x in raw],'caveat':'aggregate multimer iPTM includes target-target interfaces; not binder-specific; soft state not discrete validation'})
 return seq,logs
for seed in [42111]:
 if time.time()-started>9000:break
 folder=root/str(seed)
 af.restart(seed=seed,mode='gumbel',rm_aa='C')
 af.design_logits(50,e_soft=.9,num_models=1,num_recycles=1,ramp_recycles=False,save_best=True,verbose=1)
 af.aux=af._tmp['best'].get('aux',af.aux)
 seq,logs=save(folder,'screen50')
 if logs['plddt']<=.65:
  dump(folder/'done.json',{'seed':seed,'status':'early_rejected_low_fold_confidence','plddt':logs['plddt'],'no_binding_evidence':True});print('SEED_REJECTED',seed,logs,flush=True);continue
 af.clear_best();af.design_logits(25,e_soft=1,num_models=1,num_recycles=1,ramp_recycles=False,save_best=True,verbose=1)
 af._tmp['seq_logits']=af.aux['seq']['logits']
 af.clear_best();af.design_soft(45,e_temp=.01,num_models=1,num_recycles=1,ramp_recycles=False,save_best=True,verbose=1)
 af.aux=af._tmp['best'].get('aux',af.aux);seq,logs=save(folder,'temp')
 if logs['plddt']<=.65:
  dump(folder/'done.json',{'seed':seed,'status':'temp_rejected_low_fold_confidence','no_binding_evidence':True});continue
 af.clear_best();af.design_hard(5,temp=.01,num_models=1,num_recycles=1,dropout=False,save_best=True,verbose=1)
 af.aux=af._tmp['best'].get('aux',af.aux);seq,logs=save(folder,'discrete')
 dump(folder/'done.json',{'seed':seed,'sequence':seq,'status':'requires_independent_alone_and_intact_validation' if logs['plddt']>.65 else 'hard_rejected_low_fold_confidence','no_binding_evidence':True})
 print('SEED_DISCRETE',seed,seq,logs,flush=True)
print('COHORT_COMPLETE',flush=True)
