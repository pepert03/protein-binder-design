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
import json,time,pathlib,inspect,textwrap,hashlib,numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from colabdesign import mk_afdesign_model,clear_mem
import colabdesign.af.prep as prep_module
from folding_loss_callback import make_binder_intra_callback
from pae_scopes import pae_scopes
from audit_pdb_geometry import audit
root=pathlib.Path('output');root.mkdir(exist_ok=True)
source=textwrap.dedent(inspect.getsource(prep_module._af_prep._prep_model));assert 'subbatch_size = None' in source and 'subbatch_size = 4' in source
patched=source.replace('subbatch_size = None','subbatch_size = 1').replace('subbatch_size = 4','subbatch_size = 1');ns=dict(prep_module.__dict__)
pathlib.Path('prep_before.txt').write_text(source);pathlib.Path('prep_after.txt').write_text(patched)
def patch(on):
 exec(patched if on else source,ns);prep_module._af_prep._prep_model=ns['_prep_model']
def dump(path,obj):pathlib.Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def save(model,folder,label,seq,tL,hotspots=()):
 folder.mkdir(exist_ok=True,parents=True);aux=model.aux
 model.save_pdb(str(folder/(label+'.pdb')),get_best=False)
 logs={k:float(v) for k,v in aux['log'].items() if np.ndim(v)==0}
 allaux=aux['all'];raw=np.asarray(allaux['pae']);np.savez_compressed(folder/(label+'_pae.npz'),pae_A=raw)
 rows=[]
 for i,p in enumerate(raw):rows.append({'model':i+1,'scopes':pae_scopes(p,tL,hotspots),'binder_plddt':float(np.mean(np.asarray(allaux['plddt'])[i,tL:]))})
 # verify actual sequence in every exported model
 for m in PDBParser(QUIET=True).get_structure(label,folder/(label+'.pdb')):
  assert ''.join(seq1(r.resname) for r in list(m)[-1])==seq
 dump(folder/(label+'.json'),{'sequence':seq,'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'aggregate_log':logs,'per_model':rows})
 return logs
manifest=json.load(open('crop_manifest.json'));target=PDBParser(QUIET=True).get_structure('t','neighbour30.pdb')[0];assert sum(len(list(c)) for c in target)==366
expected=np.concatenate([np.array([r.id[1] for r in c]) for c in target])
metadata={'seed':32101,'target_length':366,'binder_length':50,'model':'model_1_ptm','optimization_recycles':0,'stages':[100,30,10],'arms':{'baseline':{'plddt':.1,'binder_intra_pae':0.},'folding':{'plddt':.5,'binder_intra_pae':.4}},'shared_weights':{'pae':.1,'i_pae':.5,'con':.5,'i_con':1.},'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pathlib.Path('.').glob('*.pdb')},'jax':jax.__version__,'devices':str(jax.devices()),'caveat':'joint folding-emphasis intervention; no causal isolation of two weights; no histidine bias, H allowed'}
dump('pilot_config.json',metadata);print('PILOT_CONFIG',json.dumps(metadata),flush=True)
for arm,weights in metadata['arms'].items():
 clear_mem();patch(True);folder=root/arm;folder.mkdir(exist_ok=True)
 af=mk_afdesign_model(protocol='binder',use_multimer=False,data_dir='.',best_metric='i_pae',model_names=['model_1_ptm'],use_remat=True,use_bfloat16=True,loss_callback=make_binder_intra_callback(366))
 af._cfg.model.embeddings_and_evoformer.evoformer.outer_product_mean.chunk_size=32
 af.prep_inputs(pdb_filename='neighbour30.pdb',chain='A,B,C',binder_len=50,hotspot='A65,A67,A113,A115,A141,A143,A146',rm_target_seq=False)
 assert af._target_len==366 and np.array_equal(af._pdb['idx']['residue'][:366],expected)
 af.opt['weights']['binder_intra_pae']=weights['binder_intra_pae'];af._opt['weights']['binder_intra_pae']=weights['binder_intra_pae']
 af.set_weights(**metadata['shared_weights'],plddt=weights['plddt'],binder_intra_pae=weights['binder_intra_pae'],set_defaults=True)
 af.restart(seed=32101,mode='gumbel',rm_aa='C');assert af.opt['weights']['binder_intra_pae']==weights['binder_intra_pae']
 t=time.time();af.run(num_models=1,num_recycles=0,backprop=True);print('FIRST_GRADIENT',arm,time.time()-t,af.aux['log'],flush=True)
 af.restart(seed=32101,mode='gumbel',rm_aa='C')
 # Stage checkpoints explicitly use current aux and sequence.
 af.design_logits(100,e_soft=1,num_models=1,num_recycles=0,ramp_recycles=False,verbose=1)
 af._tmp['seq_logits']=af.aux['seq']['logits'];s=af.get_seqs(get_best=False)[0];save(af,folder,'soft',s,366,af.opt['hotspot'])
 af.design_soft(30,e_temp=.01,num_models=1,num_recycles=0,verbose=1)
 s=af.get_seqs(get_best=False)[0];save(af,folder,'temp',s,366,af.opt['hotspot'])
 af.design_hard(10,temp=.01,dropout=False,save_best=True,num_models=1,num_recycles=0,verbose=1)
 # Best aux selected once for sequence, structure AND metrics.
 af.aux=af._tmp['best'].get('aux',af.aux);s=af.get_seqs(get_best=False)[0];lg=save(af,folder,'optimized',s,366,af.opt['hotspot']);print('ARM_DESIGN_DONE',arm,s,lg,flush=True)
 # Always inspect the fixed discrete sequence independently, even if weak.
 clear_mem();patch(False)
 for mode in ['full','alone']:
  af=mk_afdesign_model(protocol='binder' if mode=='full' else 'hallucination',use_multimer=False,data_dir='.')
  if mode=='full':af.prep_inputs(pdb_filename='full1tnf.pdb',chain='A,B,C',binder_len=50,rm_target_seq=False)
  else:af.prep_inputs(length=50)
  for seed in [32102,32103]:
   af.predict(seq=s,num_models=2,num_recycles=3,seed=seed,verbose=False);label=mode+'_'+str(seed)
   if mode=='full':
    save(af,folder,label,s,456)
    dump(folder/(label+'_geometry.json'),audit(folder/(label+'.pdb'),'B',['A'],s))
   else:
    af.save_pdb(str(folder/(label+'.pdb')),get_best=False);dump(folder/(label+'.json'),{'sequence':s,'aggregate_log':{k:float(v) for k,v in af.aux['log'].items() if np.ndim(v)==0}})
   print('VALIDATED',arm,mode,seed,af.aux['log'],flush=True)
  clear_mem()
 dump(folder/'done.json',{'arm':arm,'sequence':s,'elapsed_s':time.time()-t,'status':'completed, not binding/pH evidence'})
 print('ARM_COMPLETE',arm,flush=True)
print('PILOT_COMPLETE',flush=True)
