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
from chain_ids import set_target_protomer_ids
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

root=pathlib.Path('cohort_rescue0/42111');done=root/'done.json'
if not done.exists():raise SystemExit('NO_DONE_NO_VALIDATION')
d=json.loads(done.read_text())
if d['status']!='requires_independent_alone_and_intact_validation':raise SystemExit('REJECTED_NO_VALIDATION')
s=d['sequence'];clear_mem();patch(False)
for mode in ['alone','full']:
 af=mk_afdesign_model(protocol='hallucination' if mode=='alone' else 'binder',use_multimer=False,data_dir='.')
 if mode=='full':
  af.prep_inputs(pdb_filename='full1tnf.pdb',chain='A,B,C',binder_len=60,rm_target_seq=False)
  set_target_protomer_ids(af)
 else:af.prep_inputs(length=60)
 af.predict(seq=s,num_models=2,num_recycles=3,seed=42121,verbose=False)
 label=mode+'_independent'
 if mode=='full':save(af,root,label,s,456)
 else:
  af.save_pdb(str(root/(label+'.pdb')),get_best=False)
  dump(root/(label+'.json'),{'sequence':s,'aggregate_log':{k:float(v) for k,v in af.aux['log'].items() if np.ndim(v)==0}})
 print('VALIDATED',mode,s,af.aux['log'],flush=True);clear_mem()
print('INDEPENDENT_VALIDATION_COMPLETE',flush=True)
