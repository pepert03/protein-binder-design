"""Bounded independent prediction: controls + at most two nonidentical repairs.
Run only with no other GPU design process active. Template-free binder prediction
and target-template complex prediction, not a measurement of binding or pH.
"""
import os,inspect,json,pathlib,argparse,numpy as np
import jax
if not hasattr(jax,'tree_map'):jax.tree_map=jax.tree_util.tree_map
p=argparse.ArgumentParser();p.add_argument('repair_json');p.add_argument('--target',required=True);p.add_argument('--target-chains',required=True);p.add_argument('--params',default='.');p.add_argument('--out',required=True);a=p.parse_args()
from Bio.PDB import PDBParser
from colabdesign import mk_afdesign_model,clear_mem
from colabdesign.af.loss import get_pae_loss
out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
d=json.load(open(a.repair_json)); control=d['control_sequence']
ranked=sorted(d['samples'],key=lambda x:(x['max19_hydrophobic_count'],x['score']))
seqs=[('unmodified_control',control)];seen={control}
for r in ranked:
 if r['sequence'] not in seen:seqs.append((r['name'],r['sequence']));seen.add(r['sequence'])
 if len(seqs)==3:break
records=[]
for mode in ('alone','complex'):
 clear_mem()
 model=mk_afdesign_model(protocol='hallucination' if mode=='alone' else 'binder',use_multimer=False,data_dir=a.params)
 if mode=='alone':model.prep_inputs(length=len(control))
 else:model.prep_inputs(pdb_filename=a.target,chain=a.target_chains,binder_len=len(control),rm_target_seq=False)
 for name,seq in seqs:
  model.predict(seq=seq,num_models=2,num_recycles=3,verbose=False,seed=17001)
  lg=model.aux['log'];fn=out/f'{name}_{mode}.pdb';model.save_pdb(str(fn), get_best=False)
  r=dict(name=name,mode=mode,sequence=seq,plddt=float(lg['plddt']),ptm=float(lg['ptm']) if 'ptm' in lg else None,i_pae_normalized=float(lg['i_pae']) if 'i_pae' in lg else None,i_ptm=float(lg['i_ptm']) if 'i_ptm' in lg else None,pdb=str(fn))
  records.append(r);print('PREDICT',json.dumps(r, allow_nan=False),flush=True)
 json.dump(dict(pae_source=inspect.getsource(get_pae_loss),predictions=records),open(out/'predictions.json','w'),indent=2,allow_nan=False)
# Binder CA alignment compares shape, not interface pose. Require equal residue counts.
parser=PDBParser(QUIET=True)
for name,seq in seqs:
 ca=[]
 for mode in ('alone','complex'):
  chains=list(parser.get_structure(mode,out/f'{name}_{mode}.pdb')[0]);b=chains[-1]
  xyz=np.array([r['CA'].coord for r in b if 'CA' in r])
  if len(xyz)!=len(seq):raise ValueError(('Unexpected binder length',name,mode,len(xyz)))
  ca.append(xyz)
 x,y=ca;x=x-x.mean(0);y=y-y.mean(0);u,s,v=np.linalg.svd(x.T@y);rot=u@np.diag([1,1,np.linalg.det(u@v)])@v
 rmsd=float(np.sqrt(np.mean(np.sum((x@rot-y)**2,axis=1))))
 print('SELF_CONSISTENCY',name,'binder_CA_RMSD_A',rmsd,flush=True)
 with open(out/'self_consistency.jsonl','a') as f:f.write(json.dumps(dict(name=name,binder_CA_RMSD_A=rmsd))+'\n')
