"""Prepare noninterface/surface sequence repair. Run only after inspecting source PDB.
No model-derived binding or pH claims. Absolute SASA<30 A^2 and contact<8 A
are conservative working masks, not validated acceptance thresholds.
"""
import argparse, inspect, json, pathlib, numpy as np
p=argparse.ArgumentParser(); p.add_argument('pdb'); p.add_argument('--binder',required=True); p.add_argument('--out',required=True); a=p.parse_args()
out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af.alphafold.common import residue_constants
s=PDBParser(QUIET=True).get_structure('x',a.pdb); m=s[0]
chains=list(m); b=m[a.binder]; targets=[c for c in chains if c.id!=a.binder]
if not targets:raise ValueError('No target chains')
res=[r for r in b if r.id[0]==' ' and 'CA' in r]
if len(res)!=65:raise ValueError('Expected inspected 65-residue binder')
# Reject duplicate backbone atoms or ambiguous chain assignments before sampling.
for r in res:
 for atom in ('N','CA','C','O'):
  if atom not in r:raise ValueError(('Missing backbone',r.id,atom))
tg=np.array([x.coord for c in targets for x in c.get_atoms() if x.element not in ('H','D')])
# Binder-only solvent accessibility prevents interface burial from masquerading as core.
class Binder(Select):
 def accept_chain(self,c):return c.id==a.binder
io=PDBIO();io.set_structure(s);io.save(str(out/'binder_control.pdb'),Binder())
bs=PDBParser(QUIET=True).get_structure('b',out/'binder_control.pdb')
ShrakeRupley(n_points=100).compute(bs,level='R')
sasa={r.id:float(r.sasa) for r in bs[0][a.binder]}
fixed=[]; audit=[];seq=''
for i,r in enumerate(res):
 aa=residue_constants.restype_3to1.get(r.resname,'X');seq+=aa
 xyz=np.array([x.coord for x in r if x.element not in ('H','D')])
 dist=float(np.linalg.norm(xyz[:,None]-tg[None,:],axis=-1).min())
 keep=dist<8 or sasa[r.id]<30 or aa=='H'
 if keep:fixed.append(f'{a.binder}{r.id[1]}')
 audit.append(dict(index0=i,resnum=r.id[1],aa=aa,target_distance=dist,binder_sasa_A2=sasa[r.id],fixed=keep))
if 'X' in seq:raise ValueError('Unknown residue')
fix=','.join([c.id for c in targets]+fixed)
constructor=inspect.signature(mk_mpnn_model)
kw={'seed':12001};family='original'
if 'weights' in constructor.parameters:kw['weights']='soluble';family='soluble'
model=mk_mpnn_model(**kw)
model.prep_inputs(pdb_filename=a.pdb,chain=','.join(c.id for c in targets)+','+a.binder,fix_pos=fix,rm_aa='C')
# No additional hydrophilic bias: record exact pretrained weight family.
samples=model.sample(num=1,batch=8,temperature=0.1,rescore=True)
rows=[]
for i,sample in enumerate(samples['seq']):
 seqnew=str(sample).split('/')[-1]
 assert len(seqnew)==len(seq)
 assert all(seqnew[x['index0']]==x['aa'] for x in audit if x['fixed'])
 rows.append(dict(name=f'mpnn_repair_{i:02d}',sequence=seqnew,score=float(samples['score'][i]),weight_family=family,changes=[j for j,(x,y) in enumerate(zip(seq,seqnew)) if x!=y],max19_hydrophobic_count=max(sum(z in 'AVILMFWYC' for z in seqnew[k:k+19]) for k in range(len(seqnew)-18))))
json.dump(dict(source_pdb=a.pdb,binder_chain=a.binder,target_chains=[c.id for c in targets],control_sequence=seq,weight_family=family,fix_pos=fix,audit=audit,samples=rows),open(out/'mpnn_repair.json','w'),indent=2)
print(json.dumps(rows,indent=2),flush=True)
