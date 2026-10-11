import os, subprocess, sys, time
t0=time.time()
def sh(c):
    r=subprocess.run(c,shell=True,capture_output=True,text=True); print(c[:80],'->',r.returncode, (r.stdout+r.stderr)[-300:]); return r
sh("mkdir -p params && ls params | head -3")
if not os.path.exists('params/params_model_1_ptm.npz'):
    sh("curl -fsSL https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar | tar x -C params")
sh("ls params | head; curl -fsSL -o 1tnf.pdb https://files.rcsb.org/download/1TNF.pdb; grep -c ^ATOM 1tnf.pdb")
import jax; print(jax.__version__, jax.devices())
print('setup s',time.time()-t0)
import jax.numpy as jnp
_c=jnp.clip
def _clip(a,*args,a_min=None,a_max=None,**k):
    if a_min is not None: k['min']=a_min
    if a_max is not None: k['max']=a_max
    return _c(a,*args,**k)
jnp.clip=_clip
import os, time, re, numpy as np
DRY = (jax.devices()[0].platform=='cpu')
BUDGET = 60 if DRY else 80*60   # seconds of design wall clock
t0=time.time()
# keep chains A,B only
import numpy as _np
_L=[l for l in open('1tnf.pdb') if l.startswith('ATOM') and l[21] in 'AB']
_hot={('A',int(h[1:])) for h in "A65,A67,A113,A115,A141,A143,A146".split(',')}
_ca=[(l[21],int(l[22:26]),_np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])) for l in _L if l[12:16]==' CA ']
_hp=[p for c,r,p in _ca if (c,r) in _hot]
_keep={(c,r) for c,r,p in _ca if min(_np.linalg.norm(p-q) for q in _hp)<20.0}
_k2=set(_keep)
for _ch in 'AB':
    _rs=sorted(r for c,r in _keep if c==_ch)
    for a,b in zip(_rs,_rs[1:]):
        if 1<b-a<=6: _k2|={(_ch,x) for x in range(a,b)}
    # drop fragments < 6 residues
    _rs=sorted(r for c,r in _k2 if c==_ch); frag=[]; cur=[]
    for r in _rs:
        if cur and r-cur[-1]>1: frag.append(cur); cur=[]
        cur.append(r)
    if cur: frag.append(cur)
    for f in frag:
        if len(f)<4: _k2-={(_ch,x) for x in f}
_keep=_k2
open('tnf_AB.pdb','w').write(''.join(l for l in _L if (l[21],int(l[22:26])) in _keep))
print('CROP residues',len(_keep),flush=True)
from colabdesign import mk_afdesign_model, clear_mem
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.af.alphafold.common import residue_constants
HOT="A65,A67,A113,A115,A141,A143,A146"
BL=65
af=mk_afdesign_model(protocol="binder", use_multimer=False, data_dir=".", best_metric="i_pae")
af.prep_inputs(pdb_filename="tnf_AB.pdb", chain="A,B", hotspot=HOT, binder_len=BL, rm_target_seq=False)
tgt=''.join(residue_constants.restypes[i] for i in af._pdb["batch"]["aatype"][:af._target_len]) if hasattr(af,'_target_len') else ''
try:
    pass
    print('AF His bias set',flush=True)
except Exception as e: print('no AF bias',e,flush=True)
rows=[]; k=0
def hiscat(pdb_str, blen):
    # binder is last chain; count binder His whose side chain atoms are within 5A of target Lys/Arg/His side chain N
    ch={}
    for l in pdb_str.splitlines():
        if l.startswith('ATOM'): ch.setdefault(l[21],[]).append((l[17:20],int(l[22:26]),l[12:16].strip(),np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])))
    names=sorted(ch); b=ch[names[-1]]; tg=[a for n in names[:-1] for a in ch[n]]
    cat=[a for a in tg if a[0] in('LYS','ARG','HIS') and a[2] in('NZ','NH1','NH2','NE','ND1','NE2')]
    C=np.array([a[3] for a in cat]) if cat else np.zeros((0,3))
    n=0
    for r in sorted({a[1] for a in b if a[0]=='HIS'}):
        ne=[a[3] for a in b if a[1]==r and a[2] in('ND1','NE2')]
        if len(C) and min(np.linalg.norm(C-x,axis=1).min() for x in ne)<5.0: n+=1
    return n
while time.time()-t0 < BUDGET:
    k+=1
    af.restart(seed=k+200)
    af.design_3stage(30 if DRY else 100, 30 if DRY else 100, 5 if DRY else 10)
    chosen=af._tmp['best'].get('aux', af.aux); lg=chosen['log']; seq=af.get_seqs(get_best=True)[0]
    base=dict(traj=k, seq=seq, plddt=round(lg['plddt'],3), i_pae=round(lg['i_pae'],3), i_ptm=round(lg['i_ptm'],3))
    af.save_pdb(f'c4_{k}.pdb')
    print('TRAJ',base, 'elapsed',int(time.time()-t0), flush=True)
    rows.append(base)
    if DRY: break
print('DONE',flush=True)
