import numpy as np, json, hashlib, re

BASE='/tmp/zs2/campana_zeroshot_409'

def parse_cif(path):
    atoms=[]
    for line in open(path):
        if not line.startswith('ATOM'): continue
        f=line.split()
        atoms.append(dict(name=f[3],resname=f[5],seqid=int(f[6]),chain=f[9],
                          x=float(f[10]),y=float(f[11]),z=float(f[12]),bfac=float(f[17]),
                          elem=f[2] if len(f)>2 else ''))
    return atoms

def kabsch(P,Q):
    Pc=P-P.mean(0); Qc=Q-Q.mean(0)
    H=Pc.T@Qc; U,S,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T)); D=np.diag([1,1,d])
    R=Vt.T@D@U.T
    return R, Q.mean(0)-(R@P.mean(0))
def apply(X,R,t): return (R@X.T).T+t
def rmsd(A,B): return float(np.sqrt(((A-B)**2).sum(1).mean()))

# --- 1) YAML identity
official='689647a8c8b799bdf43c9b753e1e83ff5446b1f7'
print('== 1) YAML identity ==')
yseq={}
for n in ['gA_1','gA_2','gB_1','gB_2','nativo','aleatorio36']:
    t=open(f'{BASE}/filter_inputs/f409_{n}.yaml').read()
    seqs=re.findall(r'sequence:\s*([A-Z]+)', t)
    h=hashlib.sha1(seqs[0].encode()).hexdigest()
    yseq[n]=seqs[1]
    print(f'{n}: rec {len(seqs[0])}aa sha1==official: {h==official}; binder {len(seqs[1])}aa')
# fasta match
fa={}
name=None
for line in open(f'{BASE}/designs/sequences.fasta'):
    line=line.strip()
    if line.startswith('>'): name=line[1:]
    elif line: fa[name]=line
for n,s in fa.items():
    print(f'fasta {n}: == yaml: {yseq.get(n)==s}')
# nativo/aleatorio vs canonical repo YAMLs
can={}
for n,path in [('nativo','/tmp/b409/inputs409/glp1r409_nativo.yaml'),('aleatorio36','/tmp/b409/inputs409/glp1r409_aleatorio36.yaml')]:
    t=open(path).read(); can[n]=re.findall(r'sequence:\s*([A-Z]+)', t)[1]
print('nativo binder == canonical 409 batch:', yseq['nativo']==can['nativo'])
print('aleatorio binder == canonical 409 batch:', yseq['aleatorio36']==can['aleatorio36'])

# --- 2) epitope set
epi=set()
rows=open(f'{BASE}/docs/epitope_contacts.tsv').read().splitlines()[1:]
for line in rows:
    c=line.split('\t')
    if c and c[0].isdigit(): epi.add(int(c[2]))
print('== 2) epitope: rows', len(rows), 'unique uniprot', len(epi), 'max', max(epi))

# --- 3) contacts recomputed (5.0 A, heavy atoms, rn+23)
print('== 3) contacts (5.0 A, heavy) ==')
resj=json.load(open(f'{BASE}/filter_outputs/contacts409.json'))
mycont={}
for key in resj:
    sub,name=key.split('/',1) if key.startswith('out') else (None,None)
    if key.startswith('out_seed2/'): nm=key.split('/')[1]; cifp=f'{BASE}/filter_outputs/out_seed2/{nm}/boltz_results_f409_{nm}/predictions/f409_{nm}/f409_{nm}_model_0.cif'
    else: nm=key.split('/')[1]; cifp=f'{BASE}/filter_outputs/out/{nm}/boltz_results_f409_{nm}/predictions/f409_{nm}/f409_{nm}_model_0.cif'
    atoms=parse_cif(cifp)
    A=[a for a in atoms if a['chain']=='A']; B=[a for a in atoms if a['chain']=='B']
    bc=np.array([[b['x'],b['y'],b['z']] for b in B])
    hit=set()
    for a in A:
        pa=np.array([a['x'],a['y'],a['z']])
        if ((bc-pa)**2).sum(1).min()<5.0**2: hit.add(a['seqid']+23)
    ov=sorted(hit & epi)
    mine={'n':len(hit),'epi':len(ov)}
    theirs={'n':resj[key]['n_contacts'],'epi':resj[key]['n_epitope_overlap']}
    match = (mine==theirs) and sorted(hit)==resj[key]['contacts']
    mycont[key]=hit
    print(f'{key}: mine {mine} theirs {theirs} EXACT {match}')

# --- 4) pLDDT cross-check from cif b-factors
print('== 4) pLDDT vs results.json ==')
res1=json.load(open(f'{BASE}/filter_outputs/results.json'))
res2=json.load(open(f'{BASE}/filter_outputs/results_seed2.json'))
def plddt_all(cifp):
    atoms=parse_cif(cifp)
    return np.mean([a['bfac'] for a in atoms])/100
for nm in ['gA_1','gA_2','gB_1','gB_2','nativo','aleatorio36']:
    cifp=f'{BASE}/filter_outputs/out/{nm}/boltz_results_f409_{nm}/predictions/f409_{nm}/f409_{nm}_model_0.cif'
    p=plddt_all(cifp)
    print(f'{nm}: cif all-atom pLDDT {p:.4f} vs results.json complex_plddt {res1[nm]["complex_plddt"]:.4f} diff {abs(p-res1[nm]["complex_plddt"]):.4f}')

# --- 5) cross-seed pose stability (core TM 166 CA frame)
print('== 5) cross-seed stability ==')
TM=[]
for a,b in [(140,164),(176,201),(228,251),(266,290),(306,328),(349,370),(384,404)]:
    TM+=list(range(a,b+1))
for nm in ['gB_1','gB_2','gA_2','nativo','aleatorio36']:
    c1=parse_cif(f'{BASE}/filter_outputs/out/{nm}/boltz_results_f409_{nm}/predictions/f409_{nm}/f409_{nm}_model_0.cif')
    c2=parse_cif(f'{BASE}/filter_outputs/out_seed2/{nm}/boltz_results_f409_{nm}/predictions/f409_{nm}/f409_{nm}_model_0.cif')
    def ca(atoms,chain,off=0):
        return {a['seqid']+off:np.array([a['x'],a['y'],a['z']]) for a in atoms if a['chain']==chain and a['name']=='CA'}
    r1=ca(c1,'A',23); r2=ca(c2,'A',23)
    ctm=sorted(set(TM)&set(r1)&set(r2))
    R,t=kabsch(np.array([r1[x] for x in ctm]),np.array([r2[x] for x in ctm]))
    b1=ca(c1,'B'); b2=ca(c2,'B')
    bc=sorted(set(b1)&set(b2))
    X=apply(np.array([b1[x] for x in bc]),R,t); Y=np.array([b2[x] for x in bc])
    o1=len(mycont[f'out/{nm}']&epi); o2=len(mycont[f'out_seed2/{nm}']&epi)
    inter=len(mycont[f'out/{nm}']&mycont[f'out_seed2/{nm}']); uni=len(mycont[f'out/{nm}']|mycont[f'out_seed2/{nm}'])
    print(f'{nm}: binder CA RMSD seed1-vs-seed2 {rmsd(X,Y):.2f} A (frame core TM); epitope overlap {o1}/{o2}; contact-set Jaccard {inter}/{uni}={inter/uni:.2f}')
