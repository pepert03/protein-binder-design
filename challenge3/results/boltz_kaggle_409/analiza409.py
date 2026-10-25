import numpy as np, json

def parse_cif(path):
    atoms=[]
    for line in open(path):
        if not line.startswith('ATOM'): continue
        f=line.split()
        atoms.append(dict(name=f[3],resname=f[5],seqid=int(f[6]),chain=f[9],
                          x=float(f[10]),y=float(f[11]),z=float(f[12]),bfac=float(f[17])))
    return atoms

def parse_pdb(path):
    atoms=[]
    for line in open(path):
        if not line.startswith('ATOM'): continue
        atoms.append(dict(name=line[12:16].strip(),resname=line[17:20].strip(),
                          seqid=int(line[22:26]),chain=line[21],
                          x=float(line[30:38]),y=float(line[38:46]),z=float(line[46:54]),
                          bfac=float(line[60:66]) if line[60:66].strip() else 0))
    return atoms

def ca_by_res(atoms,chain):
    d={}
    for a in atoms:
        if a['chain']==chain and a['name']=='CA': d[a['seqid']]=np.array([a['x'],a['y'],a['z']])
    return d

def kabsch(P,Q):
    Pc=P-P.mean(0); Qc=Q-Q.mean(0)
    H=Pc.T@Qc; U,S,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T)); D=np.diag([1,1,d])
    R=Vt.T@D@U.T
    return R, Q.mean(0)-(R@P.mean(0))

def apply(X,R,t): return (R@X.T).T+t
def rmsd(A,B): return float(np.sqrt(((A-B)**2).sum(1).mean()))

# core TM FT UniProt P43220
TM=[]
for a,b in [(140,164),(176,201),(228,251),(266,290),(306,328),(349,370),(384,404)]:
    TM+=list(range(a,b+1))
assert len(TM)==166

B='/tmp/b409/out409'
preds={}
for n in ['top02','nativo','aleatorio36']:
    preds[n]=parse_cif(f'{B}/{n}/boltz_results_glp1r409_{n}/predictions/glp1r409_{n}/glp1r409_{n}_model_0.cif')

# 440 batch (cross-batch)
B4='/tmp/boltzres/out'
preds440={}
for n in ['top02','nativo','aleatorio36']:
    preds440[n]=parse_cif(f'{B4}/{n}/boltz_results_glp1r_{n}/predictions/glp1r_{n}/glp1r_{n}_model_0.cif')

em=parse_pdb('/tmp/bench/6X18_RP.pdb')
emR=ca_by_res(em,'R'); emR={r:c for r,c in emR.items() if 24<=r<=463}
emP=ca_by_res(em,'P')
cmp=parse_pdb('/tmp/bench/complejo_v6.pdb')
cmpA={r+28:c for r,c in ca_by_res(cmp,'A').items()}
cmpB=ca_by_res(cmp,'B')

out={}
for n,atoms in preds.items():
    o={}
    for ch in ['A','B']:
        vals=[a['bfac'] for a in atoms if a['chain']==ch and a['name']=='CA']
        o[f'plddt_CA_chain{ch}']=round(np.mean(vals)/100,4)
    o['len_chainA']=len({a['seqid'] for a in atoms if a['chain']=='A' and a['name']=='CA'})
    o['len_chainB']=len({a['seqid'] for a in atoms if a['chain']=='B' and a['name']=='CA'})
    # contacts
    A=[a for a in atoms if a['chain']=='A']; Bc=[a for a in atoms if a['chain']=='B']
    bcoords=np.array([[b['x'],b['y'],b['z']] for b in Bc])
    contacts=set()
    for a in A:
        pa=np.array([a['x'],a['y'],a['z']])
        if ((bcoords-pa)**2).sum(1).min()<4.5**2:
            contacts.add(a['seqid']+23)
    o['n_contact_res_UP']=len(contacts)
    o['contact_res_UP']=sorted(contacts)
    # receptor vs 6X18
    pR={a['seqid']+23:np.array([a['x'],a['y'],a['z']]) for a in atoms if a['chain']=='A' and a['name']=='CA'}
    pB=ca_by_res(atoms,'B')
    common=sorted(set(pR)&set(emR))
    P=np.array([pR[r] for r in common]); Q=np.array([emR[r] for r in common])
    R,t=kabsch(P,Q)
    o['rec_CA_vs_6X18_global_RMSD']=round(rmsd(apply(P,R,t),Q),3); o['rec_CA_common']=len(common)
    # core TM alignment
    commonTM=sorted(set(TM)&set(pR)&set(emR))
    Ptm=np.array([pR[r] for r in commonTM]); Qtm=np.array([emR[r] for r in commonTM])
    Rtm,ttm=kabsch(Ptm,Qtm)
    o['rec_CA_coreTM166_RMSD']=round(rmsd(apply(Ptm,Rtm,ttm),Qtm),3); o['coreTM_common']=len(commonTM)
    # global RMSD computed on core-TM-aligned? no: also ECD displacement after TM alignment (descriptive)
    # binder analyses
    if n=='nativo':
        # 1:1 vs EM peptide (same sequence, 30 res) under both frames
        ref=np.array([emP[r] for r in sorted(emP)])
        pB_al_g=apply(np.array([pB[r] for r in sorted(pB)]),R,t)
        pB_al_tm=apply(np.array([pB[r] for r in sorted(pB)]),Rtm,ttm)
        o['peptido_vs_EM_RMSD_globalframe']=round(rmsd(pB_al_g,ref),3)
        o['peptido_vs_EM_RMSD_coreTMframe']=round(rmsd(pB_al_tm,ref),3)
    else:
        common2=sorted(set(pR)&set(cmpA))
        P2=np.array([pR[r] for r in common2]); Q2=np.array([cmpA[r] for r in common2])
        R2,t2=kabsch(P2,Q2)
        o['rec_CA_vs_complejo_RMSD']=round(rmsd(apply(P2,R2,t2),Q2),3)
        pB_al2=apply(np.array([pB[r] for r in sorted(pB)]),R2,t2)
        ref2=np.array([cmpB[r] for r in sorted(cmpB)])
        o['binder_CA_vs_diseno_RMSD']=round(rmsd(pB_al2,ref2),3)
        # nearest-CA to native EM peptide under coreTM frame
        pB_al_tm=apply(np.array([pB[r] for r in sorted(pB)]),Rtm,ttm)
        refp=np.array([emP[r] for r in sorted(emP)])
        d=np.sqrt(((pB_al_tm[:,None,:]-refp[None,:,:])**2).sum(2))
        md=d.min(1)
        o['binderCA_a_nativoEM_nearest_media']=round(float(md.mean()),2)
        o['binderCA_a_nativoEM_nearest_max']=round(float(md.max()),2)
    out[n]=o

# cross-batch: receptor 409 vs 440 (common UP 24-432), core-TM frame
for n in ['top02','nativo','aleatorio36']:
    a409={a['seqid']+23:np.array([a['x'],a['y'],a['z']]) for a in preds[n] if a['chain']=='A' and a['name']=='CA'}
    a440={a['seqid']+23:np.array([a['x'],a['y'],a['z']]) for a in preds440[n] if a['chain']=='A' and a['name']=='CA'}
    common=sorted(set(a409)&set(a440))
    P=np.array([a409[r] for r in common]); Q=np.array([a440[r] for r in common])
    R,t=kabsch(P,Q)
    out[n]['xbatch_rec_RMSD']=round(rmsd(apply(P,R,t),Q),3)
    # binder pose across batches, receptor aligned on core TM common
    ctm=sorted(set(TM)&set(a409)&set(a440))
    Rc,tc=kabsch(np.array([a409[r] for r in ctm]),np.array([a440[r] for r in ctm]))
    b409=ca_by_res(preds[n],'B'); b440=ca_by_res(preds440[n],'B')
    bc=sorted(set(b409)&set(b440))
    X=apply(np.array([b409[r] for r in bc]),Rc,tc); Y=np.array([b440[r] for r in bc])
    out[n]['xbatch_binder_RMSD_coreTMframe']=round(rmsd(X,Y),3)

print(json.dumps(out,indent=1))
