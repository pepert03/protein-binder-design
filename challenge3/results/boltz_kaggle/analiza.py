import numpy as np, json, re

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
    # alinea P sobre Q; P,Q Nx3
    Pc=P-P.mean(0); Qc=Q-Q.mean(0)
    H=Pc.T@Qc; U,S,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T)); D=np.diag([1,1,d])
    R=Vt.T@D@U.T
    return R, Q.mean(0)-(R@P.mean(0))

def apply(X,R,t): return (R@X.T).T+t

def rmsd(A,B): return float(np.sqrt(((A-B)**2).sum(1).mean()))

out={}
preds={}
for n in ['top02','nativo','aleatorio36']:
    p=f'out/{n}/boltz_results_glp1r_{n}/predictions/glp1r_{n}/glp1r_{n}_model_0.cif'
    preds[n]=parse_cif(p)

# 1) pLDDT por cadena (B-factor CA, escala 0-100)
for n,atoms in preds.items():
    for ch in ['A','B']:
        vals=[a['bfac'] for a in atoms if a['chain']==ch and a['name']=='CA']
        out.setdefault(n,{})[f'plddt_CA_chain{ch}']=round(np.mean(vals)/100,4)

# 2) contactos de interfaz B-A < 4.5 A (atomos pesados), receptor en numeracion UP (+23)
for n,atoms in preds.items():
    A=[a for a in atoms if a['chain']=='A']; B=[a for a in atoms if a['chain']=='B']
    bcoords=np.array([[b['x'],b['y'],b['z']] for b in B])
    contacts=set()
    for a in A:
        pa=np.array([a['x'],a['y'],a['z']])
        if ((bcoords-pa)**2).sum(1).min()<4.5**2:
            contacts.add(a['seqid']+23)
    out[n]['n_contact_res_UP']=len(contacts)
    out[n]['contact_res_UP']=sorted(contacts)

# 3) receptor predicho vs EM 6X18_RP (UP 29-423), y poses
em=parse_pdb('/tmp/bench/6X18_RP.pdb')
emR=ca_by_res(em,'R'); emR={r:c for r,c in emR.items() if 24<=r<=463}
emP=ca_by_res(em,'P')  # GLP-1 7-36
cmp=parse_pdb('/tmp/bench/complejo_v6.pdb')
cmpA={r+28:c for r,c in ca_by_res(cmp,'A').items()}  # file->UP
cmpB=ca_by_res(cmp,'B')

for n in ['top02','nativo']:
    atoms=preds[n]
    pR={a['seqid']+23:np.array([a['x'],a['y'],a['z']]) for a in atoms if a['chain']=='A' and a['name']=='CA'}
    pB=ca_by_res(atoms,'B')
    common=sorted(set(pR)&set(emR))
    P=np.array([pR[r] for r in common]); Q=np.array([emR[r] for r in common])
    R,t=kabsch(P,Q)
    out[n]['rec_CA_vs_6X18_RMSD']=round(rmsd(apply(P,R,t),Q),3)
    out[n]['rec_CA_common']=len(common)
    pB_al=apply(np.array(list(pB.values())),R,t)
    if n=='nativo':
        ref=np.array([emP[r] for r in sorted(emP)])
        out[n]['binder_CA_vs_nativoEM_RMSD']=round(rmsd(pB_al,ref),3)
    else:
        # binder disenado en complejo_v6: superponer pred-receptor sobre complejo-receptor
        common2=sorted(set(pR)&set(cmpA))
        P2=np.array([pR[r] for r in common2]); Q2=np.array([cmpA[r] for r in common2])
        R2,t2=kabsch(P2,Q2)
        out[n]['rec_CA_vs_complejo_RMSD']=round(rmsd(apply(P2,R2,t2),Q2),3)
        pB_al2=apply(np.array(list(pB.values())),R2,t2)
        ref2=np.array([cmpB[r] for r in sorted(cmpB)])
        out[n]['binder_CA_vs_diseno_RMSD']=round(rmsd(pB_al2,ref2),3)

# 4) binder predicho top02 vs pose EM del nativo (misma referencia receptor)
atoms=preds['top02']
pR={a['seqid']+23:np.array([a['x'],a['y'],a['z']]) for a in atoms if a['chain']=='A' and a['name']=='CA'}
pB=ca_by_res(atoms,'B')
common=sorted(set(pR)&set(emR))
R,t=kabsch(np.array([pR[r] for r in common]),np.array([emR[r] for r in common]))
# binder top02 tiene 36 res, nativo EM 30: comparar solo si misma longitud -> no; medir solapamiento espacial por vecinos mas cercanos
pB_al=apply(np.array(list(pB.values())),R,t)
ref=np.array([emP[r] for r in sorted(emP)])
d=np.sqrt(((pB_al[:,None,:]-ref[None,:,:])**2).sum(2))
out['top02']['binderCA_a_nativoEM_min_dist_media']=round(float(d.min(1).mean()),2)

print(json.dumps(out,indent=1))
