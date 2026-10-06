import gemmi,numpy as np,pathlib,json,csv
from scipy.optimize import minimize
refs=gemmi.read_structure('/tmp/1IVO.pdb')[0]['A'];rr=[r for r in refs if r.het_flag=='A'];Y=np.array([[r.find_atom('CA','*').pos.x,r.find_atom('CA','*').pos.y,r.find_atom('CA','*').pos.z] for r in rr]);out=[]
for row in csv.DictReader(open('/downloads/pepe_sequences-ff9b9e50.csv')):
 name=row['name'];c=gemmi.read_structure(str(next(pathlib.Path('/tmp/egfr20/predictions/'+name).rglob('*.cif'))))[0]['A'];X=np.array([[r.find_atom('CA','*').pos.x,r.find_atom('CA','*').pos.y,r.find_atom('CA','*').pos.z] for r in c]);L=len(X);d0=max(.5,1.24*(max(L,16)-15)**(1/3)-1.8);best=None
 for i in range(len(Y)-L+1):
  Z=Y[i:i+L];xc=X.mean(0);zc=Z.mean(0);u,s,v=np.linalg.svd((X-xc).T@(Z-zc));d=np.eye(3);d[-1,-1]=np.linalg.det(u@v);R=u@d@v;dist=np.linalg.norm((X-xc)@R+zc-Z,axis=1);score=float(np.mean(1/(1+(dist/d0)**2)));rms=float(np.sqrt(np.mean(dist**2)))
  if best is None or score>best['qTM_like']:best={'name':name,'qTM_like':round(score,4),'rmsd_A':round(rms,3),'n_CA':L,'target_start_UniProt':rr[i].seqid.num+24,'target_end_UniProt':rr[i+L-1].seqid.num+24}
 out.append(best)
json.dump(out,open('/tmp/deep-research/egfr-audit/EGFR_narrow_structure_screen.json','w'),indent=2)
for r in sorted(out,key=lambda x:-x['qTM_like']):print(r)
