import gemmi,numpy as np,json,csv,pathlib
from scipy.spatial.distance import cdist
from Bio.Align import PairwiseAligner
from Bio import SeqIO
OUT=pathlib.Path('/tmp/deep-research/egfr-audit');OUT.mkdir(exist_ok=True)
al=json.load(open('/tmp/alignment.json'));hm={int(k):v for k,v in al['mapping'].items()};h,m=al['human'],al['mouse']
rows=list(csv.DictReader(open('/downloads/pepe_sequences-ff9b9e50.csv')))
def aa(r):return gemmi.find_tabulated_residue(r.name).one_letter_code
def chain_atoms(c,poly=True):
 return [(r,a) for r in c for a in r if a.element.name!='H' and (not poly or r.het_flag=='A')]
def fit(x,y):
 x=np.asarray(x);y=np.asarray(y);xc=x.mean(0);yc=y.mean(0);u,s,v=np.linalg.svd((x-xc).T@(y-yc));d=np.eye(3);d[-1,-1]=np.linalg.det(u@v);R=u@d@v;t=yc-xc@R;return R,t,float(np.sqrt(np.mean(np.sum((x@R+t-y)**2,axis=1))))
def xyz(a):return [a.pos.x,a.pos.y,a.pos.z]
refs=[]
for pdb in ['1IVO','1NQL']:
 s=gemmi.read_structure('/tmp/'+pdb+'.pdb'); c=s[0]['A'];rs={r.seqid.num+24:r for r in c if r.het_flag=='A' and r.find_atom('CA','*')}; gly=[(r,a) for ch in s[0] if ch.name in (['A','E'] if pdb=='1IVO' else ['A','C','D','E']) for r in ch if r.name in ['NAG','MAN','BMA','FUC','GAL'] for a in r if a.element.name!='H'];refs.append((pdb,rs,gly))
res=[]
for row in rows:
 name=row['name'];p=next(pathlib.Path('/tmp/egfr20/predictions/'+name).rglob('*.cif'));s=gemmi.read_structure(str(p));a,b=s[0]['A'],s[0]['B'];assert ''.join(aa(r) for r in a)==row['sequence'];assert ''.join(aa(r) for r in b)==h[24:645]
 at,bt=chain_atoms(a),chain_atoms(b);X=np.array([xyz(t) for r,t in at]);Y=np.array([xyz(t) for r,t in bt]);D=cdist(X,Y);ix=np.where(D.min(0)<5)[0];iface=sorted({bt[i][0].seqid.num+24 for i in ix});con=[i for i in iface if i in hm and h[i-1]==m[hm[i]-1]];dif=[f'{h[i-1]}{i}>{m[hm[i]-1]}{hm[i]}' for i in iface if i in hm and h[i-1]!=m[hm[i]-1]]
 cf=json.load(open(next(p.parent.glob('confidence*.json'))));r={'name':name,'iptm':cf['iptm'],'interface_uniprot':iface,'conserved':len(con),'interface_n':len(iface),'mutations':dif,'min_interface_A':round(float(D.min()),3),'contacts_under_2.5':int((D<2.5).sum()),'contacts_under_1.5':int((D<1.5).sum()),'glycan_context':[]}
 # flags for sequon sites within 10A even when no sugars resolved
 sites=[i+25 for i in range(619) if h[24+i]=='N' and h[25+i]!='P' and h[26+i] in 'ST'];near=[]
 for z in sites:
  rr=b[z-25];d=float(cdist(X,np.array([xyz(t) for t in rr if t.element.name!='H'])).min());
  if d<10:near.append({'uniprot_N':z,'distance_A':round(d,2)})
 r['near_sequons']=near
 for pdb,rs,gly in refs:
  for mode in ['DI','DII','DIII','DIV','local']:
   dom={'DI':range(25,190),'DII':range(190,335),'DIII':range(335,506),'DIV':range(506,646)}; majority=max(dom,key=lambda d:sum(i in dom[d] for i in iface));r['majority_domain']=majority; nums=[i for i in (dom[mode] if mode!='local' else [i for i in iface if i in dom[majority]]) if i in rs and aa(rs[i])==h[i-1]]
   if len(nums)<8:continue
   x=[xyz(b[i-25].find_atom('CA','*')) for i in nums];y=[xyz(rs[i].find_atom('CA','*')) for i in nums];R,t,rms=fit(x,y);XX=X@R+t;other=[(i,at) for i,rr in rs.items() if i not in dom[mode if mode!='local' else majority] for at in rr if at.element.name!='H'];O=np.array([xyz(at) for i,at in other]);dg=float(cdist(XX,np.array([xyz(at) for rr,at in gly])).min()) if gly else None
   rc={'pdb':pdb,'fit':mode,'n_CA':len(nums),'rmsd_A':round(rms,2),'resolved_glycan_min_A':round(dg,2) if dg else None,'other_domain_min_A':round(float(cdist(XX,O).min()),2),'other_domain_pairs_under_2.5':int((cdist(XX,O)<2.5).sum())};r['glycan_context'].append(rc)
 res.append(r)
json.dump(res,open(OUT/'audit.json','w'),indent=2)
with open(OUT/'audit.csv','w') as f:
 w=csv.writer(f);w.writerow(['name','ipTM','interface_conserved','interface_n','human_mouse_changes','min_interface_A','pairs_under_2.5','near_N_sequons','1IVO_DIII_RMSD','1IVO_glycan_min','1IVO_other_domain_min','1NQL_glycan_min']);
 for r in res:
  g={z['pdb']:z for z in r['glycan_context'] if z['fit']=='DIII'};w.writerow([r['name'],r['iptm'],r['conserved'],r['interface_n'],';'.join(r['mutations']),r['min_interface_A'],r['contacts_under_2.5'],str(r['near_sequons']),g.get('1IVO',{}).get('rmsd_A'),g.get('1IVO',{}).get('resolved_glycan_min_A'),g.get('1IVO',{}).get('other_domain_min_A'),g.get('1NQL',{}).get('resolved_glycan_min_A')])
for r in sorted(res,key=lambda r:-r['iptm']): print(r['name'],round(r['iptm'],3),f"{r['conserved']}/{r['interface_n']}",r['mutations'],r['contacts_under_2.5'],r['near_sequons'],[(x['pdb'],x['fit'],x['rmsd_A'],x['resolved_glycan_min_A'],x['other_domain_min_A']) for x in r['glycan_context']])
