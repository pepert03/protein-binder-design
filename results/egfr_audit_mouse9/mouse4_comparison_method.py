import gemmi,numpy as np,json,pathlib,csv,hashlib
from scipy.spatial.distance import cdist
base=pathlib.Path('/tmp/mouse4-results');al=json.load(open('/tmp/alignment.json'));mp={int(k):v for k,v in al['mapping'].items()};h,m=al['human'],al['mouse'];orig={r['name']:r['sequence'] for r in csv.DictReader(open('/downloads/pepe_sequences-ff9b9e50.csv'))};human={r['name']:r for r in json.load(open('/tmp/egfr-report-page/src/data.json'))};out=[]
def xyz(t):return [t.pos.x,t.pos.y,t.pos.z]
def seq(c):return ''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in c)
def ca(r):return xyz(r.find_atom('CA','*'))
def fit(X,Y):
 X=np.array(X);Y=np.array(Y);xc=X.mean(0);yc=Y.mean(0);u,s,v=np.linalg.svd((X-xc).T@(Y-yc));d=np.eye(3);d[-1,-1]=np.linalg.det(u@v);R=u@d@v;t=yc-xc@R;return R,t,float(np.sqrt(np.mean(np.sum((X@R+t-Y)**2,axis=1))))
for rec in json.load(open(base/'run_records.json')):
 assert rec['completed'] and rec['returncode']==0;n=rec['name'];p=next((base/'predictions'/n).rglob('*.cif'));S=gemmi.read_structure(str(p))[0];H=gemmi.read_structure(str(next(pathlib.Path('/tmp/egfr20/predictions/'+n).rglob('*.cif'))))[0];assert seq(S['A'])==orig[n];assert seq(S['B'])==m[24:647]
 a=[(r,t) for r in S['A'] for t in r if t.element.name!='H'];b=[(r,t) for r in S['B'] for t in r if t.element.name!='H'];D=cdist([xyz(t) for r,t in a],[xyz(t) for r,t in b]);iface=sorted({b[j][0].seqid.num+24 for j in np.where(D.min(0)<5)[0]});rev={v:k for k,v in mp.items()};mapped=sorted({rev[i] for i in iface if i in rev});hi=set(human[n]['interface_uniprot']);dom=human[n]['majority_domain'];ranges={'DI':range(25,190),'DII':range(190,335),'DIII':range(335,506),'DIV':range(506,646)};nums=[i for i in ranges[dom] if i in mp and 25<=mp[i]<=647];R,t,rms=fit([ca(S['B'][mp[i]-25]) for i in nums],[ca(H['B'][i-25]) for i in nums]);X=np.array([ca(r) for r in S['A']])@R+t;Y=np.array([ca(r) for r in H['A']]);brms=float(np.sqrt(np.mean(np.sum((X-Y)**2,axis=1))));own=fit([ca(r) for r in S['A']],[ca(r) for r in H['A']])[2];i,j=np.unravel_index(D.argmin(),D.shape)
 r={'name':n,'human_ipTM':human[n]['iptm'],'mouse_ipTM':rec['confidence']['iptm'],'seconds':rec['seconds'],'exact_binder_sequence':True,'exact_mouse623':True,'human_fit_domain':dom,'target_domain_CA_RMSD_A':round(rms,3),'binder_pose_CA_RMSD_A':round(brms,3),'binder_self_fit_CA_RMSD_A':round(own,3),'mouse_interface_UniProt':iface,'mouse_interface_mapped_human':mapped,'human_interface_n':len(hi),'mouse_interface_n':len(iface),'interface_common_n':len(hi&set(mapped)),'interface_jaccard':round(len(hi&set(mapped))/len(hi|set(mapped)),3),'mouse_self_min_A':round(float(D.min()),3),'mouse_pairs_under1.5':int((D<1.5).sum()),'shortest_pair':[a[i][0].name,a[i][0].seqid.num,a[i][1].name,b[j][0].name,b[j][0].seqid.num+24,b[j][1].name],'mouse_majority_domain':max(ranges,key=lambda dom:sum(i in ranges[dom] for i in mapped))};out.append(r)
json.dump(out,open('/tmp/deep-research/egfr-audit/mouse4_comparison.json','w'),indent=2)
for r in out:print(r)
# verify all SHA256 recorded paths
bad=[]
for line in (base/'SHA256SUMS.txt').read_text().splitlines():
 sha,p=line.split('  ',1);f=base/p
 if not f.exists() or hashlib.sha256(f.read_bytes()).hexdigest()!=sha:bad.append(p)
print('HASH_CHECK',len(bad),'bad',bad[:5]);assert not bad
