import pathlib,json,gemmi,csv
out=pathlib.Path('/downloads/Proteinbase_pH_analysis');r=json.load(open(out/'analysis.json'))
for x in r:
 s=gemmi.read_structure(str(next((pathlib.Path('/downloads/Proteinbase_Boltz_results/predictions')/x['name']).rglob('*.cif')))); cs={c.name:c for c in s[0]};ring=[]
 for ch,other in [('A','B'),('B','A')]:
  for res in cs[ch]:
   if res.name!='HIS':continue
   for rr in cs[other]:
    pairs=[(a.pos.dist(b.pos),a,b) for a in res if a.name in ['ND1','CD2','CE1','NE2','CG'] for b in rr if b.element.name!='H']
    if not pairs:continue
    d,a,b=min(pairs,key=lambda p:p[0])
    if d<=5:ring.append(dict(his_chain=ch,his_res=res.seqid.num,partner_res=rr.seqid.num,partner_name=rr.name,distance=round(d,3),atoms=[a.name,b.name],his_plddt=round(a.b_iso,2)))
 x['histidine_ring_contacts']=ring
(out/'analysis.json').write_text(json.dumps(r,indent=2))
with open(out/'fold_confidence.csv','w') as f:
 w=csv.writer(f);w.writerow(['candidate','ipTM','binder_His_count','EGFR_His_ring_contacts','binder_His_ring_contacts']);
 for x in r:w.writerow([x['name'],x['iptm'],len(x['his_positions']['A']),sum(y['his_chain']=='B' for y in x['histidine_ring_contacts']),sum(y['his_chain']=='A' for y in x['histidine_ring_contacts'])])
for x in r[:5]:print(x['name'],x['histidine_ring_contacts'])
