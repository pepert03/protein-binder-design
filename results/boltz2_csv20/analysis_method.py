import pathlib,json,math,contextlib,io,gemmi,propka.run,csv
root=pathlib.Path('/downloads/Proteinbase_Boltz_results');out=pathlib.Path('/downloads/Proteinbase_pH_analysis');out.mkdir(exist_ok=True)
records=json.load(open(root/'run_records.json'));results=[]
for rec in sorted(records,key=lambda x:x['confidence']['iptm'],reverse=True):
 name=rec['name'];s=gemmi.read_structure(str(next((root/'predictions'/name).rglob('*.cif'))));model=s[0]
 atoms={c.name:[(r,a) for r in c for a in r if a.element.name!='H'] for c in model};contacts={};interface={'A':set(),'B':set()}
 for r,a in atoms['A']:
  for rr,aa in atoms['B']:
   d=a.pos.dist(aa.pos)
   if d<=5:
    interface['A'].add(r.seqid.num);interface['B'].add(rr.seqid.num)
    if r.name=='HIS' or rr.name=='HIS':
     key=(r.seqid.num,r.name,rr.seqid.num,rr.name)
     if key not in contacts or d<contacts[key]['distance']:
      contacts[key]={'distance':round(d,3),'atoms':[a.name,aa.name],'plddt':[round(a.b_iso,2),round(aa.b_iso,2)]}
 hs={c.name:[r.seqid.num for r in c if r.name=='HIS'] for c in model}
 row={'name':name,'iptm':rec['confidence']['iptm'],'his_positions':hs,'interface_residues':{k:sorted(v) for k,v in interface.items()},'his_contacts':[dict(binder_res=k[0],binder_name=k[1],egfr_res=k[2],egfr_name=k[3],**v) for k,v in contacts.items()]}
 if len(results)<5:
  pkas={}
  for label,chains in [('complex',['A','B']),('A',['A']),('B',['B'])]:
   t=s.clone()
   for c in list(t[0]):
    if c.name not in chains:del t[0][c.name]
   p=out/(name+'_'+label+'.pdb');t.write_pdb(str(p))
   with contextlib.redirect_stdout(io.StringIO()):m=propka.run.single(str(p),write_pka=False)
   pkas[label]={(g.atom.chain_id,g.atom.res_num,g.type):g.pka_value for g in m.conformations['AVR'].groups if g.titratable and not g.exclude_cys_from_results}
  shifts=[]
  for (chain,num,typ),pka in pkas['complex'].items():
   if num not in interface[chain]:continue
   free=pkas[chain].get((chain,num,typ))
   def frac(p,pH):return 1/(1+10**(pH-p))
   shifts.append({'chain':chain,'residue':num,'type':typ,'pka_complex':pka,'pka_separated_fixed_coordinates':free,'delta_pka':pka-free if free is not None else None,'protonated_fraction_complex_6.5':frac(pka,6.5),'protonated_fraction_complex_7.4':frac(pka,7.4)})
  row['interface_pka']=shifts
 results.append(row);print(name,round(row['iptm'],3),'Hiscontacts',len(contacts),'Hisinterface',[(k,len(set(hs[k])&interface[k])) for k in hs],flush=True)
 (out/'analysis.json').write_text(json.dumps(results,indent=2))
