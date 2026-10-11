import pathlib,json,numpy as np
from Bio.PDB import PDBParser
P=PDBParser(QUIET=True);base=pathlib.Path('/tmp/tnf');full=P.get_structure('f',base/'1tnf.pdb')[0];crop=P.get_structure('c',base/'regenerated/tnf_crop.pdb')[0]
cr=[r for c in crop for r in c if 'CA' in r];ref=np.array([r['CA'].coord for r in cr]);keep={(r.parent.id,r.id) for r in cr}
omitted=[a for c in full for r in c if (c.id,r.id) not in keep for a in r if a.element not in ('H','D')]
records=[]
for name in ['unmodified_control','mpnn_repair_03','mpnn_repair_02']:
 m=P.get_structure(name,base/f'compare_results/output/{name}_complex.pdb')[0];target=list(m)[0];binder=list(m)[-1];x=np.array([r['CA'].coord for r in target if 'CA' in r]);assert len(x)==len(ref)
 xc=x-x.mean(0);yc=ref-ref.mean(0);u,_,v=np.linalg.svd(xc.T@yc);rot=u@np.diag([1,1,np.linalg.det(u@v)])@v;shift=ref.mean(0)-x.mean(0)@rot
 rmsd=np.sqrt(np.mean(np.sum((x@rot+shift-ref)**2,axis=1)))
 ba=[a for a in binder.get_atoms() if a.element not in ('H','D')];coords=np.array([a.coord for a in ba])@rot+shift
 dist=np.linalg.norm(coords[:,None]-np.array([a.coord for a in omitted])[None],axis=-1);pairs=np.argwhere(dist<2.0)
 unique_b=set(int(i) for i,j in pairs);res=set(ba[i].parent.id[1] for i in unique_b)
 nearest=np.unravel_index(np.argmin(dist),dist.shape)
 r=dict(name=name,target_CA_fit_RMSD_A=float(rmsd),omitted_atom_pairs_below2A=len(pairs),binder_residues_with_omitted_clash=sorted(res),min_distance_A=float(dist[nearest]),nearest_binder_res=ba[nearest[0]].parent.id[1],nearest_target_chain=omitted[nearest[1]].parent.parent.id,nearest_target_res=omitted[nearest[1]].parent.id[1])
 records.append(r)
 print(r)
json.dump(dict(method='Fit predicted crop target CA to original crop CA ordered by chain/residue; inspect binder heavy atoms vs omitted full-trimer atoms. <2 A is geometric overlap screen, not relaxed energy.',source='https://files.rcsb.org/download/1TNF.pdb',results=records),open(base/'full_trimer_clash_audit.json','w'),indent=2)
