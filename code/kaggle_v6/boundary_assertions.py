"""Prepared runtime hook. No model instantiation or network; caller supplies model/batch.
Call after loading verified checkpoint and before Trainer.predict. Never edits features.
"""
import numpy as np

def check_arrays(before,after,expected_pairs,design_mask,resolved_mask):
 design_mask=np.asarray(design_mask,bool);resolved_mask=np.asarray(resolved_mask,bool)
 assert len(design_mask)==len(resolved_mask)
 for k in ['token_bonds','type_bonds','template_mask','template_mask_frame','template_mask_cb','design_mask','chain_design_mask','msa']:
  assert np.array_equal(before[k],after[k]),f'Model boundary changed {k}'
 b=np.asarray(after['token_bonds']).squeeze();assert b.shape==(len(design_mask),len(design_mask))
 assert np.array_equal(b,b.T) and np.all(np.diag(b)==0)
 assert np.count_nonzero(b)//2==expected_pairs
 assert np.array_equal(np.asarray(after['design_mask']).squeeze().astype(bool),design_mask)
 assert np.array_equal(np.asarray(after['chain_design_mask']).squeeze().astype(bool),design_mask)
 assert np.array_equal(np.asarray(after['template_mask']).squeeze().astype(bool),~design_mask)
 missing=(~resolved_mask)&(~design_mask)
 for k in ['template_mask_frame','template_mask_cb']:assert not np.asarray(after[k]).squeeze()[missing].any()
 assert np.asarray(after['msa']).shape[-2:]==(1,len(design_mask))
 return {'postmask_assertions_pass':True,'covalent_pairs':expected_pairs,'tokens':len(design_mask),'prediction_executed':False}

def inspect_actual_model_boundary(model,batch,metadata,mode):
 assert model.masker.mask is False,'Unexpected real masker setting: stop before prediction'
 assert model.bond_type_feature is True,'Unexpected bond type setting'
 assert model.use_templates is True,'Unexpected template setting'
 keys=['token_bonds','type_bonds','template_mask','template_mask_frame','template_mask_cb','design_mask','chain_design_mask','msa']
 before={k:batch[k].detach().cpu().numpy().copy()for k in keys}
 masked=model.masker(batch)
 after={k:masked[k].detach().cpu().numpy().copy()for k in keys}
 return check_arrays(before,after,28 if mode else 0,metadata['design_mask'],metadata['token_resolved_mask'])

if __name__=='__main__':
 # Tiny synthetic regression fixture: only verifies assertions, not protein features.
 n=6;dm=np.array([1,1,0,0,0,0],bool);rm=np.array([1,1,1,1,0,1],bool);b=np.zeros((n,n));b[0,1]=b[1,0]=1
 x={'token_bonds':b,'type_bonds':b.astype(int),'template_mask':~dm,'template_mask_frame':(~dm)&rm,'template_mask_cb':(~dm)&rm,'design_mask':dm,'chain_design_mask':dm,'msa':np.zeros((1,n),int)}
 assert check_arrays(x,x,1,dm,rm)['postmask_assertions_pass']
 y={k:v.copy()for k,v in x.items()};y['token_bonds']*=0
 try:check_arrays(x,y,1,dm,rm)
 except AssertionError:pass
 else:raise AssertionError('Failed to reject changed bonds')
 y={k:v.copy()for k,v in x.items()};y['template_mask_frame'][4]=True
 try:check_arrays(y,y,1,dm,rm)
 except AssertionError:pass
 else:raise AssertionError('Failed to reject invented missing geometry')
 print('Synthetic assertion regression passed. Actual model boundary not executed.')
