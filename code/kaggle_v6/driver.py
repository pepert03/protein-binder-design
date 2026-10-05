"""Corrected public reference control only; four isolated arms, one sample each."""
from pathlib import Path
import json,hashlib,importlib.util,gc,shutil,os,time
import numpy as np,torch
from hydra.utils import instantiate
from omegaconf import OmegaConf
from huggingface_hub import hf_hub_download
from boltzgen.task.predict.data_from_generated import collate
from boundary_assertions import inspect_actual_model_boundary
from importlib.metadata import version
root=Path('/kaggle/working/corrected_controls');inputs=Path('/kaggle/working/corrected_inputs');root.mkdir(exist_ok=True)
assert version('boltzgen')=='0.3.2'
pkg=Path(importlib.util.find_spec('boltzgen').origin).parent
expected=json.loads((inputs/'expected_source_hashes.json').read_text())
for n,h in expected.items():assert hashlib.sha256((pkg/n).read_bytes()).hexdigest()==h,n
assert 'msa={}'in(pkg/'task/predict/data_from_generated.py').read_text()
assert torch.cuda.is_available();print('GPU_RESOURCES',torch.cuda.get_device_name(0),torch.cuda.mem_get_info(),flush=True)
moldir=hf_hub_download(repo_id='boltzgen/inference-data',filename='mols.zip',repo_type='dataset',revision='c3d36fd276e9caf098c75d4113c6d5eb320b1a4c')
checkpoint=hf_hub_download(repo_id='boltzgen/boltzgen-1',filename='boltz2_conf_final.ckpt',revision='c1be29e1f82ffcc72264f64b993c43fb4e0d17f0')
(root/'runtime_provenance.json').write_text(json.dumps({'checkpoint':checkpoint,'source_hashes':expected,'boltzgen':version('boltzgen'),'torch':version('torch'),'samples_per_arm':1,'original_filters_changed':False},indent=2))
shutil.copy2(inputs/'driver.py',root/'driver.py');shutil.copy2(inputs/'boundary_assertions.py',root/'boundary_assertions.py')
rows=[]
# Fab first closes the actual resource blocker. No shortened chains.
for name in ['Fab_6ARU','EGF_1NQL']:
 for mode in [False,True]:
  if mode and any(r['id']==name and not r['prediction_complete'] for r in rows):
   rows.append({'id':name,'disulfide_on':mode,'samples':1,'skipped':True,'reason':'Earlier arm failed; no repeat resource configuration','model_execution_started':False});continue
  arm=f'{name}_SS_{str(mode).lower()}';d=root/arm;d.mkdir(exist_ok=True);task=None
  for ext in ['.cif','.npz']:shutil.copy2(inputs/(name+ext),d/(name+ext))
  cfg=OmegaConf.load(inputs/f'fold_one_sample_SS_{str(mode).lower()}_PREPARED.yaml');cfg.data.design_dir=str(d);cfg.output=str(d);cfg.data.cfg.moldir=moldir;cfg.checkpoint=checkpoint
  OmegaConf.save(cfg,d/'resolved_folding.yaml');r={'arm':arm,'id':name,'disulfide_on':mode,'samples':1,'model_execution_started':False,'prediction_complete':False,'biological_conclusion':None};rows.append(r)
  try:
   task=instantiate(cfg);assert len(task.data.predict_set)==1
   batch=collate([task.data.predict_set[0]]);assert batch['id']==[name]
   # Instantiate exact checkpoint and trainer only; not a prediction call.
   task.run(cfg,run_prediction=False)
   meta=np.load(d/(name+'.npz'));boundary=inspect_actual_model_boundary(task.model_module,batch,meta,mode)
   (d/'actual_model_boundary.json').write_text(json.dumps(boundary,indent=2));r['actual_model_boundary_pass']=True
   print('ACTUAL_MODEL_BOUNDARY_PASS',arm,boundary,flush=True)
   r['model_execution_started']=True;torch.cuda.reset_peak_memory_stats()
   task.lightning_trainer.predict(task.model_module,datamodule=task.data,return_predictions=False)
   r['cuda_peak_allocated']=torch.cuda.max_memory_allocated();r['cuda_peak_reserved']=torch.cuda.max_memory_reserved()
   zpath=d/'fold_out_npz'/f'{name}.npz';assert zpath.exists()
   with np.load(zpath)as z:
    assert all(k in z for k in ['coords','iptm','ptm']);assert z['coords'].shape[0]==1
    assert np.isfinite(z['coords']).all()and np.isfinite(z['iptm']).all()
    r['iptm']=z['iptm'].tolist();r['ptm']=z['ptm'].tolist()
   r['prediction_complete']=True;print('CORRECTED_ONE_SAMPLE_FINISHED',r,flush=True)
  except Exception as e:
   r['failure']=repr(e);(d/'failure.txt').write_text(repr(e));print('ARM_FAILED',arm,repr(e),flush=True)
   # Stop the workflow if actual boundary failed; never override model parameters.
   if not r.get('actual_model_boundary_pass'):
    (root/'arm_results.json').write_text(json.dumps(rows,indent=2));raise
  finally:
   (root/'arm_results.json').write_text(json.dumps(rows,indent=2))
   del task;gc.collect();torch.cuda.empty_cache()
print('FOUR_ARM_WORKFLOW_ENDED',json.dumps(rows),flush=True)
