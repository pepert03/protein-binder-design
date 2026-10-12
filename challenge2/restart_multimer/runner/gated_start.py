"""Local process wait, no account/site polling. Start cohort only on both pass files."""
import pathlib,json,time,subprocess
started=time.time()
while time.time()-started<1200:
 files=[pathlib.Path('preflight_recycles'+str(n)+'.json') for n in [0,1]]
 if all(p.exists() for p in files) and 'PREFLIGHT_COMPLETE' in pathlib.Path('preflight.log').read_text():
  data=[json.loads(p.read_text()) for p in files]
  for d in data:
   assert d['models']==['model_1_multimer_v3']
   assert d['target_len']==277 and d['binder_len']==60
   assert len(set(d['chain_ids_seen']['asym_id']))==4
   assert d['elapsed_s']>0
  assert 'PREFLIGHT_COMPLETE' in pathlib.Path('preflight.log').read_text()
  print('GATES_PASSED_START_COHORT',flush=True)
  raise SystemExit(subprocess.call(['timeout','10000','../venv/bin/python','cohort.py']))
 text=pathlib.Path('preflight.log').read_text()
 if 'Traceback' in text or 'RESOURCE_EXHAUSTED' in text:
  print('PREFLIGHT_FAILED_NO_COHORT',flush=True);raise SystemExit(1)
 time.sleep(10)
print('PREFLIGHT_TIMEOUT_NO_COHORT',flush=True)
