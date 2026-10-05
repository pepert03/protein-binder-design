import os,subprocess,urllib.request,json,time,hashlib,zipfile,pathlib,sys
root=pathlib.Path('/kaggle/working/egfr');root.mkdir(exist_ok=True);os.chdir(root)
def stage(name,cmd):
 import threading,queue
 print('STAGE_START',name,time.time(),flush=True)
 q=queue.Queue()
 with open(name+'.log','w') as f:
  proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
  def reader():
   for line in proc.stdout: q.put(line)
   q.put(None)
  worker=threading.Thread(target=reader,daemon=True);worker.start()
  last=time.monotonic();done=False
  while not done:
   try:
    line=q.get(timeout=5)
    if line is None:done=True
    else:f.write(line);f.flush();print(line,end='',flush=True)
   except queue.Empty:pass
   if time.monotonic()-last>=300:
    print('HEARTBEAT',name,'pid',proc.pid,'poll',proc.poll(),'time',time.time(),flush=True)
    diag=subprocess.run(['nvidia-smi','--query-gpu=memory.used,utilization.gpu','--format=csv,noheader'],capture_output=True,text=True,timeout=30)
    print('GPU_SNAPSHOT',diag.returncode,diag.stdout,diag.stderr,flush=True)
    files=[p for p in root.rglob('*') if p.is_file()]
    print('OUTPUT_SNAPSHOT','cif',sum(p.suffix=='.cif' for p in files),'csv',sum(p.suffix=='.csv' for p in files),'bytes',sum(p.stat().st_size for p in files),flush=True)
    pathlib.Path('heartbeat.json').write_text(json.dumps({'stage':name,'pid':proc.pid,'time':time.time(),'files':len(files)}))
    last=time.monotonic()
  rc=proc.wait();worker.join(timeout=5)
 print('STAGE_EXIT',name,rc,flush=True)
 return rc

try:
 if stage('install',[sys.executable,'-m','pip','install','-q','boltzgen'])!=0:raise RuntimeError('installation failed')
 stage('uninstall',[sys.executable,'-m','pip','uninstall','-y','torchvision','torchaudio','torchcodec'])
 urllib.request.urlretrieve('https://files.rcsb.org/download/6ARU.pdb','6ARU.pdb')
 pathlib.Path('chainA.pdb').write_text(''.join(l for l in open('6ARU.pdb') if l.startswith('ATOM') and l[21]=='A')+'END\n')
 # Objetivo: parte DIII del parche de novo 396 a 479 UniProt (17 residuos).
 # label_seq_id = resnum 6ARU - 3 = UniProt - 27, verificado residuo a residuo (code 126).
 LABELS='369,395,397,400,402,403,406,407,408,427,428,431,446,447,449,451,452'
 pathlib.Path('patchDIII.yaml').write_text("entities:\n  - protein:\n      id: B\n      sequence: 50..90\n  - file:\n      path: chainA.pdb\n      include:\n        - chain:\n            id: A\n      binding_types:\n        - chain:\n            id: A\n            binding: "+LABELS+"\n")
 import importlib.util
 p=pathlib.Path(importlib.util.find_spec('boltzgen').origin).parent
 for n in ['task/predict/data_from_yaml.py','task/predict/data_from_generated.py']:
  s=(p/n).read_text();assert 'msa={}' in s and 'msa_server' not in s and 'mmseqs' not in s
  print('LOCAL_MSA_EMPTY_VERIFIED',n,flush=True)
 if stage('check',['boltzgen','check','patchDIII.yaml'])!=0:raise RuntimeError('invalid specification')
 rc=stage('runDIII',['boltzgen','run','patchDIII.yaml','--output','oc_DIII','--protocol','protein-anything','--num_designs','16','--budget','8','--config','filtering','filter_bindingsite=true'])
 pathlib.Path('runDIII.exit').write_text(str(rc))
 print('RUN_FINISHED',rc,flush=True)
 pathlib.Path('post.py').write_text('''import pandas as pd,glob,gemmi,numpy as np,os,hashlib
D='/kaggle/working/egfr/oc_DIII/final_ranked_designs/'
P=[369,395,397,400,402,403,406,407,408,427,428,431,446,447,449,451,452]
out=open('/kaggle/working/egfr/post_DIII.txt','w')
def p(*a): print(*a,file=out,flush=True)
a=pd.read_csv(D+'all_designs_metrics.csv')
cols=[c for c in ['id','design_to_target_iptm','min_design_to_target_pae','filter_rmsd','designfolding-filter_rmsd','liability_score','num_design','final_rank','pass_filters'] if c in a.columns]
p('M>> ALL'); p(a[cols].round(3).to_string())
for fn in sorted(glob.glob(D+'final_8_designs/*.cif')):
    m=gemmi.read_structure(fn)[0];A=m[0];B=m[1]
    p('CHAINS>>',os.path.basename(fn),[c.name for c in m])
    seq=gemmi.one_letter_code([r.name for r in A])
    xa=np.array([[x.pos.x,x.pos.y,x.pos.z] for r in A for x in r if x.element.name!='H'])
    res=[]
    for r in B:
        xb=np.array([[x.pos.x,x.pos.y,x.pos.z] for x in r if x.element.name!='H'])
        res.append((r.seqid.num,float(np.sqrt(((xb[:,None]-xa[None])**2).sum(-1)).min())))
    c4=[n for n,d in res if d<4]
    p('M>>',os.path.basename(fn),len(seq),'| contactos<4A target',len(c4),'en parche',[n for n in c4 if n in P],'| min dist a parche',round(min(d for n,d in res if n in P),1))
    pf=a.loc[a['id']==os.path.basename(fn).split('_',1)[1].replace('.cif',''),'pass_filters'].tolist()
    p('PASS>>',os.path.basename(fn),pf)
    bb=open(fn,'rb').read()
    p('SHA>>',os.path.basename(fn),len(bb),hashlib.sha256(bb).hexdigest())
p('DONE')
''')
 stage('post',[sys.executable,'post.py'])
 print(pathlib.Path('post_DIII.txt').read_text()[-6000:],flush=True)
except Exception as e:
 pathlib.Path('failure.txt').write_text(repr(e));print('FAILED',repr(e),flush=True)
finally:
 files=sorted(p for p in root.rglob('*') if p.is_file() and p.name!='manifest.json')
 manifest=[{'path':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
 pathlib.Path('manifest.json').write_text(json.dumps(manifest,indent=2))
 with zipfile.ZipFile('/kaggle/working/DIII_all_outputs.zip','w',zipfile.ZIP_DEFLATED) as z:
  for p in root.rglob('*'):
   if p.is_file():z.write(p,p.relative_to(root))
 print('ARCHIVE_BYTES',os.path.getsize('/kaggle/working/DIII_all_outputs.zip'),'FILES',len(manifest),flush=True)
 print('ALL_OUTPUTS_PRESERVED',flush=True)