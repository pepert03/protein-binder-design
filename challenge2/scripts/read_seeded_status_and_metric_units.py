import pathlib, subprocess, inspect, json, re
print('PROCESS', subprocess.getoutput("ps -eo pid,etime,args | grep '[p]ython -u run.py'"))
p=pathlib.Path('/kaggle/working/tnf/run.log')
print('LOG',p.read_text()[-18000:] if p.exists() else 'NO_LOG')
from colabdesign.af.loss import get_pae_loss
print('INSTALLED_PAE_SOURCE',inspect.getsource(get_pae_loss))
r=pathlib.Path('/kaggle/working/tnf/run.py').read_text()
print('SCRIPT_TAIL',r[-6500:])
print('FILES',[(x.name,x.stat().st_size) for x in pathlib.Path('/kaggle/working/tnf').glob('*') if x.suffix in ('.pdb','.json')])
