import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from his_contacts import model_counts,conservative_count

def atom(n,name,res,chain,r,x):
 return f'ATOM  {n:5d} {name:^4s} {res:3s} {chain}{r:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{80.:6.2f}           N\n'
def model(i,target,binder):
 return f'MODEL     {i:4d}\n'+atom(1,'NZ','LYS','A',1,target)+atom(2,'ND1','HIS','B',1,binder)+atom(3,'NE2','HIS','B',1,binder+1)+'ENDMDL\n'
# Same-model separation100Å, but target1 near binder2 and vice versa.
text=model(1,0,100)+model(2,100,0)+'END\n'
rows=model_counts(text);assert [r['his_count_within_cutoff'] for r in rows]==[0,0]
assert conservative_count(text,1)==0
# One model has a real contact; other does not. Conservative summary must be0.
text=model(1,0,2)+model(2,0,100)+'END\n'
assert [r['his_count_within_cutoff'] for r in model_counts(text)]==[1,0]
assert conservative_count(text,1)==0
assert conservative_count(model(1,0,2)+'END\n',1)==1
try:conservative_count(text,2)
except ValueError:pass
else:raise AssertionError('Length guard')
print('PASS: cross-model false contacts blocked, min-model rule and length guard')
