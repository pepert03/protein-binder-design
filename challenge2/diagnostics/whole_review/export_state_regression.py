"""CPU mock regression: catches current-prediction vs saved-best export ambiguity."""
from pathlib import Path
class MockAF:
    def __init__(self):self.best='OLD_DESIGN';self.current='NEW_PREDICTION'
    def save_pdb(self,filename,get_best=True):return self.best if get_best else self.current
m=MockAF()
assert m.save_pdb('tmp')=='OLD_DESIGN'
assert m.save_pdb('tmp',get_best=False)=='NEW_PREDICTION'
root=Path(__file__).resolve().parents[2]/'scripts'
for name in ('run3.py','run5.py'):
    assert "af.save_pdb('tmp.pdb', get_best=False)" in (root/name).read_text()
assert 'model.save_pdb(str(fn), get_best=False)' in (root/'repredict_repair.py').read_text()
print('PASS mock export-state regression; AF inference not run')
