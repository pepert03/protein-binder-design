import gemmi,numpy as np,pathlib,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
s=gemmi.read_structure('/tmp/1IVO.pdb');c=s[0]['A'];rs={r.seqid.num+24:r for r in c if r.het_flag=='A'};g=np.array([[a.pos.x,a.pos.y,a.pos.z] for ch in s[0] if ch.name in ['A','E'] for r in ch if r.name=='NAG' for a in r])
fig,axs=plt.subplots(2,2,figsize=(9,9)); names=['r03_l07','n03_rank03_native']
for k,n in enumerate(names):
 st=gemmi.read_structure(str(next(pathlib.Path('/tmp/egfr20/predictions/'+n).rglob('*.cif'))))[0];nums=[i for i in range(335,506) if i in rs];x=np.array([[st['B'][i-25].find_atom('CA','*').pos[j] for j in range(3)] for i in nums]);y=np.array([[rs[i].find_atom('CA','*').pos[j] for j in range(3)] for i in nums]);xc=x.mean(0);yc=y.mean(0);u,ss,v=np.linalg.svd((x-xc).T@(y-yc));d=np.eye(3);d[-1,-1]=np.linalg.det(u@v);R=u@d@v;t=yc-xc@R;A=np.array([[r.find_atom('CA','*').pos[j] for j in range(3)] for r in st['A']])@R+t;Y=np.array([[r.find_atom('CA','*').pos[j] for j in range(3)] for r in c if r.het_flag=='A'])
 for z,(i,j) in enumerate([(0,1),(0,2)]):
  ax=axs[k,z];ax.plot(Y[:,i],Y[:,j],color='#c0bebf',lw=1,label='EGFR experimental');ax.plot(A[:,i],A[:,j],color='#73a89a',lw=2.5,label='Binder predicted');ax.scatter(g[:,i],g[:,j],color='#ed313e',s=9,label='Resolved sugar atoms');ax.set_aspect('equal');ax.set_title(n+(' / XY' if z==0 else ' / XZ'),fontsize=12);ax.set_xlabel('X (A)');ax.set_ylabel(('Y' if z==0 else 'Z')+' (A)');ax.grid(alpha=.15)
axs[0,0].legend(fontsize=8);fig.suptitle('1IVO context after domain III alignment',fontsize=18);fig.text(.08,.01,'Partial experimental glycans, one predicted pose. Geometric flags are not a binding assay.',fontsize=10);fig.tight_layout(rect=[0,.04,1,.96]);fig.savefig('/downloads/egfr-context.png',dpi=170)
