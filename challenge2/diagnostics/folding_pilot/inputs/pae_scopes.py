"""Explicit PAE scopes from one model's raw Angstrom matrix. No binding claims."""
import numpy as np

def pae_scopes(raw_pae_A,target_length,hotspot_indices=()):
    p=np.asarray(raw_pae_A,dtype=float)
    if p.ndim!=2 or p.shape[0]!=p.shape[1] or not np.isfinite(p).all():raise ValueError('Finite square per-model PAE required')
    if not 0<target_length<len(p):raise ValueError('Nonempty target and binder required')
    hs=list(hotspot_indices)
    if any(not isinstance(i,(int,np.integer)) or i<0 or i>=target_length for i in hs):raise ValueError('Hotspots must be zero-based target indices')
    if len(set(hs))!=len(hs):raise ValueError('Duplicate hotspot indices')
    sym=(p+p.T)/2;b=np.arange(target_length,len(p));t=np.arange(target_length)
    scopes={'binder_intra':sym[np.ix_(b,b)],'binder_to_entire_target':sym[np.ix_(b,t)]}
    if hs:scopes['binder_to_hotspots']=sym[np.ix_(b,hs)]
    return {k:{'mean_A':float(v.mean()),'normalized_by31':float(v.mean()/31),'pair_count':int(v.size),'symmetrized':True} for k,v in scopes.items()}
