import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import jax,jax.numpy as jnp,numpy as np
from folding_loss_callback import make_binder_intra_callback
fn=make_binder_intra_callback(2);base=jnp.zeros((4,4,3));inputs={'seq_mask':jnp.ones(4)}
def loss(x):return fn(inputs,{'predicted_aligned_error':{'logits':x,'breaks':jnp.array([0.,10.])}})['binder_intra_pae']
value=loss(base);g=jax.grad(loss)(base)
assert np.isclose(float(value),15/31)
assert np.allclose(np.asarray(g)[:2],0) and np.allclose(np.asarray(g)[:,:2],0)
assert np.max(np.abs(np.asarray(g)[2:,2:]))>0
assert np.isclose(float(jax.jit(loss)(base)),15/31)
print('PASS CPU-JAX differentiable callback: expected units, target/interface zero-gradient, binder gradient, jit')
