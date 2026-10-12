"""ColabDesign1.1.1-compatible optional pure binder-intra PAE callback.
Not a complete runner; must be explicitly weighted and GPU validated.
"""
import jax
import jax.numpy as jnp

def expected_pae(outputs):
    head=outputs['predicted_aligned_error'];breaks=head['breaks'];step=breaks[1]-breaks[0]
    centers=jnp.append(breaks+step/2,breaks[-1]+step*1.5)
    return (jax.nn.softmax(head['logits'],axis=-1)*centers).sum(-1)

def make_binder_intra_callback(target_length):
    if not isinstance(target_length,int) or target_length<=0:raise ValueError('Positive target length required')
    def binder_intra(inputs,outputs):
        p=expected_pae(outputs)/31.;p=(p+p.T)/2
        mask=inputs['seq_mask'];binder=(jnp.arange(mask.shape[0])>=target_length)*mask
        pair=binder[:,None]*binder[None,:]
        return {'binder_intra_pae':jnp.sum(p*pair)/(jnp.sum(pair)+1e-8)}
    return binder_intra
