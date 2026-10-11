# Pure binder-intra PAE callback prototype,11October21:08

Version1.1.1 model source accepts loss_callback(inputs,outputs) and merges its returned losses. Added callback uses only binder-binder pairs, all valid seq_mask positions, symmetric PAE normalized by31. It follows the versioned expected-PAE bin-center definition.

CPU JAX0.4.35 test passed: equal logits over5/15/25Å bins gives15/31, target/interface logits have zero gradient, binder-binder logits have nonzero gradient, jit matches. This tests the differentiable mathematical term, NOT AF inference, API wiring onDGX/JAX0.6.2, parameter compatibility, memory fit or design quality.

Proposed wiring only: make callback with known target_length before mk_afdesign_model(loss_callback=callback), then explicitly include binder_intra_pae weight. Versioned model assigns weight1 to callback losses absent from opt weights, so omission is dangerous. Baseline must set it0; intervention must record its chosen weight. Refuse target_length mismatch and empty binder in runner preflight. Do not launch based on this note.

This callback does not change stock pae/i_pae semantics or reconstruct historical raw PAE. All reporting should continue explicit scopes. NopH mechanism or affinity inference.

Source read:
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/model.py
https://github.com/sokrypton/ColabDesign/blob/v1.1.1/colabdesign/af/loss.py
