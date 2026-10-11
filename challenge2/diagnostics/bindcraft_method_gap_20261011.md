# Method gap check,11October19:05

Freshly inspected BindCraft maintainer documentation supports changing the design protocol, not simply running the same366 trajectory longer.

- BindCraft designs through AF2 multimer and swaps five model weights, then uses monomer complex and binder-alone predictions as separate self-consistency checks. Our memory-limited fresh run used one monomer model and zero recycles for gradient design.
- Documentation notes relaxed-to-discrete confidence collapse can indicate adversarial or non-globular designs. Our140-step trajectory was weak throughout, so this explanation is plausible context, not a diagnosis of its cause.
- Maintainers explicitly say current scores are not affinity predictors. Their ranking should not be translated into pH response or mouse affinity.
- Their hardware reference says32GBGPU fits about550 total residues,80GB about950. Full456+50 is506 residues, so a32GB class GPU is a documented plausible starting point, not a verified available resource or fit guarantee. The actual university slice is16GB.40GB host-RAM recommendation is a different quantity.
- Their typical globular binder range is60-180aa. A50aa run is not prohibited, but reducing binder length solely to fit memory changes the design space. Their optimal-size recommendation does not guarantee TNF success.
- PyRosetta licensing and environment requirements apply. We have not installed or run full BindCraft and have not verified a permitted free larger-memory resource. No admin request or paid resource proposed.

A bounded next pilot could compare current monomer baseline with a folding-focused objective at identical target/crop/seed and metric masks. Multimer model diversity is a separate resource-dependent branch. Neither should be promised to deliver by the competition deadline. First retain a stable folded candidate and intact-trimer interface before histidine engineering.

Sources fetched11October2026:
https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft
https://github.com/martinpacesa/BindCraft/blob/main/README.md
