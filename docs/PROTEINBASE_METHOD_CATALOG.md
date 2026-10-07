# Proteinbase design method options supplied by Jose

Recorded from the options Jose saw in the submission form on 7 October 2026. This is a reference catalogue, not a list of methods used in this EGFR batch. Names, authors and descriptions below are supplied form text, not independently verified scientific claims. Truncated URLs and descriptions are left truncated rather than guessed.

## Listed methods

AntiBERTY; BindCraft; BLOSUM62; BoltzGen; co-TIMED; DSM Synteract; ESM2 Optimization; Evo Diff; Foundry; Llama 3.1 405B; Mosaic; PPIFlow; Pro-1; ProTrek; ProtRL; PXDesign; RFdiffusion; RFpeptides.

## Supplied form text

The standalone letters and prefixes in this copy appear as pasted by Jose; they are not treated as additional method names.

```text
Los que me salen son:
AntiBERTY
by Jeff Ruffolo Antibody-specific transformer language model pre-trained on 558M natural antibody sequences.
BindCraft
by EPFL/MIT Protein binder design pipeline that can be used to design miniprotein and peptide binders
BLOSUM62
by agitter BLOSUM (BLOcks Substitution Matrix) matrix is a substitution matrix used for sequence alignment of proteins
M
BoltzGen
by Massachusetts Institute of Technology All-atom diffusion-based generative model developed for de novo design of protein and peptide binders across a wide array of target types
co-TIMED
by Mert Unal (Wells Wood Lab). A workflow for designing protein binders using co-TIMED, a method for co-designing sequences in a protein complex
DSM Synteract
by Gleghorn Lab X Synthyra pLM trained with masked diffusion to enable both high-quality representation learning and generative protein design
ESM2 Optimization
by Meta. ESM2 is a masked language model trained on UniRef protein sequences. Used for structure prediction, property prediction, functional annotation, and clustering similar proteins.
Evo Diff
by Microsoft Research Motif-scaffolding pipeline that designs protein sequences around a given motif using the EvoDiff model.
w Foundry
by Baker Lab Foundry provides tooling and infrastructure for using and training all classes of models for protein design, including design (RFD3), inverse folding (ProteinMPNN) and protein folding (RF3).
https://github.com/RosettaCommons/foundry
Llama 3.1 405B
by Meta Llama 3 is a family of multilingual large language models (LLMs) which are pretrained and instruction tuned generative models
Mosaic
by Escalante Bio Open-source framework for combining and composing a variety of structure prediction and other models for design, filtering, or guidance.
PPIFlow
by Mingchen Chen PPIFlow demonstrated a flow-matching-based framework and integrated design workflow for the de novo generation of high-affinity binders targeting precise epitopes. It supports a diverse range of tasks including protein-protein interaction (...
Pro-1
by Michael Hla. A reasoning model that proposes modifications to improve protein stability.
ProTrek
by Westlake Representation Learning Lab ProTrek is a tri-modal protein language model that jointly models protein sequence, structure and function (SSF)
E
ProtRL
by CRG Ferruz Lab Fine-tuning pipeline combining reinforcement learning algorithms (GRPO, DPO) with the ZymCTRL protein language model.
(PXDesign
by ByteDance PXDesign is a model suite for de novo protein-binder design a diffusion generator (PXDesign-d) paired with Protenix and AF2-IG confidence models for selection. Across seven targets, PXDesign delivers 17-82% nanomolar hits on six. https://...
RFdiffusion
by Baker Lab. A flexible de novo binder design workflow using RFdiffusion for backbone generation
RFpeptides
by Baker Lab Extended RFdiffusion that incorporates cyclic relative positional encodings to enable the generation of the macrocyclic peptide backbones. https://www.nature.com/articles/s41589-025-01929-w
https://github.com/RosettaCommons/RFdiffusion/blob...

Y tambien deja poner Create new
Method name
e.g. RFDiffusion + AlphaFold2 refinement
Description
Explain the general pipeline or reasoning behind this method.

Configuration (optional)
Explain the specific workflow, prompts, parameters, or settings used to produce the selected proteins.
```

## Relation to the current EGFR work

The recorded campaign uses BoltzGen for initial candidate generation, ProtonPottsMPNN for sequence redesign, and Boltz 2 for structural screening. The original C5 generation manifest remains incomplete. The presence of a tool in this catalogue does not mean it was run, tested, adopted or validated. The form also allows a new method entry with a name, description and optional configuration. No new computation or submission was performed when recording this catalogue.
