# Proton-PottsMPNN: technical brief for the EGFR challenge

Date: 2026-10-05. Source: repo https://github.com/christian-creator/ProtonPottsMPNN
(README + code reviewed; the bioRxiv PDF could not be fetched, so this brief is
code-grounded, not paper-grounded). Paper: Jacobsen et al. 2026, bioRxiv
10.64898/2026.09.30.755438. Extends PottsMPNN (Birnbaum & Keating, PNAS 2026).

## What it is

An inverse-folding (sequence-design) model for **pH-switchable binders**. It takes
PottsMPNN — a Potts-energy model (single-site fields + pairwise couplings) over
backbone geometry — and adds an **explicit protonation-state alphabet** (30-token
"v6" vocab):

- `HIS-P` = protonated histidine (+1) vs `HIS-S` = neutral histidine
- `ASP-P` / `GLU-P` = protonated acids (neutral COOH) vs `ASP-D` / `GLU-D` = deprotonated (-1)

Because the learned Potts energy is protonation-aware, one energy function scores
alternative protonation assignments on a fixed backbone. That is what makes
pH-conditioned design possible without molecular dynamics.

## How pH dependence is modeled

1. **Label**: a FLAML labeller (EV6, shipped in `transforms/ev6/`) assigns a
   protonation state to every titratable residue of the input complex (binder AND
   target; target states are then fixed context).
2. **Pin + redesign**: `PottsMPNNPHEngine` pins protonated centres on the binder
   (e.g. `HIS-P`, `ASP-P`, `GLU-P`) and redesigns their neighbourhood with
   block-descent to maximize the **selectivity gap**
   `E_selective = E_P - E_D` (energy with centre protonated minus deprotonated),
   trading against stability:
   `J = (1-lambda) * z(H_stab) + lambda * z(sum sel)` (manuscript Eq. 6).
3. **Pareto front**: sweeping `combined_lambda` 0 (pure stability) -> 1 (pure
   selectivity) yields designs in (stability, selectivity) space; the notebook
   marks the Pareto front and exports the top designs to fold.

Validated by the authors on PD-L1 (yeast display, pH 5.0 vs 7.4, flow cytometry).

## Inputs / outputs

- **Input**: a binder-target complex backbone (PDB/CIF or AtomArray) + binder chain
  ID. The backbone is fixed; only sequence (and protonation states) change.
- **Outputs** (`inference/outputs/`): `designs.fasta` (canonical 1-letter, foldable),
  `designs_states.fasta` (3-letter + protonation states, e.g. `... ASP-P ... HIS-P ...`),
  `designs.tsv`/`.json` (energies, pinned centres), `trajectory.tsv` (every
  optimisation step), `pareto_front.png`, `pareto_fold_manifest.json` (designs
  selected for folding).

## Compute requirements

- **Design is CPU-only.** Checkpoint is 21 MB
  (`checkpoints/potts_v6_afdb_edge_his0.3_acid0.06`); model loads with
  `map_location="cpu"`; one forward pass per backbone, then block-descent on CPU.
  `n_jobs > 1` fans independent (criteria x seed) solves across a fork pool.
- **No GPU needed for design.** GPU is only needed if you use their packaged RF3
  folder (weights ~3 GB **not shipped**) — we skip that and fold with Boltz-2,
  which already works on the DGX (~3 min per EGF-sized fold).
- **HBPLUS**: the *external C binary* is needed only for the standalone labeller
  script and fold scoring (`scoring/`). The design path uses foundry's built-in
  H-bond calculator, so **design runs without it**.
- **Install**: needs Python 3.12 (`uv python install 3.12` works anywhere),
  `uv venv`, `uv pip install -e ./foundry` (torch, lightning, atomworks).
  Note: `install.sh` is **empty (0 bytes)** in the current snapshot — install
  manually with uv; do not rely on it.

## Install and run (exact commands)

```bash
cd /tmp && git clone https://github.com/christian-creator/ProtonPottsMPNN
cd ProtonPottsMPNN
uv python install 3.12
uv venv --python 3.12 --clear
source .venv/bin/activate
uv pip install -e ./foundry          # NOT ./install.sh (empty); add extras only if needed
python -c "import mpnn; print('ok')"
```

```python
from mpnn.inference_engines.potts_mpnn_ph import PottsMPNNPHEngine, PHDesignCriteria

CKPT = "checkpoints/potts_v6_afdb_edge_his0.3_acid0.06"
engine = PottsMPNNPHEngine(checkpoint_path=CKPT, extended_vocab="v6")  # must be "v6"
crit = PHDesignCriteria(
    method="block_descent", backend="potts",
    combined_lambda=0.3,                 # sweep 0.0..1.0 for the Pareto front
    center_types=["HIS-P", "ASP-P", "GLU-P"],
    placement_by="scan_potts",           # automatic placement of centres
    dep_map={"HIS-P": ["HIS-S"], "ASP-P": ["ASP-D"], "GLU-P": ["GLU-D"]},
    block_size=3, temperature=0.05, neighbour_k=16, max_mutations=20,
    record_trajectory=True,
)
design_set = engine.run_ph_redesign(
    atom_array=aa,                       # binder-EGFR complex (e.g. a Boltz-2 fold of a candidate)
    binder_chain="B",                    # <- set to the binder's chain in our files
    criteria_list=[crit], seed=0, n_jobs=4)
```

There is also a placement-scan script (`inference/design_placement_scan.py`) and a
`method="converged_mcmc"` / `"autoregressive"` alternative in `PHDesignCriteria`.

## How to apply it to our EGFR case

**Path A — rescore existing candidates (cheap, do first).**
Featurize the 7 passing candidates' binder-EGFR complexes with
`prepare_potts_input(..., extended_vocab="v6")` and run the placement scan /
selectivity-gap scoring. This ranks which backbones already sit near a
pH-switchable solution and where protonated centres want to go. CPU-only,
minutes per backbone.

**Path B — pH-redesign of the best backbones (the real play).**
Take 2-3 backbones (prefer cC4_00 / cC_12; **avoid cC6** — unresolved PDB novelty
risk), run `run_ph_redesign` with a lambda sweep (e.g. 0.1, 0.3, 0.5, 0.7),
collect the Pareto front, take the top ~10 designs per backbone, fold each with
Boltz-2 on the DGX, and re-apply the original filters (RMSD, interface
confidence, contacts) + novelty checks. De novo constraint is preserved: the
starting backbones are our own de novo BoltzGen designs, not known binders
(disclose the method in the submission as required).

**The His433 angle.** EGFR domain-III His433 is the natural pH sensor: protonated
(+1) at pH 6.5, neutral at 7.4 (consistent with the PMC9703009 precedent, where
Asp/Glu on an antibody facing receptor His433 increased acidic-pH binding, and
with His433's role in the active 1IVO structure). Proton-PottsMPNN can exploit
it directly: placing a binder `ASP-D`/`GLU-D` against target `HIS-P` creates a
salt bridge at pH 6.5 that is lost at 7.4. The model's protonation-aware pairwise
couplings score exactly this microstate contrast, and `placement_by="scan_potts"`
should discover His433-facing positions automatically; if not, bias
`center_types=["ASP-P","GLU-P"]` with placement restricted to interface positions
near His433. Note the selectivity gap varies the *binder* centre's state while
the target's labelled state is fixed context — the His433 switch is captured
through the pairwise coupling term, which is the correct physical channel.

## Feasibility in ~40 h (Kaggle T4 / 16 GiB DGX slice)

**Feasible.** The critical path is CPU-bound and light:

1. Install on DGX (uv + Python 3.12 + foundry): ~1-2 h (heaviest step; watch
   atomworks build).
2. Featurize 2-3 candidate backbones + placement scans: <1 h CPU.
3. Lambda sweeps, ~20-40 designs per backbone, `n_jobs=4`: a few hours CPU.
4. Fold top ~10/backbone with Boltz-2 on DGX (~3 min each observed): ~1.5 h.
5. Filters + novelty + submission CSV: remainder.

Risks: (1) `install.sh` empty — manual uv install; (2) Python 3.12 required
(Kaggle default is older — `uv python install 3.12` fixes it); (3) no GPU needed
for design, so the 16 GiB slice and the 5.86 GiB RAM limit are non-issues here;
(4) the method is validated on PD-L1, not EGFR — treat outputs as hypotheses and
keep the Boltz-2/Chai-1 cross-checks. Do **not** spend time on RF3 (weights not
shipped) or retraining anything.

## Caveats

- Paper PDF could not be fetched (access failed); method description above is
  grounded in the repo README + source, not the manuscript text.
- `inference/design_ph.py` is 0 bytes in this snapshot; the notebook
  `design_ph.ipynb` would not parse as JSON here — use the README's Python API
  directly.
- Labeller AUPR: HIS 0.70, acids 0.33 (from README) — acid protonation labels are
  noisier; prefer HIS-P centres when the scan offers a choice.
