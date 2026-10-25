#!/usr/bin/env bash
# setup_dgx_glp1r.sh -- bootstrap for the GLP1R campaign on the Comillas DGX.
# Run inside a Terminal-DGX session (GPU profile). Only ~ persists.
#
# Expected tree (transfer step places the two scripts via the browser session):
#   ~/glp1r/scripts/prepare_target.py
#   ~/glp1r/scripts/run_glp1r_colabdesign.py
# This script then creates ~/glp1r/{target,docs} and fills them.
set -euxo pipefail
mkdir -p ~/glp1r/target ~/glp1r/docs && cd ~/glp1r

# 1) reference structure + UniProt sequence (node has internet; verified 24 Oct).
#    prepare_target.py prefers the challenge-page segment FASTA if transferred;
#    otherwise it uses this full-length UniProt download (numbering base 1).
[ -f target/6X18.cif ] || curl -fsSL -o target/6X18.cif https://files.rcsb.org/download/6X18.cif
[ -f target/glp1r_P43220_24-432.fasta ] || [ -f target/P43220.fasta ] || \
  curl -fsSL -o target/P43220.fasta https://rest.uniprot.org/uniprotkb/P43220.fasta

# 2) extract receptor chain R + contacts. prepare_target.py resolves its outputs
#    relative to its own location: scripts/../target and scripts/../docs = the tree above.
python3 scripts/prepare_target.py

# 3) validate the crop (no GPU, no downloads) and STOP for human review of fragments
python3 scripts/run_glp1r_colabdesign.py --crop-only

# 4) python env: reuse the challenge2 venv if present (colabdesign + jax-cuda installed)
if [ -f ~/tnf/venv/bin/activate ]; then
  . ~/tnf/venv/bin/activate
else
  python3 -m venv venv && . venv/bin/activate
  pip install -q -U pip
  pip install -q "jax[cuda12]" git+https://github.com/sokrypton/ColabDesign.git biopython pandas
fi
python -c "import jax; print(jax.__version__, jax.devices())"

# 5) launch campaign in background; tmux keeps it alive while the session lives
nohup python -u scripts/run_glp1r_colabdesign.py > g3_run.log 2>&1 &
echo "campaign started; tail -f ~/glp1r/g3_run.log"
