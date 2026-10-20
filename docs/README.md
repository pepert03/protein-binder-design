# Protein design website

A static page. No build, framework, CDN or external JavaScript dependencies.

## GitHub Pages

Copy index.html, assets/ and .nojekyll into the directory you want to publish. In the repository settings, choose Pages, Deploy from a branch, then the branch and directory containing index.html. If using /docs, place the whole site there. All runtime paths are relative.

For local preview, serve this folder with `python3 -m http.server 8000` and open http://localhost:8000. JavaScript modules need an HTTP server, so opening index.html directly as a file is not supported.

## What the views show

Challenge 1 (EGFR): a grey surface of human EGFR (the extracellular domain, 621 residues) and, one at a time, the Boltz 2 predicted binder for each of the 20 draft sequences in results/boltz2_csv20 of this repository. The page cycles to the next design every second, ordered by predicted ipTM from highest (r03_l07, 0.79) to lowest (r01_l07, 0.15). The caption names the design and its ipTM. Dragging pauses the cycle and the rotation; they resume 2.5 seconds after release. The Pause button stops both.

How the 20 poses share one surface: each prediction is its own complex with its own EGFR conformation. The grey surface is the one from r02_l07. For every other design, the EGFR residues within 15 Angstrom of the binder are superimposed onto the same residues of r02_l07 (rigid fit, local RMSD 0.25 to 4.0 Angstrom, see assets/designs.json), and the binder is drawn at the resulting position. A global superposition of the whole receptor would be worse, because EGFR domains move between predictions (global CA RMSD up to 20.7 Angstrom for rank07_native). A few binders sit within about 0.5 to 1 Angstrom of the r02_l07 atoms, so some contacts look tighter than in their own prediction.

The ribbon cartoon follows the exact CA trace of each Boltz 2 structure. Only CA positions are used, so secondary structure is derived from CA geometry (virtual dihedral and CA i to i+3 distances, a DSSP-like heuristic), not from full-atom DSSP or experiment. Helices are flattened coiled ribbons, strands are flat arrows, loops are thin tubes. These are predicted coordinates, not evidence of experimental binding. Original research repository: https://github.com/pepert03/protein-binder-design.

Challenge 2 (TNF alpha): the grey molecular envelope of human TNF alpha, chains A, B and C from RCSB PDB 1TNF, fetched 20 October 2026. The original coordinates are included as assets/1TNF.pdb. The envelope is a Gaussian atomic density isosurface, not an exact solvent excluded surface. Source: https://www.rcsb.org/structure/1TNF

No candidate is shown for Challenge 2. The repository was searched in full (challenge2/ and results/): it has no coordinates for the submitted session 2 sequence or for the other submitted sequences. The Challenge 2 PDB files it holds are inputs, controls and diagnostic runs of other sequences (an unmodified control, MPNN repair variants, 50 aa crops), and their own audit found no real interface. The sequence alone does not define a pose, so nothing is drawn.

The session 2 TNF candidate sequence is:
EKEPEFLLDYLVEEARKRKESDPEAAEELKRDLYYLLKLLKEEDDEEKVEKEVEERVEGK

Its earlier intact trimer prediction was weak (iPTM about 0.31) with steric clashes. No binding, affinity, pH response or mouse binding was experimentally demonstrated.

## Adding the Challenge 2 candidate later

Recover session2_evidence2.tar.gz, containing session2_out/refilter/02/complex.pdb and session2_out/intact/00/intact.pdb. Use the complex coordinates without altering or inventing the binder pose, then add a designs entry for the tnf key in assets/structures.js and update the caption text in assets/app.js. Keep the in silico and geometry caveats.

## Controls

Choose a challenge at top left. Drag with mouse or one finger to rotate. Scroll with a mouse or trackpad to zoom. Pause stops the slow rotation and the one second design cycle. Reset view restores the initial orientation. With the browser setting reduce motion, the page starts paused.

The supplied Challenge 2 submission CSV is preserved, unchanged, at assets/challenge2-submission.csv. Sequences are research candidates, not structural coordinates.
