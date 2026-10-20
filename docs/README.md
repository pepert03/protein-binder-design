# Protein design website

A static page. No build, framework, CDN or external JavaScript dependencies.

## GitHub Pages

Copy index.html, assets/ and .nojekyll into the directory you want to publish. In the repository settings, choose Pages, Deploy from a branch, then the branch and directory containing index.html. If using /docs, place the whole site there. All runtime paths are relative.

For local preview, serve this folder with `python3 -m http.server 8000` and open http://localhost:8000. JavaScript modules need an HTTP server, so opening index.html directly as a file is not supported.

## What the views show

Challenge 1 (EGFR): a grey surface of human EGFR (the extracellular domain, 621 residues) and, one at a time, the Boltz 2 predicted binder for each of the 20 draft sequences in results/boltz2_csv20 of this repository. The page cycles to the next design every second, ordered by predicted ipTM from highest (r03_l07, 0.79) to lowest (r01_l07, 0.15). The caption names the design and its ipTM. Dragging pauses the cycle and the rotation; they resume 2.5 seconds after release. The Pause button stops both.

How the 20 poses share one surface: each prediction is its own complex with its own EGFR conformation. The grey surface is the one from r02_l07. For every other design, the EGFR residues within 15 Angstrom of the binder are superimposed onto the same residues of r02_l07 (rigid fit, local RMSD 0.25 to 4.0 Angstrom, see assets/designs.json), and the binder is drawn at the resulting position. A global superposition of the whole receptor would be worse, because EGFR domains move between predictions (global CA RMSD up to 20.7 Angstrom for rank07_native). A few binders sit within about 0.5 to 1 Angstrom of the r02_l07 atoms, so some contacts look tighter than in their own prediction.

Each binder is drawn as a ribbon cartoon built from the predicted backbone: hydrogen-bond secondary-structure assignment, a Catmull-Rom spline through the alpha carbons, and carbonyl-guided ribbon orientation, so helices appear as wide flat spiraling ribbons and loops as a thin rope. These are predicted coordinates, not evidence of experimental binding. Original research repository: https://github.com/pepert03/protein-binder-design.

Challenge 2 (TNF alpha): the grey molecular envelope of human TNF alpha, chains A, B and C from RCSB PDB 1TNF, fetched 20 October 2026. The original coordinates are included as assets/1TNF.pdb. The envelope is a Gaussian atomic density isosurface, not an exact solvent excluded surface. Source: https://www.rcsb.org/structure/1TNF

The recovered original AF2-Multimer crop complexes are shown: tnf_01 (the submitted 60-aa candidate) and 24 session-3 MPNN redesigns. The session-3 redesigns are research diagnostics, not the submitted ranks 2-20. All failed the screening gate; the best binder-to-entire-target normalized iPAE is 0.597.

For context, the tnf_01 crop receptor is rigidly fitted to the 1TNF trimer using sequence-matched receptor fragments (about 1.94 A CA RMSD). Session-3 crop receptors are fitted to the tnf_01 crop first (about 0.18-0.19 A RMSD), then use that same transform, avoiding trimer symmetry ambiguity. These are contextual superpositions of crop predictions, not predictions against the intact trimer. No binder sequence or internal coordinates have been redesigned for display.

The separate intact-trimer re-prediction of tnf_01 is not shown: it had a severe steric clash, minimum heavy-atom distance 0.27 A and 36 pairs below 2 A. The submitted ranks 2-20 from the 65-aa-family work are not shown; the intact audit found no real interface. No binding, affinity, pH response or mouse binding was experimentally demonstrated.

Original source archives recovered from the DGX on 20 October 2026: session2_evidence2.tar.gz (tnf_01 is session2_out/refilter/02/complex.pdb) and session3_evidence.tar.gz (session3_out/refilter/00-23/complex.pdb). The candidate sequence matches assets/challenge2-submission.csv exactly.

## Controls

Choose a challenge at top left. Drag with mouse or one finger to rotate. Scroll with a mouse or trackpad to zoom. Pause stops the slow rotation and the one second design cycle. Reset view restores the initial orientation. With the browser setting reduce motion, the page starts paused.

The supplied Challenge 2 submission CSV is preserved, unchanged, at assets/challenge2-submission.csv. Sequences are research candidates, not structural coordinates.
