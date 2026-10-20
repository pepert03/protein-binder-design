# Protein design website

A static page. No build, framework, CDN or external JavaScript dependencies.

## GitHub Pages

Copy index.html, assets/ and .nojekyll into the directory you want to publish. In the repository settings, choose Pages, Deploy from a branch, then the branch and directory containing index.html. If using /docs, place the whole site there. All runtime paths are relative.

For local preview, serve this folder with `python3 -m http.server 8000` and open http://localhost:8000. JavaScript modules need an HTTP server, so opening index.html directly as a file is not supported.

## What the views show

Challenge 1: a grey surface of human EGFR and the blue ribbon cartoon of the candidate r02_l07 (helices as flattened coiled ribbons, strands as flat arrows, loops as thin tubes). The cartoon follows the exact CA trace of the Boltz 2 prediction (r02_l07.cif, 71 residues). Because only CA positions were carried over, the secondary structure shown is derived from CA geometry (CA i to i+4 distances and local curvature, DSSP-like heuristic), not from full-atom DSSP or any experimental assignment. These are predicted coordinates, not evidence of experimental binding. Original research repository: https://github.com/pepert03/protein-binder-design (may require access).

Challenge 2: the grey molecular envelope of human TNF alpha, chains A, B and C from RCSB PDB 1TNF, fetched 20 October 2026. The original coordinates are included as assets/1TNF.pdb. The envelope is a Gaussian atomic density isosurface, not an exact solvent excluded surface. Source: https://www.rcsb.org/structure/1TNF

No blue candidate is shown for Challenge 2 because the earlier candidate complex coordinates were not recoverable. The sequence alone does not define its pose. Do not use the EGFR binder mesh or an invented helix as a replacement.

The session 2 TNF candidate sequence is:
EKEPEFLLDYLVEEARKRKESDPEAAEELKRDLYYLLKLLKEEDDEEKVEKEVEERVEGK

Its earlier intact trimer prediction was weak (iPTM about 0.31) with steric clashes. No binding, affinity, pH response or mouse binding was experimentally demonstrated.

## Replacing the TNF view later

Recover session2_evidence2.tar.gz, containing session2_out/refilter/02/complex.pdb and session2_out/intact/00/intact.pdb. Use the complex coordinates without altering or inventing the binder pose. Generate grey target and blue binder meshes, replace the tnf entry in assets/structures.js, and change the target only caption in assets/app.js. Keep the in silico and geometry caveats.

## Controls

Choose a challenge at top left. Drag with mouse or one finger to rotate. Scroll with a mouse or trackpad to zoom. Reset view restores the initial orientation. Motion is not automatic.

The supplied Challenge 2 submission CSV is preserved, unchanged, at assets/challenge2-submission.csv. Sequences are research candidates, not structural coordinates.
