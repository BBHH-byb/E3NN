# Codex Agent Rules

## General rules

- Always follow SPEC.md when implementing dataset builders
- Do not invent file paths
- Always verify files exist before reading
- Never modify the dataset schema without updating SPEC.md

## Coding rules

- Language: Python
- Use pandas and numpy only unless necessary
- Prefer simple procedural scripts over complex class structures
- Always print warnings when files are missing

## Dataset rules

The dataset schema must always follow SPEC.md exactly.

Columns must be:

pdbid
chain
nresid
resid_type
resid_chain
coordinate
PC
var

Never add or remove columns.

## File locations

Base directory:

../../../domenico/27_NMR_screening/

Master list:

final_list_RMSD.txt

PCA files:

pca_results/{pdbid}_{chain}_PCs.txt

Variance files:

pca_results/{pdbid}_{chain}_pc_variances.txt

Centroid structures:

pdb_centroids/{pdbid}_{chain}_centroid.pdb
