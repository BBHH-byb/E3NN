# NMR Protein Dataset

This dataset contains structural and dynamical information derived from NMR protein ensembles.

Each row corresponds to a single `(pdbid, chain)` entry.

The dataset was generated using centroid structures and precomputed PCA results.

---

# Dataset Specification

The dataset contains **8155 PDB entries**.

For each entry, a representative structure is selected as the **centroid frame**, defined as the structure whose coordinates are closest to the average structure of all frames in the ensemble.

The order of proteins in the dataset follows the list defined in:

`final_list_RMSD.txt`

# Dataset Columns

| Column | Description |
|------|-------------|
| pdbid | PDB identifier |
| chain | Chain identifier (e.g. A, AB) |
| nresid | Total number of residues |
| resid_type | Residue types in sequence order (3-letter code) |
| resid_chain | Chain identity for each residue |
| coordinate | Flattened CA coordinates (x1,y1,z1,...,xN,yN,zN) |
| PC | Principal component vector (PC1 only) |
| var | Variance corresponding to PC1 |

---

# Structural Representation

Only **Cα atoms** are used to represent the protein structure.

Coordinates are stored as:

x1,y1,z1,x2,y2,z2,...,xN,yN,zN

where N is the number of residues.

---

# PCA Representation

The dataset includes only **PC1** from PCA analysis.

PC vectors are stored as a flattened vector of length:

3 × nresid

The variance column contains the eigenvalue associated with PC1.

---
