# NMR Protein Dataset Builder (PC1-based)

## 1. Goal

Build a Python pipeline that aggregates structural and PCA-derived information from an NMR screening dataset into a structured dataset for machine learning.

The dataset only includes **monomer structures**

The pipeline **does not compute PCA**.  
All PCA results already exist and will simply be collected and merged.

The output dataset will contain **one row per (pdbid, chain)** entry.

The final dataset will be saved as:

- CSV file

---

# 2. Base Data Directory

All data are stored under the directory:

../../../domenico/27_NMR_screening/

All input paths described below are relative to this directory.

---

# 3. Master List File

Dataset entries are defined by the file:

final_list_RMSD.txt

Location:

../../../domenico/27_NMR_screening/final_list_RMSD.txt

This file contains **no header**.

Each line contains:

pdbid  chain  num_of_frames  length  max_RMSD

Example:

1abc A 20 156 4.83

From this file the pipeline extracts:

pdbid  
chain  
length  

Where:

pdbid → first dataset column  
chain → second dataset column  
length → number of residues (nresid)

---

# 4. Output Dataset Schema

Each row in the dataset will contain the following columns:

1. pdbid  
2. chain  
3. nresid  
4. resid_type  
5. resid_chain  
6. coordinate  
7. PC  
8. var  

Each row corresponds to one `(pdbid, chain)` entry.

---

# 5. Column Definitions

## pdbid

PDB identifier read from `final_list_RMSD.txt`.

Example:

1abc

---

## chain

Chain identifier read from `final_list_RMSD.txt`.

Examples:

A  
AB  

Here, the single chain cases are only selected.
---

## nresid

Total number of residues taken from the `length` field in `final_list_RMSD.txt`.

---

## resid_type

Residue types extracted from the centroid PDB file using **only CA atoms**.

Source file:

pdb_centroids/{pdbid}_{chain}_centroid.pdb

Residue names should be stored as a comma-separated list in residue order.

Example:

ALA,GLY,LEU,ASP,THR,...

---

## resid_chain

Chain identity for each residue extracted from the centroid PDB file.

Source file:

pdb_centroids/{pdbid}_{chain}_centroid.pdb

Only CA atoms are used.

Chain IDs should be recorded in residue order as a continuous string.

Example:

AAAAAAABBBBBBBB

---

## coordinate

Flattened CA coordinates from the centroid PDB file.

Source file:

pdb_centroids/{pdbid}_{chain}_centroid.pdb

Coordinates are stored as:

x1,y1,z1,x2,y2,z2,...,xN,yN,zN

Only CA coordinates are used.

Heavy atoms are not used in order to maintain consistent dimensionality across proteins.

---

## PC

Principal component vector corresponding to **PC1 only**.

Source file:

pca_results/{pdbid}_{chain}_PCs.txt

The pipeline must:

1. read the file as a matrix  
2. extract the **first column only**  

The vector should be stored as a comma-separated list.

Example:

0.0031,-0.0054,0.0012,...

Expected vector length:

3 × nresid

---

## var

Variance corresponding to **PC1**.

Source file:

pca_results/{pdbid}_{chain}_pc_variances.txt

The pipeline must read **only the first value** from this file.

Example:

0.35291

---

# 6. Structural Input Files

Centroid structures are stored in:

pdb_centroids/

File naming convention:

pdb_centroids/{pdbid}_{chain}_centroid.pdb

Example:

pdb_centroids/1abc_AB_centroid.pdb

From this file the pipeline extracts:

- residue types
- residue chain IDs
- CA coordinates

Only CA atoms are used.

---

# 7. PCA Input Files

Principal component files are stored in:

pca_results/

File naming convention:

pca_results/{pdbid}_{chain}_PCs.txt

Example:

pca_results/1abc_AB_PCs.txt

The pipeline reads:

- first column only (PC1)

---

# 8. Variance Files

Variance files are stored in:

pca_results/

File naming convention:

pca_results/{pdbid}_{chain}_pc_variances.txt

Example:

pca_results/1abc_AB_pc_variances.txt

The pipeline reads:

- first value only (PC1 variance)

---

# 9. Processing Workflow

Step 1  
Read `final_list_RMSD.txt` line by line.

Step 2  
Extract:

pdbid  
chain  
length  

Step 3  
Open centroid structure:

pdb_centroids/{pdbid}_{chain}_centroid.pdb

Parse CA atoms and collect:

residue types  
chain IDs  
coordinates  

Step 4  
Open PCA file:

pca_results/{pdbid}_{chain}_PCs.txt

Extract first column (PC1).

Step 5  
Open variance file:

pca_results/{pdbid}_{chain}_pc_variances.txt

Extract first value (PC1 variance).

Step 6  
Combine all collected data into one row.

Step 7  
Repeat for all entries.

---

# 10. Missing File Handling

If any required file is missing:

- print warning
- skip the entry

Example:

Warning: missing file for 1abc_AB

---

# 11. Output Files

The pipeline should generate:

csv_input/input_dataset.csv  
csv_input/input_dataset.pkl  

Optional:

csv_input/input_dataset.xlsx

---

# 12. Implementation Requirements

Language:

Python

Recommended libraries:

pandas  
numpy  
pathlib  
os  

Optional:

biopython

---

# 13. Suggested Project Structure

project/

SPEC.md

scripts/
    build_input_dataset.py

src/
    dataset_builder.py
    parse_centroid_pdb.py
    load_pc.py
    load_variance.py

csv_input/

---

# 14. Dataset Summary

Each dataset row contains:

pdbid  
chain  
number of residues  
residue sequence  
residue chain mapping  
CA coordinates  
PC1 vector  
PC1 variance  

The resulting dataset represents a **residue-level structural + dynamical dataset** suitable for machine learning models.
