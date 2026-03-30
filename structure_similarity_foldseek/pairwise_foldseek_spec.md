# Pairwise Foldseek + TM-align Specification

## Purpose

This specification defines the structure-similarity stage of the NMR dataset OOD pipeline using:

- Foldseek for pair discovery
- TM-align for pair rescoring

The goal of this stage is:

1. reconstruct CA-only pseudo-structures from the input CSV
2. run Foldseek all-vs-all retrieval on valid structures
3. collect all observed Foldseek hits
4. collapse those hits into unique unordered pairs
5. run TM-align on every observed unique pair
6. use TM-align scores as the primary downstream structural similarity values

This specification is the source of truth for implementing:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/run_pairwise_foldseek.py`

All implementation code for this stage must remain inside:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/`

---

## Input

### Input CSV

Path:

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/csv_input/input_monomer_dataset.csv`

Each row corresponds to a single `(pdbid, chain)` entry.

---

## Required Columns

The script must read the following columns:

- `pdbid`
- `chain`
- `nresid`
- `coordinate`

Other columns may exist but are not required for this stage.

---

## Structural Representation

The `coordinate` column stores flattened CA coordinates:

`x1,y1,z1,x2,y2,z2,...,xN,yN,zN`

The script must reshape this into:

`(nresid, 3)`

Validation rule:

- number of coordinate values must equal `3 * nresid`
- otherwise mark structure invalid

Only CA coordinates are available from the CSV, so each structure must be reconstructed as a CA-only pseudo-structure.

---

## Structure Reconstruction

Each valid row must be converted into a temporary PDB file containing CA atoms only.

### PDB generation rules

For each residue index `i`:

- atom name: `CA`
- residue name: `ALA`
- chain ID: use input `chain` if valid as a single character; otherwise use `A`
- residue number: sequential starting from 1
- coordinates: from reshaped `(nresid, 3)` array
- occupancy: `1.00`
- B-factor: `0.00`
- element: `C`

### Output pseudo-structure naming

Each generated PDB file must be named deterministically as:

`{pdbid}_{chain}.pdb`

All generated PDB files must be written into a dedicated structure directory.

---

## Invalid Structure Handling

A structure must be marked invalid and skipped if any of the following holds:

- coordinate parsing fails
- number of coordinate values does not equal `3 * nresid`
- `nresid <= 0`
- any coordinate is non-finite (`NaN`, `inf`)
- PDB generation fails

Invalid structures must be logged and excluded from downstream Foldseek and TM-align processing.

Processing must continue after invalid entries.

---

## Pair Discovery Strategy

This stage does not run dense all-vs-all TM-align directly.

Instead:

1. Foldseek is used as the scalable all-vs-all retrieval layer on valid structures.
2. Every observed Foldseek pair is kept.
3. TM-align is run on every unique unordered pair observed by Foldseek.

This means:

- pair discovery depends on Foldseek hit coverage
- pair scoring for downstream use depends on TM-align

Missing Foldseek pairs remain absent from the final pair table.

---

## Foldseek Database Preparation

The script must:

1. generate CA-only PDB files for all valid structures
2. create a Foldseek database from these PDB files
3. run Foldseek all-vs-all search on that database
4. parse the resulting directional Foldseek hits
5. collapse them into unique unordered observed pairs
6. run TM-align on those unique pairs

---

## Executable Requirements

A working Foldseek executable and a working TM-align executable must be installed before running this stage.

The script must accept:

- `--foldseek_path`
- `--tmalign_path`

If either executable cannot be found, the script must fail clearly before the expensive scoring stages begin.

---

## Foldseek Output Fields

The script must capture Foldseek output fields sufficient for auditing and downstream diagnostics.

Required parsed fields:

- query structure ID
- target structure ID
- chosen Foldseek primary score field
- alignment length if available
- e-value if available
- bits if available
- sequence identity if available
- structural similarity field if available

In the current implementation, the primary Foldseek field used for raw hit storage is:

- `alntmscore`

This Foldseek score is not the final downstream similarity score. It is retained for provenance and comparison.

---

## TM-align Rescoring

Every unique unordered pair observed by Foldseek must be rescored with TM-align.

TM-align is run once per unique pair:

- `TMalign pdb1 pdb2`

TM-align outputs two TM-scores, one normalized by each structure length.

For each pair, retain:

- `tm_score_max = max(tm_score_1, tm_score_2)`
- `tm_score_min = min(tm_score_1, tm_score_2)`

The primary structural similarity for downstream use must be:

- `tm_score_max`

This rule must be documented in the run log.

---

## Pair Construction

Foldseek search results are directional.

Convert them to unique unordered pairs:

A pair `(i, j)` is included if:

- at least one valid Foldseek hit exists between `i` and `j`

Final pair list must:

- contain unique unordered pairs
- exclude self-pairs
- remove duplicates

Each pair must appear exactly once.

Canonical pair ordering rule:

- sort first by `pdbid`
- then by `chain`

So pair `(i, j)` must always be stored with the lexicographically smaller `(pdbid, chain)` first.

If both Foldseek directions `(i -> j)` and `(j -> i)` exist and differ numerically, retain:

- `foldseek_score_max = max(score_ij, score_ji)`
- `foldseek_score_min = min(score_ij, score_ji)`

If multiple raw Foldseek hits exist for the same direction, keep the highest primary Foldseek score for that direction before computing:

- `foldseek_score_max`
- `foldseek_score_min`
- `hit_direction_count`

TM-align rescoring is then run on the collapsed unique pair table.

---

## Output

Default output directory:

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score`

### 1. Valid structure manifest

File:

`valid_structures.tsv`

Columns:

- `pdbid`
- `chain`
- `pdb_path`
- `nresid`

---

### 2. Invalid structures

File:

`invalid_structures.tsv`

Columns:

- `pdbid`
- `chain`
- `reason`

---

### 3. Raw Foldseek hits

File:

`foldseek_raw_hits.tsv`

This file stores parsed directional Foldseek results before undirected collapsing.

Required columns:

- `query_pdbid`
- `query_chain`
- `target_pdbid`
- `target_chain`
- `primary_score`
- additional Foldseek fields if available

Self-hits may be retained in this raw file for auditability, but must not appear in the final pair file.

---

### 4. Unique pairwise table

File:

`pairwise_foldseek.tsv`

Columns (exact order):

- `pdbid1`
- `chain1`
- `pdbid2`
- `chain2`
- `foldseek_score_max`
- `foldseek_score_min`
- `tm_score_max`
- `tm_score_min`
- `hit_direction_count`

Rules:

- each pair appears once
- unordered
- no self-pairs
- no duplicate pairs
- `hit_direction_count` is 1 or 2 depending on whether one or both Foldseek directions were observed
- `tm_score_max` is the primary downstream similarity score

---

### 5. Run log

File:

`run.log`

Must include:

- input CSV path
- total rows
- valid structures
- invalid structures
- structure output directory
- Foldseek executable path
- TM-align executable path
- Foldseek database path
- Foldseek commands used
- Foldseek primary score field
- pair aggregation rule
- total raw hits
- total unique pairs
- runtime

---

### 6. Score correlation summary

File:

`score_correlation_summary.txt`

This file must summarize the one-to-one comparison between:

- `foldseek_score_max`
- `tm_score_max`

computed over the final unique pair table.

At minimum it must include:

- number of pairs included in the comparison
- Pearson correlation between `foldseek_score_max` and `tm_score_max`

It may also include simple descriptive statistics such as:

- mean and median of `foldseek_score_max`
- mean and median of `tm_score_max`

---

## Output Directory Structure

Example:

- `results_foldseek_based_TM_score/structures/`
- `results_foldseek_based_TM_score/foldseek_db/`
- `results_foldseek_based_TM_score/valid_structures.tsv`
- `results_foldseek_based_TM_score/invalid_structures.tsv`
- `results_foldseek_based_TM_score/foldseek_raw_hits.tsv`
- `results_foldseek_based_TM_score/pairwise_foldseek.tsv`
- `results_foldseek_based_TM_score/score_correlation_summary.txt`
- `results_foldseek_based_TM_score/run.log`

The exact internal subdirectory names may vary, but the final output files above must exist.

---

## Execution

Script location:

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/run_pairwise_foldseek.py`

Example:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/run_pairwise_foldseek.py \
    --input_csv /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/csv_input/input_monomer_dataset.csv \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score \
    --foldseek_path /path/to/foldseek \
    --tmalign_path /path/to/TMalign
```

---

## Terminal Progress Reporting

The script must print progress updates to the terminal during execution.

Required visibility:

- input CSV path
- total CSV rows
- current processed rows during structure generation
- valid structure count
- invalid structure count
- estimated remaining time during structure generation
- Foldseek stage start and finish messages for:
  - `createdb`
  - `search`
  - `convertalis`
- elapsed time per Foldseek stage
- TM-align rescoring progress:
  - processed pair count
  - total pair count
  - ETA
- final counts for raw hits and unique pairs

If Foldseek emits progress lines to stdout or stderr, those lines should be streamed to the terminal when practical.
