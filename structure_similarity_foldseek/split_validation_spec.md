# Split Validation Specification (Foldseek-Based OOD Splits)

## Objective

Validate generated OOD split TSV files by quantifying how similar each test sample is to the final training set, and by visualizing split structure.

This validation stage must produce:

- sample-level maximum similarity tables
- split-level summary statistics
- violin and ECDF plots across effective percentiles
- UMAP visualizations colored by `train` and `test`

---

## Inputs

### 1. Valid structure manifest

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/valid_structures.tsv`

Used to define the complete node list.

### 2. Pairwise Foldseek scores

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/pairwise_foldseek.tsv`

Use:

- `foldseek_score_max`

as the similarity value for validation.

### 3. OOD split TSV files

Directory:

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/`

For a given seed, split files are expected under:

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_<seed>/`

Examples include:

- `random_seed_0/split_p90.tsv`
- `random_seed_0/split_p95.tsv`
- `random_seed_0/split_p97.tsv`
- `random_seed_0/split_p99.tsv`

The split TSV format is the new train/test format produced by `generate_ood_splits.py`, including:

- `split`
- `component_id`
- `component_size`
- `effective_percentile`

---

## Primary Validation Quantity

For each test sample in a given split file:

1. take all samples assigned to `train`
2. compute similarity to each training sample
3. record the maximum similarity

This quantity is:

- `max_similarity_to_train`

Missing Foldseek pairs are treated as similarity `0` for this validation stage.

---

## Sample-Level Output

For each split file, write one TSV.

Example:

- `split_p95_test_vs_train.tsv`

Columns:

- `effective_percentile`
- `pdbid`
- `chain`
- `component_id`
- `component_size`
- `is_singleton`
- `max_similarity_to_train`
- `nearest_train_pdbid`
- `nearest_train_chain`

Rows:

- one row per test sample

---

## Summary Statistics

Write:

- `test_max_similarity_summary.tsv`
- `test_max_similarity_summary.txt`

These validation outputs must be written under:

- `.../ood_splits/random_seed_<seed>/validation/`

For each effective percentile, record:

- `percentile`
- `n_test`
- `min`
- `q1`
- `median`
- `mean`
- `q3`
- `p90`
- `p95`
- `max`
- `n_singletons_in_test`

---

## Plots

### 1. Violin Plot

Write:

- `test_max_similarity_violin.png`

Show, for each effective percentile:

- the distribution of `max_similarity_to_train`

### 2. ECDF Plot

Write:

- `test_max_similarity_ecdf.png`

Show, for each effective percentile:

- the empirical cumulative distribution of `max_similarity_to_train`

### 3. UMAP Visualization

Build a global 2D UMAP embedding from the Foldseek similarity graph using:

- observed Foldseek similarities from `pairwise_foldseek.tsv`
- missing pairs treated as similarity `0`
- distance defined from similarity for embedding

Write:

- `umap_embedding.tsv`

Columns:

- `pdbid`
- `chain`
- `umap1`
- `umap2`

Then for each split file write one plot:

- `umap_split_p90.png`
- `umap_split_p95.png`
- `umap_split_p97.png`
- `umap_split_p99.png`

Each UMAP plot must:

- color `train` and `test` differently
- show the full dataset

---

## Singleton Definition

A singleton is a connected component of size `1`.

For validation reporting:

- `is_singleton = 1` if `component_size == 1`
- otherwise `0`

---

## Execution

Script location:

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/validate_ood_splits.py`

Example:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/validate_ood_splits.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/pairwise_foldseek.tsv \
    --split_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits \
    --percentiles 90 95 97 99 \
    --random_seed 0 \
    --skip_umap
```

If `umap-learn` is installed, omit `--skip_umap` to also generate UMAP plots.

This command writes validation artifacts to:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/validation/`

### Seed Sweep Wrapper

Wrapper script:

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/run_ood_split_seed_sweep.py`

This wrapper runs:

- `generate_ood_splits.py`
- `validate_ood_splits.py`

for multiple seeds in one command, and stores outputs under:

- `.../ood_splits/random_seed_<seed>/`
- `.../ood_splits/random_seed_<seed>/validation/`

It also writes a seed-comparison summary and plot to the base output directory:

- `seed_sweep_test_fraction.tsv`
- `seed_sweep_test_fraction.png`

The plot shows:

- x-axis = `effective_percentile`
- y-axis = `final_test_fraction`
- one line per seed
