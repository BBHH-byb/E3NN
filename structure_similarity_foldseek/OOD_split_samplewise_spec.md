# OOD Split Samplewise Specification (TM-score-Based, No Graph Components)

## Objective

Define a samplewise OOD split procedure that:

- uses `tm_score_max` from `pairwise_foldseek.tsv`
- does not use graph construction or connected components
- selects test samples directly at the sample level
- targets a final test fraction of approximately 10%
- then draws a validation set equal to 10% of the full dataset size from the remaining non-test pool
- repeats the procedure across multiple random seeds

This specification is the source of truth for implementing:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_samplewise_ood_splits.py`

---

## Input

### 1. Valid structure manifest

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv`

### 2. Pairwise score table

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv`

For this stage, the primary similarity score must be:

- `tm_score_max`

Do not use graph tau or connected-component structure.

---

## Effective Tau

The script must support two ways to specify effective thresholds:

- `--effective_percentiles`
- `--taus`

If `--effective_percentiles` is used:

- compute `effective_tau = percentile(S_obs, p)`

where:

- `S_obs` = observed `tm_score_max` values from `pairwise_foldseek.tsv`

If `--taus` is used:

- the provided absolute tau values are used directly

---

## Split Policy

### High-level idea

For each seed and each effective tau:

1. shuffle all samples with the given random seed
2. start with every sample in train
3. scan the shuffled order
4. move a sample from train to test only if it satisfies the effective tau criterion with respect to the current train set
5. stop when the target test size is reached, or when no more eligible samples remain
6. from the remaining non-test pool, draw a validation set uniformly at random
7. assign all remaining samples to train

This creates a samplewise OOD split without any component constraint, with a random validation split taken after the OOD test split has been finalized.

### Effective tau criterion

For candidate sample `i`, define:

- `max_similarity_to_train(i)`

where the maximum is taken over all other samples currently still assigned to train using `tm_score_max`.

Sample `i` can be moved into test if:

- `max_similarity_to_train(i) <= effective_tau`

Because the train set only shrinks as samples are moved to test, once a sample satisfies the criterion and is accepted, it remains valid.

### Target size

Use:

- `target_test_fraction = 0.10`
- `target_validation_fraction = 0.10`

Let:

- `target_test_size = round(target_test_fraction * n_proteins)`
- `target_validation_size = round(target_validation_fraction * n_proteins)`

The script should keep accepting eligible samples until:

- `final_test_size == target_test_size`

or until the shuffled list is exhausted.

Therefore:

- exact 10% is achieved when enough eligible samples exist
- otherwise the final test set may be smaller than target
- validation is then drawn from the remaining pool with target size equal to 10% of the full dataset size
- if the remaining non-test pool is smaller than the requested validation size, validation uses the entire remaining pool

---

## Seed Sweep

The script must accept:

- `--seeds`

For each seed:

- generate independent sample orders
- run the full effective tau sweep
- write outputs to:

`.../ood_splits_samplewise/random_seed_<seed>/`

---

## Output Format

For each seed and each effective tau, write one TSV:

- `split_p90.tsv`
- `split_p95.tsv`
- or, if absolute taus are used:
- `split_tau_0p500000.tsv`

Columns:

- `pdbid`
- `chain`
- `split`
- `effective_tau`
- `effective_percentile`
- `max_similarity_to_train`

Where:

- `split ∈ {train, validation, test}`
- `max_similarity_to_train` is populated for accepted test samples

---

## Per-Seed Summary

For each seed, write:

- `split_summary.tsv`

Required fields:

- `seed`
- `effective_percentile`
- `effective_tau`
- `n_proteins`
- `target_test_size`
- `target_validation_size`
- `final_train_size`
- `final_validation_size`
- `final_test_size`
- `final_train_fraction`
- `final_validation_fraction`
- `final_test_fraction`
- `median_max_similarity_to_train`
- `target_test_fraction`
- `target_validation_fraction`

`median_max_similarity_to_train` must be computed only over accepted test samples for that `(seed, effective_tau)` run, using the `max_similarity_to_train` values available from the observed pairwise `tm_score_max` table.

---

## Aggregate Cross-Seed Summary

The script must also write to the base output directory:

- `all_seed_split_summary.tsv`
- `all_seed_split_summary.txt`

The aggregate TSV must contain one row per `(seed, effective_tau)` run.

It must include `median_max_similarity_to_train` so seeds can be compared not only by final test size/fraction but also by the central tendency of available nearest-train similarity among accepted test samples.

The aggregate TXT must summarize seed-by-seed differences across all effective tau settings in a human-readable format, including:

- `final_validation_fraction`
- `final_validation_size`
- `final_test_fraction`
- `final_test_size`
- `median_max_similarity_to_train`

---

## Execution

Script location:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_samplewise_ood_splits.py`

Example using percentiles:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_samplewise_ood_splits.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv \
    --score_column tm_score_max \
    --effective_percentiles 90 95 97 99 \
    --target_test_fraction 0.10 \
    --target_validation_fraction 0.10 \
    --seeds 0 1 2 3 \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits_samplewise
```

Example using absolute taus:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_samplewise_ood_splits.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv \
    --score_column tm_score_max \
    --taus 0.4 0.5 0.6 0.7 \
    --target_test_fraction 0.10 \
    --target_validation_fraction 0.10 \
    --seeds 0 1 2 3 \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits_samplewise
```

---

## Final Summary

This procedure removes graph and component constraints entirely. Instead, it directly selects OOD test samples at the sample level by enforcing `max_similarity_to_train <= effective_tau` while greedily growing the test set toward 10% for each seed. After the test split is fixed, it draws a random validation split equal to 10% of the full dataset size from the remaining pool. This makes the split definition simpler and lets the user study how absolute or percentile-based effective tau values affect the final train/validation/test proportions across seeds.
