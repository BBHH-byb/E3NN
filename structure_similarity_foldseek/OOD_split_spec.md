# OOD Split Specification (TM-score-Based Graph Tau Exploration + Split Generation)

## Objective

This stage now uses `tm_score_max` from `pairwise_foldseek.tsv` as the structural similarity score.

The work must be split into two scripts:

- `subagent1`: summarize graph structure across many `graph_tau` values
- `subagent2`: generate OOD splits across many `effective_tau` values and many random seeds after fixing the chosen `graph_tau`

This specification is the source of truth for implementing:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/summarize_graph_tau_components.py`
- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_ood_splits.py`

All implementation code for this stage must remain inside:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/`

---

## Input

### 1. Valid structure manifest

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv`

### 2. Pairwise score table

`/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv`

This file now contains both Foldseek-derived and TM-align-derived values.

For graph construction and OOD splitting, the score that must be used is:

- `tm_score_max`

Do not use `foldseek_score_max` as the primary similarity for this stage.

---

## Percentile Set

The percentile sweep for both graph tau exploration and effective tau splitting must be:

- `0`
- `20`
- `40`
- `60`
- `80`
- `90`
- `92`
- `94`
- `95`
- `96`
- `97`
- `98`
- `99`

All percentile thresholds are computed over observed `tm_score_max` values only.

Missing pairs must not be filled with zeros during percentile estimation.

Let:

- `S_obs` = observed `tm_score_max` values from `pairwise_foldseek.tsv`

Then for percentile `p`:

- `tau(p) = percentile(S_obs, p)`

---

## Subagent1: Graph Tau Exploration Summary

### Purpose

This script explores many candidate `graph_tau` values and summarizes the connected-component structure induced by each threshold.

It is intended to help choose a good `graph_tau` before running the actual OOD splitting stage.

### Graph definition

For a given `graph_tau`:

- node = protein
- edge between `i, j` if observed `tm_score_max(i, j) >= graph_tau`

Compute connected components on this graph.

### Required summary statistics

For each percentile in the full list above, compute:

- `graph_percentile`
- `graph_tau`
- `n_proteins`
- `observed_pairs`
- `retained_edges`
- `n_components`
- `n_singletons`
- `largest_component_size`

Optional but useful extra fields may also be included, such as:

- `mean_component_size`
- `median_component_size`
- `non_singleton_component_count`

### Output

This script must write:

- `graph_tau_component_summary.tsv`

Default output directory:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/graph_tau_exploration/`

### Purpose of this output

This table is used to choose a single final `graph_tau` for downstream OOD split generation.

---

## Subagent2: OOD Split Generation From a Fixed Graph Tau

### Purpose

After selecting one `graph_tau` from the Subagent1 summary, generate OOD splits by sweeping `effective_tau` values while keeping the graph definition fixed, and do this across multiple random seeds in one run.

### Fixed graph tau

The script must accept either:

- `--graph_percentile`
- `--graph_tau`

That fixed value defines the graph used for connected-component computation.

### Graph construction

Build `G(graph_tau)`:

- node = protein
- edge between `i, j` if observed `tm_score_max(i, j) >= graph_tau`

Compute connected components once from this graph.

### Candidate test selection and seed sweep

Use the same policy as the current OOD split logic.

For each seed:

1. the largest connected component is always assigned to the non-test side
2. randomly sample non-largest whole components
3. stop when candidate test size reaches approximately 30% of the full dataset

Use:

- `candidate_test_fraction = 0.30`

All proteins not selected into the candidate test pool belong to the non-test pool.

The script must accept multiple seeds, and for each seed it must independently:

- sample candidate test components
- run the full effective tau sweep
- write outputs to that seed's own directory

### Effective tau sweep

The script must accept either:

- `--effective_percentiles`
- `--taus`

If `--effective_percentiles` is used, then for each percentile in the list:

- compute `effective_tau = percentile(S_obs, p)`

where `S_obs` is the observed `tm_score_max` distribution.

If `--taus` is used, then the provided absolute tau values must be used directly without percentile conversion.

### Initial test acceptance

For each protein `i` in the candidate test pool, compute:

- `max_similarity_to_non_test(i)`

where the maximum is taken over all proteins in the non-test pool using `tm_score_max`.

In the initial pass:

- accept `i` into provisional test if `max_similarity_to_non_test(i) <= effective_tau`
- otherwise return `i` to train

### Iterative consistency filtering

After the initial pass:

1. define current training set as:
   - original non-test proteins
   - plus rejected candidate-test proteins
2. for each currently accepted test protein `i`, compute:
   - `max_similarity_to_train(i)`
3. if `max_similarity_to_train(i) <= effective_tau`, keep `i` in test
4. otherwise move `i` from test to train
5. repeat until no accepted test protein violates the criterion

At convergence:

- final `test` contains only proteins satisfying the final criterion with respect to final `train`
- final `train` contains all remaining proteins

### Output format

For each effective percentile, generate one TSV:

Schema:

`pdbid\tchain\tsplit\tcomponent_id\tcomponent_size\tgraph_tau\tcandidate_test\teffective_tau\teffective_percentile\tmax_similarity_to_non_test`

Where:

- `split ∈ {train, test}`
- `component_id` and `component_size` come from the fixed graph defined by the chosen `graph_tau`
- `graph_tau` is the fixed graph-building threshold
- `effective_tau` changes across the effective percentile sweep

File naming:

- `split_p0.tsv`
- `split_p20.tsv`
- `split_p40.tsv`
- `split_p60.tsv`
- `split_p80.tsv`
- `split_p90.tsv`
- `split_p92.tsv`
- `split_p94.tsv`
- `split_p95.tsv`
- `split_p96.tsv`
- `split_p97.tsv`
- `split_p98.tsv`
- `split_p99.tsv`

### Split summary

The script must also write:

- `split_summary.tsv`

Required fields:

- `graph_percentile`
- `graph_tau`
- `effective_percentile`
- `effective_tau`
- `n_proteins`
- `observed_pairs`
- `retained_graph_edges`
- `n_components`
- `n_singletons`
- `largest_component_size`
- `candidate_test_size`
- `candidate_test_fraction`
- `provisional_test_size`
- `initial_rejected_candidate_size`
- `iterative_removals`
- `final_train_size`
- `final_test_size`
- `final_train_fraction`
- `final_test_fraction`

All outputs for a given seed must be written under:

- `.../ood_splits/random_seed_<seed>/`

In addition, Subagent2 must write aggregate cross-seed summaries to the base output directory:

- `all_seed_split_summary.tsv`
- `all_seed_split_summary.txt`

The aggregate TSV must contain one row per `(seed, effective_percentile)` pair and include:

- `seed`
- all fields from `split_summary.tsv`

The aggregate TXT must summarize the same table in a human-readable way so that seed-by-seed differences across all effective tau levels can be inspected at once.

---

## Execution

### Subagent1 script

Script location:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/summarize_graph_tau_components.py`

This script must accept at minimum:

- `--valid_structures`
- `--pairwise_scores`
- `--score_column tm_score_max`
- `--graph_percentiles 0 20 40 60 80 90 92 94 95 96 97 98 99` or `--graph_taus`
- `--output_dir`

### Subagent2 script

Script location:

- `/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_ood_splits.py`

This script must accept at minimum:

- `--valid_structures`
- `--pairwise_scores`
- `--score_column tm_score_max`
- `--graph_percentile` or `--graph_tau`
- `--candidate_test_fraction 0.30`
- `--effective_percentiles 0 20 40 60 80 90 92 94 95 96 97 98 99` or `--taus`
- `--target_test_fraction 0.10`
- `--seeds`

Example:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_ood_splits.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv \
    --score_column tm_score_max \
    --graph_percentile 95 \
    --candidate_test_fraction 0.30 \
    --effective_percentiles 0 20 40 60 80 90 92 94 95 96 97 98 99 \
    --target_test_fraction 0.10 \
    --seeds 0 1 2 3 \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits
```

Alternative:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_ood_splits.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv \
    --score_column tm_score_max \
    --graph_tau 0.657210 \
    --candidate_test_fraction 0.30 \
    --taus 0.50 0.60 0.70 \
    --target_test_fraction 0.10 \
    --seeds 0 1 2 3 \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits
```

---

## Final Summary

This stage no longer defines OOD splits directly from Foldseek scores. Instead, it uses `tm_score_max` from `pairwise_foldseek.tsv` as the similarity score. The workflow is intentionally split into two code paths:

- `subagent1` explores many candidate `graph_tau` values and summarizes component structure
- `subagent2` fixes the selected `graph_tau` and generates OOD splits across many `effective_tau` values and many seeds using the existing candidate-test and iterative filtering policy

This separates graph-threshold selection from final split generation and makes the chosen graph tau easier to justify from summary statistics.
