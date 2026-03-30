#!/usr/bin/env python3
"""
Example run:

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
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


sys.dont_write_bytecode = True


DEFAULT_VALID_STRUCTURES = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv"
)
DEFAULT_PAIRWISE_SCORES = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits_samplewise"
)
SPLIT_COLUMNS = [
    "pdbid",
    "chain",
    "split",
    "effective_tau",
    "effective_percentile",
    "max_similarity_to_train",
]
SUMMARY_COLUMNS = [
    "seed",
    "effective_percentile",
    "effective_tau",
    "n_proteins",
    "target_test_size",
    "target_validation_size",
    "final_train_size",
    "final_validation_size",
    "final_test_size",
    "final_train_fraction",
    "final_validation_fraction",
    "final_test_fraction",
    "median_max_similarity_to_train",
    "target_test_fraction",
    "target_validation_fraction",
]
AGGREGATE_SUMMARY_COLUMNS = SUMMARY_COLUMNS


@dataclass(frozen=True)
class Node:
    pdbid: str
    chain: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate samplewise TM-score-based OOD splits with seed sweeps and no graph/component constraints."
    )
    parser.add_argument("--valid_structures", type=Path, default=DEFAULT_VALID_STRUCTURES)
    parser.add_argument("--pairwise_scores", type=Path, default=DEFAULT_PAIRWISE_SCORES)
    parser.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--score_column", default="tm_score_max")
    threshold_group = parser.add_mutually_exclusive_group()
    threshold_group.add_argument("--effective_percentiles", type=float, nargs="+")
    threshold_group.add_argument("--taus", type=float, nargs="+")
    parser.add_argument("--target_test_fraction", type=float, default=0.10)
    parser.add_argument("--target_validation_fraction", type=float, default=0.10)
    parser.add_argument("--seeds", type=int, nargs="+", default=[0])
    return parser.parse_args()


def format_float(value: float) -> str:
    return format(float(value), ".6f")


def percentile_suffix(percentile: float) -> str:
    return format(percentile, "g").replace(".0", "")


def tau_suffix(value: float) -> str:
    return format_float(value).replace(".", "p")


def resolve_seed_output_dir(base_output_dir: Path, random_seed: int) -> Path:
    seed_name = f"random_seed_{random_seed}"
    if base_output_dir.name == seed_name:
        return base_output_dir
    return base_output_dir / seed_name


def load_valid_nodes(path: Path) -> list[Node]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pdbid", "chain"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing valid_structures columns: {', '.join(sorted(missing))}")

        nodes: list[Node] = []
        seen: set[tuple[str, str]] = set()
        for row in reader:
            pdbid = (row.get("pdbid") or "").strip()
            chain = (row.get("chain") or "").strip()
            if not pdbid or not chain:
                continue
            key = (pdbid, chain)
            if key in seen:
                raise ValueError(f"duplicate valid structure entry: {pdbid}:{chain}")
            seen.add(key)
            nodes.append(Node(pdbid=pdbid, chain=chain))
    if not nodes:
        raise ValueError("valid_structures.tsv contains no valid rows")
    return nodes


def load_similarity_adjacency(
    path: Path,
    node_to_index: dict[tuple[str, str], int],
    score_column: str,
) -> tuple[list[list[tuple[int, float]]], np.ndarray]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pdbid1", "chain1", "pdbid2", "chain2", score_column}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing pairwise score columns: {', '.join(sorted(missing))}")

        adjacency: list[list[tuple[int, float]]] = [[] for _ in range(len(node_to_index))]
        scores: list[float] = []
        seen_pairs: set[tuple[int, int]] = set()
        for row in reader:
            left_key = ((row.get("pdbid1") or "").strip(), (row.get("chain1") or "").strip())
            right_key = ((row.get("pdbid2") or "").strip(), (row.get("chain2") or "").strip())
            if left_key not in node_to_index or right_key not in node_to_index:
                continue
            left_index = node_to_index[left_key]
            right_index = node_to_index[right_key]
            if left_index == right_index:
                continue
            pair_key = (min(left_index, right_index), max(left_index, right_index))
            if pair_key in seen_pairs:
                raise ValueError(
                    f"duplicate unordered pair in pairwise scores: {left_key[0]}:{left_key[1]} / {right_key[0]}:{right_key[1]}"
                )
            seen_pairs.add(pair_key)
            score = float((row.get(score_column) or "").strip())
            adjacency[pair_key[0]].append((pair_key[1], score))
            adjacency[pair_key[1]].append((pair_key[0], score))
            scores.append(score)
    return adjacency, np.asarray(scores, dtype=float)


def compute_percentile(scores: np.ndarray, percentile: float) -> float:
    if scores.size == 0:
        return 0.0
    return float(np.percentile(scores, percentile))


def resolve_effective_thresholds(args: argparse.Namespace, observed_scores: np.ndarray) -> list[tuple[str, float, str]]:
    if args.taus is not None:
        return [("", float(value), f"tau_{tau_suffix(float(value))}") for value in args.taus]
    percentiles = args.effective_percentiles
    if percentiles is None:
        percentiles = [90, 95, 97, 99]
    return [
        (percentile_suffix(percentile), compute_percentile(observed_scores, percentile), f"p{percentile_suffix(percentile)}")
        for percentile in percentiles
    ]


def max_similarity_to_train(node_index: int, train_mask: np.ndarray, adjacency: list[list[tuple[int, float]]]) -> float:
    best = 0.0
    for neighbor_index, score in adjacency[node_index]:
        if train_mask[neighbor_index] and score > best:
            best = score
    return best


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_aggregate_summary_text(path: Path, aggregate_rows: list[dict[str, str]]) -> None:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in aggregate_rows:
        key = row["effective_percentile"] if row["effective_percentile"] else f"tau={row['effective_tau']}"
        grouped.setdefault(key, []).append(row)

    def sort_key(value: str) -> float:
        if value.startswith("tau="):
            return float(value.split("=", 1)[1])
        return float(value)

    lines = ["Seed-by-Seed Samplewise Effective Tau Summary", ""]
    for effective_key in sorted(grouped, key=sort_key):
        rows = sorted(grouped[effective_key], key=lambda item: int(item["seed"]))
        fractions = np.asarray([float(row["final_test_fraction"]) for row in rows], dtype=float)
        medians = np.asarray(
            [
                float(row["median_max_similarity_to_train"])
                for row in rows
                if row["median_max_similarity_to_train"]
            ],
            dtype=float,
        )
        if effective_key.startswith("tau="):
            lines.append(f"effective_tau: {effective_key.split('=', 1)[1]}")
        else:
            lines.append(f"effective_percentile: {effective_key}")
        lines.append(f"seed_count: {len(rows)}")
        validation_fractions = np.asarray([float(row["final_validation_fraction"]) for row in rows], dtype=float)
        lines.append(f"final_test_fraction_min: {format_float(np.min(fractions))}")
        lines.append(f"final_test_fraction_median: {format_float(np.median(fractions))}")
        lines.append(f"final_test_fraction_max: {format_float(np.max(fractions))}")
        lines.append(f"final_validation_fraction_min: {format_float(np.min(validation_fractions))}")
        lines.append(f"final_validation_fraction_median: {format_float(np.median(validation_fractions))}")
        lines.append(f"final_validation_fraction_max: {format_float(np.max(validation_fractions))}")
        if medians.size > 0:
            lines.append(f"median_max_similarity_to_train_min: {format_float(np.min(medians))}")
            lines.append(f"median_max_similarity_to_train_median: {format_float(np.median(medians))}")
            lines.append(f"median_max_similarity_to_train_max: {format_float(np.max(medians))}")
        else:
            lines.append("median_max_similarity_to_train_min: NA")
            lines.append("median_max_similarity_to_train_median: NA")
            lines.append("median_max_similarity_to_train_max: NA")
        for row in rows:
            lines.append(
                "  "
                f"seed={row['seed']} "
                f"effective_tau={row['effective_tau']} "
                f"final_validation_fraction={row['final_validation_fraction']} "
                f"final_validation_size={row['final_validation_size']} "
                f"final_test_fraction={row['final_test_fraction']} "
                f"final_test_size={row['final_test_size']} "
                f"median_max_similarity_to_train={row['median_max_similarity_to_train'] or 'NA'}"
            )
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def run_for_seed(
    *,
    seed: int,
    nodes: list[Node],
    adjacency: list[list[tuple[int, float]]],
    threshold_specs: list[tuple[str, float, str]],
    target_test_fraction: float,
    target_validation_fraction: float,
    base_output_dir: Path,
) -> list[dict[str, str]]:
    output_dir = resolve_seed_output_dir(base_output_dir, seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    order = rng.permutation(len(nodes)).tolist()
    target_test_size = int(round(target_test_fraction * len(nodes)))
    target_validation_size = int(round(target_validation_fraction * len(nodes)))
    summary_rows: list[dict[str, str]] = []

    for effective_percentile, effective_tau, file_suffix in threshold_specs:
        train_mask = np.ones(len(nodes), dtype=bool)
        selected_scores: dict[int, float] = {}
        current_test_size = 0

        for node_index in order:
            if current_test_size >= target_test_size:
                break
            candidate_score = max_similarity_to_train(node_index, train_mask, adjacency)
            if candidate_score <= effective_tau:
                train_mask[node_index] = False
                selected_scores[node_index] = candidate_score
                current_test_size += 1

        split_rows: list[dict[str, str]] = []
        final_train_size = 0
        final_validation_size = 0
        final_test_size = 0
        remaining_train_indices = [index for index in range(len(nodes)) if train_mask[index]]
        validation_size = min(target_validation_size, len(remaining_train_indices))
        validation_indices = set(rng.permutation(remaining_train_indices).tolist()[:validation_size])
        for node_index, node in enumerate(nodes):
            if not train_mask[node_index]:
                split = "test"
                final_test_size += 1
                max_score = format_float(selected_scores[node_index])
            elif node_index in validation_indices:
                split = "validation"
                final_validation_size += 1
                max_score = ""
            else:
                split = "train"
                final_train_size += 1
                max_score = ""
            split_rows.append(
                {
                    "pdbid": node.pdbid,
                    "chain": node.chain,
                    "split": split,
                    "effective_tau": format_float(effective_tau),
                    "effective_percentile": effective_percentile,
                    "max_similarity_to_train": max_score,
                }
            )

        write_tsv(output_dir / f"split_{file_suffix}.tsv", SPLIT_COLUMNS, split_rows)
        selected_score_values = np.asarray(list(selected_scores.values()), dtype=float)
        median_max_similarity_to_train = (
            format_float(np.median(selected_score_values)) if selected_score_values.size > 0 else ""
        )

        summary_rows.append(
            {
                "seed": str(seed),
                "effective_percentile": effective_percentile,
                "effective_tau": format_float(effective_tau),
                "n_proteins": str(len(nodes)),
                "target_test_size": str(target_test_size),
                "target_validation_size": str(target_validation_size),
                "final_train_size": str(final_train_size),
                "final_validation_size": str(final_validation_size),
                "final_test_size": str(final_test_size),
                "final_train_fraction": format_float(final_train_size / len(nodes)),
                "final_validation_fraction": format_float(final_validation_size / len(nodes)),
                "final_test_fraction": format_float(final_test_size / len(nodes)),
                "median_max_similarity_to_train": median_max_similarity_to_train,
                "target_test_fraction": format_float(target_test_fraction),
                "target_validation_fraction": format_float(target_validation_fraction),
            }
        )

    write_tsv(output_dir / "split_summary.tsv", SUMMARY_COLUMNS, summary_rows)
    return summary_rows


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_valid_nodes(args.valid_structures)
    node_to_index = {(node.pdbid, node.chain): index for index, node in enumerate(nodes)}
    adjacency, observed_scores = load_similarity_adjacency(
        args.pairwise_scores,
        node_to_index=node_to_index,
        score_column=args.score_column,
    )
    threshold_specs = resolve_effective_thresholds(args, observed_scores)

    aggregate_rows: list[dict[str, str]] = []
    for seed in args.seeds:
        aggregate_rows.extend(
            run_for_seed(
                seed=seed,
                nodes=nodes,
                adjacency=adjacency,
                threshold_specs=threshold_specs,
                target_test_fraction=args.target_test_fraction,
                target_validation_fraction=args.target_validation_fraction,
                base_output_dir=args.output_dir,
            )
        )

    write_tsv(args.output_dir / "all_seed_split_summary.tsv", AGGREGATE_SUMMARY_COLUMNS, aggregate_rows)
    write_aggregate_summary_text(args.output_dir / "all_seed_split_summary.txt", aggregate_rows)


if __name__ == "__main__":
    main()
