#!/usr/bin/env python3
"""
Example run:

PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_ood_splits.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits \
    --graph_percentile 95 \
    --candidate_test_fraction 0.30 \
    --score_column tm_score_max \
    --effective_percentiles 0 20 40 60 80 90 92 94 95 96 97 98 99 \
    --target_test_fraction 0.10 \
    --seeds 0 1 2 3

Alternative:

PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/generate_ood_splits.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits \
    --graph_tau 0.657210 \
    --candidate_test_fraction 0.30 \
    --score_column tm_score_max \
    --taus 0.50 0.60 0.70 \
    --target_test_fraction 0.10 \
    --seeds 0 1 2 3

This will generate:

* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/split_p0.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/split_p20.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/split_p40.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/split_p60.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/split_p80.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/split_summary.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_1/split_summary.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_2/split_summary.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_3/split_summary.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/all_seed_split_summary.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/all_seed_split_summary.txt
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
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits"
)
SPLIT_COLUMNS = [
    "pdbid",
    "chain",
    "split",
    "component_id",
    "component_size",
    "graph_tau",
    "candidate_test",
    "effective_tau",
    "effective_percentile",
    "max_similarity_to_non_test",
]
SUMMARY_COLUMNS = [
    "graph_percentile",
    "graph_tau",
    "effective_percentile",
    "effective_tau",
    "n_proteins",
    "total_possible_pairs",
    "observed_pairs",
    "retained_graph_edges",
    "n_components",
    "largest_component_size",
    "n_singletons",
    "candidate_test_size",
    "candidate_test_fraction",
    "provisional_test_size",
    "initial_rejected_candidate_size",
    "iterative_removals",
    "final_train_size",
    "final_test_size",
    "final_train_fraction",
    "final_test_fraction",
    "target_test_fraction",
]
AGGREGATE_SUMMARY_COLUMNS = ["seed", *SUMMARY_COLUMNS]


@dataclass(frozen=True)
class Node:
    pdbid: str
    chain: str


class DisjointSet:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            self.parent[left_root] = right_root
        elif self.rank[left_root] > self.rank[right_root]:
            self.parent[right_root] = left_root
        else:
            self.parent[right_root] = left_root
            self.rank[left_root] += 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate TM-score-based OOD train/test splits using a fixed graph threshold, "
            "random candidate-test component selection, and iterative "
            "max-similarity-to-train filtering."
        )
    )
    parser.add_argument("--valid_structures", type=Path, default=DEFAULT_VALID_STRUCTURES)
    parser.add_argument("--pairwise_scores", type=Path, default=DEFAULT_PAIRWISE_SCORES)
    parser.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--score_column", default="tm_score_max")
    graph_group = parser.add_mutually_exclusive_group()
    graph_group.add_argument("--graph_percentile", type=float)
    graph_group.add_argument("--graph_tau", type=float)
    parser.add_argument("--candidate_test_fraction", type=float, default=0.30)
    effective_group = parser.add_mutually_exclusive_group()
    effective_group.add_argument("--effective_percentiles", type=float, nargs="+")
    effective_group.add_argument("--taus", type=float, nargs="+")
    parser.add_argument("--target_test_fraction", type=float, default=0.10)
    parser.add_argument("--seeds", type=int, nargs="+", default=[0])
    return parser.parse_args()


def format_float(value: float) -> str:
    return format(float(value), ".6f")


def percentile_suffix(percentile: float) -> str:
    return format(percentile, "g").replace(".0", "")


def resolve_graph_threshold(args: argparse.Namespace, observed_scores: np.ndarray) -> tuple[str, float]:
    if args.graph_tau is not None:
        return "", float(args.graph_tau)
    graph_percentile = 99.0 if args.graph_percentile is None else args.graph_percentile
    return percentile_suffix(graph_percentile), compute_percentile(observed_scores, graph_percentile)


def resolve_effective_thresholds(args: argparse.Namespace, observed_scores: np.ndarray) -> list[tuple[str, float]]:
    if args.taus is not None:
        return [("", float(value)) for value in args.taus]
    percentiles = args.effective_percentiles
    if percentiles is None:
        percentiles = [0, 20, 40, 60, 80, 90, 92, 94, 95, 96, 97, 98, 99]
    return [(percentile_suffix(percentile), compute_percentile(observed_scores, percentile)) for percentile in percentiles]


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


def load_pairwise_scores(
    path: Path,
    node_to_index: dict[tuple[str, str], int],
    score_column: str,
) -> tuple[list[tuple[int, int, float]], list[list[tuple[int, float]]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pdbid1", "chain1", "pdbid2", "chain2", score_column}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing pairwise score columns: {', '.join(sorted(missing))}")

        adjacency: list[list[tuple[int, float]]] = [[] for _ in range(len(node_to_index))]
        edges: list[tuple[int, int, float]] = []
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
            score = float((row.get(score_column) or "").strip())
            seen_pairs.add(pair_key)
            edges.append((pair_key[0], pair_key[1], score))
            adjacency[pair_key[0]].append((pair_key[1], score))
            adjacency[pair_key[1]].append((pair_key[0], score))
    return edges, adjacency


def compute_percentile(scores: np.ndarray, percentile: float) -> float:
    if scores.size == 0:
        return 0.0
    return float(np.percentile(scores, percentile))


def compute_components(n_nodes: int, retained_edges: list[tuple[int, int]]) -> list[list[int]]:
    dsu = DisjointSet(n_nodes)
    for left, right in retained_edges:
        dsu.union(left, right)

    groups: dict[int, list[int]] = {}
    for node_index in range(n_nodes):
        root = dsu.find(node_index)
        groups.setdefault(root, []).append(node_index)

    components = list(groups.values())
    components.sort(key=lambda comp: (-len(comp), min(comp)))
    return components


def build_component_lookup(components: list[list[int]]) -> dict[int, tuple[int, int]]:
    node_to_component: dict[int, tuple[int, int]] = {}
    for component_id, members in enumerate(components):
        component_size = len(members)
        for node_index in members:
            node_to_component[node_index] = (component_id, component_size)
    return node_to_component


def select_candidate_test_components(
    components: list[list[int]],
    total_nodes: int,
    candidate_test_fraction: float,
    random_seed: int,
) -> tuple[set[int], set[int], int]:
    if not components:
        return set(), set(), 0

    non_test_component_ids = {0}
    candidate_component_ids: set[int] = set()
    desired_candidate_size = candidate_test_fraction * total_nodes

    remaining_component_ids = list(range(1, len(components)))
    rng = np.random.default_rng(random_seed)
    if remaining_component_ids:
        shuffled = rng.permutation(remaining_component_ids).tolist()
    else:
        shuffled = []

    current_candidate_size = 0
    for component_id in shuffled:
        if current_candidate_size >= desired_candidate_size:
            break
        candidate_component_ids.add(component_id)
        current_candidate_size += len(components[component_id])

    for component_id in remaining_component_ids:
        if component_id not in candidate_component_ids:
            non_test_component_ids.add(component_id)

    return non_test_component_ids, candidate_component_ids, current_candidate_size


def max_similarity_to_members(
    node_index: int,
    member_mask: np.ndarray,
    adjacency: list[list[tuple[int, float]]],
) -> float:
    best = 0.0
    for neighbor_index, score in adjacency[node_index]:
        if member_mask[neighbor_index] and score > best:
            best = score
    return best


def compute_candidate_max_similarity_to_non_test(
    candidate_indices: list[int],
    non_test_mask: np.ndarray,
    adjacency: list[list[tuple[int, float]]],
) -> dict[int, float]:
    return {
        node_index: max_similarity_to_members(node_index=node_index, member_mask=non_test_mask, adjacency=adjacency)
        for node_index in candidate_indices
    }


def iterative_filter_test_set(
    candidate_indices: list[int],
    candidate_max_to_non_test: dict[int, float],
    effective_tau: float,
    non_test_mask: np.ndarray,
    adjacency: list[list[tuple[int, float]]],
) -> tuple[set[int], set[int], int]:
    accepted_test: set[int] = {
        node_index
        for node_index in candidate_indices
        if candidate_max_to_non_test[node_index] <= effective_tau
    }
    initially_rejected: set[int] = set(candidate_indices).difference(accepted_test)

    train_mask = non_test_mask.copy()
    for node_index in initially_rejected:
        train_mask[node_index] = True

    iterative_removals = 0
    while accepted_test:
        violating_nodes: list[int] = []
        for node_index in sorted(accepted_test):
            max_similarity_to_train = max_similarity_to_members(
                node_index=node_index,
                member_mask=train_mask,
                adjacency=adjacency,
            )
            if max_similarity_to_train > effective_tau:
                violating_nodes.append(node_index)
        if not violating_nodes:
            break
        for node_index in violating_nodes:
            accepted_test.remove(node_index)
            train_mask[node_index] = True
        iterative_removals += len(violating_nodes)

    final_train: set[int] = {index for index, is_train in enumerate(train_mask) if is_train}
    return final_train, accepted_test, iterative_removals


def build_split_rows(
    nodes: list[Node],
    node_to_component: dict[int, tuple[int, int]],
    candidate_test_mask: np.ndarray,
    candidate_max_to_non_test: dict[int, float],
    final_test: set[int],
    graph_tau: float,
    effective_tau: float,
    effective_percentile_label: str,
) -> tuple[list[dict[str, str]], dict[str, float]]:
    rows: list[dict[str, str]] = []
    final_test_size = 0
    final_train_size = 0

    for node_index, node in enumerate(nodes):
        component_id, component_size = node_to_component[node_index]
        split = "test" if node_index in final_test else "train"
        if split == "test":
            final_test_size += 1
        else:
            final_train_size += 1

        max_to_non_test = candidate_max_to_non_test.get(node_index)
        rows.append(
            {
                "pdbid": node.pdbid,
                "chain": node.chain,
                "split": split,
                "component_id": str(component_id),
                "component_size": str(component_size),
                "graph_tau": format_float(graph_tau),
                "candidate_test": "1" if candidate_test_mask[node_index] else "0",
                "effective_tau": format_float(effective_tau),
                "effective_percentile": effective_percentile_label,
                "max_similarity_to_non_test": "" if max_to_non_test is None else format_float(max_to_non_test),
            }
        )

    return rows, {
        "final_train_size": final_train_size,
        "final_test_size": final_test_size,
        "final_train_fraction": final_train_size / len(nodes),
        "final_test_fraction": final_test_size / len(nodes),
    }


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

    lines = ["Seed-by-Seed Effective Tau Summary", ""]
    def sort_key(value: str) -> float:
        if value.startswith("tau="):
            return float(value.split("=", 1)[1])
        return float(value)

    for effective_key in sorted(grouped, key=sort_key):
        rows = sorted(grouped[effective_key], key=lambda item: int(item["seed"]))
        fractions = np.asarray([float(row["final_test_fraction"]) for row in rows], dtype=float)
        if effective_key.startswith("tau="):
            lines.append(f"effective_tau: {effective_key.split('=', 1)[1]}")
        else:
            lines.append(f"effective_percentile: {effective_key}")
        lines.append(f"seed_count: {len(rows)}")
        lines.append(f"final_test_fraction_min: {format_float(np.min(fractions))}")
        lines.append(f"final_test_fraction_median: {format_float(np.median(fractions))}")
        lines.append(f"final_test_fraction_max: {format_float(np.max(fractions))}")
        for row in rows:
            lines.append(
                "  "
                f"seed={row['seed']} "
                f"graph_tau={row['graph_tau']} "
                f"effective_tau={row['effective_tau']} "
                f"candidate_test_fraction={row['candidate_test_fraction']} "
                f"final_test_fraction={row['final_test_fraction']} "
                f"final_test_size={row['final_test_size']}"
            )
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def run_for_seed(
    *,
    args: argparse.Namespace,
    nodes: list[Node],
    scored_edges: list[tuple[int, int, float]],
    adjacency: list[list[tuple[int, float]]],
    observed_scores: np.ndarray,
    seed: int,
) -> list[dict[str, str]]:
    output_dir = resolve_seed_output_dir(args.output_dir, seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_nodes = len(nodes)
    total_pairs = n_nodes * (n_nodes - 1) // 2
    graph_percentile_label, graph_tau = resolve_graph_threshold(args, observed_scores)
    retained_edges = [(left, right) for left, right, score in scored_edges if score >= graph_tau]
    components = compute_components(n_nodes=n_nodes, retained_edges=retained_edges)
    node_to_component = build_component_lookup(components)
    component_sizes = [len(component) for component in components]
    largest_component_size = max(component_sizes) if component_sizes else 0
    n_singletons = sum(1 for size in component_sizes if size == 1)

    non_test_component_ids, candidate_component_ids, candidate_test_size = select_candidate_test_components(
        components=components,
        total_nodes=n_nodes,
        candidate_test_fraction=args.candidate_test_fraction,
        random_seed=seed,
    )

    candidate_test_mask = np.zeros(n_nodes, dtype=bool)
    non_test_mask = np.zeros(n_nodes, dtype=bool)
    candidate_indices: list[int] = []
    for node_index in range(n_nodes):
        component_id, _ = node_to_component[node_index]
        if component_id in candidate_component_ids:
            candidate_test_mask[node_index] = True
            candidate_indices.append(node_index)
        elif component_id in non_test_component_ids:
            non_test_mask[node_index] = True
        else:
            raise RuntimeError(f"node {node_index} belongs to an unassigned component")

    candidate_max_to_non_test = compute_candidate_max_similarity_to_non_test(
        candidate_indices=candidate_indices,
        non_test_mask=non_test_mask,
        adjacency=adjacency,
    )

    summary_rows: list[dict[str, str]] = []
    for effective_percentile_label, effective_tau in resolve_effective_thresholds(args, observed_scores):
        final_train, final_test, iterative_removals = iterative_filter_test_set(
            candidate_indices=candidate_indices,
            candidate_max_to_non_test=candidate_max_to_non_test,
            effective_tau=effective_tau,
            non_test_mask=non_test_mask,
            adjacency=adjacency,
        )
        split_rows, split_stats = build_split_rows(
            nodes=nodes,
            node_to_component=node_to_component,
            candidate_test_mask=candidate_test_mask,
            candidate_max_to_non_test=candidate_max_to_non_test,
            final_test=final_test,
            graph_tau=graph_tau,
            effective_tau=effective_tau,
            effective_percentile_label=effective_percentile_label,
        )

        accepted_initial_size = sum(
            1 for node_index in candidate_indices if candidate_max_to_non_test[node_index] <= effective_tau
        )
        initial_rejected_size = len(candidate_indices) - accepted_initial_size

        file_suffix = effective_percentile_label if effective_percentile_label else format_float(effective_tau).replace(".", "p")
        write_tsv(output_dir / f"split_p{file_suffix}.tsv", SPLIT_COLUMNS, split_rows)

        summary_rows.append(
            {
                "graph_percentile": graph_percentile_label,
                "graph_tau": format_float(graph_tau),
                "effective_percentile": effective_percentile_label,
                "effective_tau": format_float(effective_tau),
                "n_proteins": str(n_nodes),
                "total_possible_pairs": str(total_pairs),
                "observed_pairs": str(len(scored_edges)),
                "retained_graph_edges": str(len(retained_edges)),
                "n_components": str(len(components)),
                "largest_component_size": str(largest_component_size),
                "n_singletons": str(n_singletons),
                "candidate_test_size": str(candidate_test_size),
                "candidate_test_fraction": format_float(candidate_test_size / n_nodes),
                "provisional_test_size": str(accepted_initial_size),
                "initial_rejected_candidate_size": str(initial_rejected_size),
                "iterative_removals": str(iterative_removals),
                "final_train_size": str(split_stats["final_train_size"]),
                "final_test_size": str(split_stats["final_test_size"]),
                "final_train_fraction": format_float(split_stats["final_train_fraction"]),
                "final_test_fraction": format_float(split_stats["final_test_fraction"]),
                "target_test_fraction": format_float(args.target_test_fraction),
            }
        )

        if len(final_train) + len(final_test) != n_nodes:
            raise RuntimeError("final train/test partition does not cover all nodes")

    write_tsv(output_dir / "split_summary.tsv", SUMMARY_COLUMNS, summary_rows)
    return summary_rows


def main() -> None:
    args = parse_args()
    nodes = load_valid_nodes(args.valid_structures)
    node_to_index = {(node.pdbid, node.chain): index for index, node in enumerate(nodes)}
    scored_edges, adjacency = load_pairwise_scores(
        args.pairwise_scores,
        node_to_index=node_to_index,
        score_column=args.score_column,
    )
    observed_scores = np.asarray([score for _, _, score in scored_edges], dtype=float)
    aggregate_rows: list[dict[str, str]] = []
    for seed in args.seeds:
        summary_rows = run_for_seed(
            args=args,
            nodes=nodes,
            scored_edges=scored_edges,
            adjacency=adjacency,
            observed_scores=observed_scores,
            seed=seed,
        )
        for row in summary_rows:
            aggregate_rows.append({"seed": str(seed), **row})

    write_tsv(args.output_dir / "all_seed_split_summary.tsv", AGGREGATE_SUMMARY_COLUMNS, aggregate_rows)
    write_aggregate_summary_text(args.output_dir / "all_seed_split_summary.txt", aggregate_rows)


if __name__ == "__main__":
    main()
