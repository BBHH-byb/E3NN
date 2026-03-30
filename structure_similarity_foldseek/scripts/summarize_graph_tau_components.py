#!/usr/bin/env python3
"""
Example run:

PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/summarize_graph_tau_components.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv \
    --score_column tm_score_max \
    --graph_percentiles 0 20 40 60 80 90 92 94 95 96 97 98 99 \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/graph_tau_exploration

Alternative:

PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/summarize_graph_tau_components.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score/pairwise_foldseek.tsv \
    --score_column tm_score_max \
    --graph_taus 0.50 0.60 0.70 \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/graph_tau_exploration
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
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/graph_tau_exploration"
)
SUMMARY_COLUMNS = [
    "graph_percentile",
    "graph_tau",
    "n_proteins",
    "observed_pairs",
    "retained_edges",
    "n_components",
    "n_singletons",
    "largest_component_size",
    "mean_component_size",
    "median_component_size",
    "non_singleton_component_count",
]


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
        description="Summarize connected-component structure across graph tau percentiles or explicit graph tau values using tm_score_max."
    )
    parser.add_argument("--valid_structures", type=Path, default=DEFAULT_VALID_STRUCTURES)
    parser.add_argument("--pairwise_scores", type=Path, default=DEFAULT_PAIRWISE_SCORES)
    parser.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--score_column", default="tm_score_max")
    threshold_group = parser.add_mutually_exclusive_group()
    threshold_group.add_argument("--graph_percentiles", type=float, nargs="+")
    threshold_group.add_argument("--graph_taus", type=float, nargs="+")
    return parser.parse_args()


def format_float(value: float) -> str:
    return format(float(value), ".6f")


def percentile_suffix(percentile: float) -> str:
    return format(percentile, "g").replace(".0", "")


def resolve_graph_thresholds(args: argparse.Namespace, observed_scores: np.ndarray) -> list[tuple[str, float]]:
    if args.graph_taus is not None:
        return [("", float(value)) for value in args.graph_taus]
    percentiles = args.graph_percentiles
    if percentiles is None:
        percentiles = [0, 20, 40, 60, 80, 90, 92, 94, 95, 96, 97, 98, 99]
    return [(percentile_suffix(percentile), compute_percentile(observed_scores, percentile)) for percentile in percentiles]


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
) -> list[tuple[int, int, float]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pdbid1", "chain1", "pdbid2", "chain2", score_column}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing pairwise score columns: {', '.join(sorted(missing))}")

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
            seen_pairs.add(pair_key)
            score = float((row.get(score_column) or "").strip())
            edges.append((pair_key[0], pair_key[1], score))
    return edges


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


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_valid_nodes(args.valid_structures)
    node_to_index = {(node.pdbid, node.chain): index for index, node in enumerate(nodes)}
    scored_edges = load_pairwise_scores(args.pairwise_scores, node_to_index=node_to_index, score_column=args.score_column)
    observed_scores = np.asarray([score for _, _, score in scored_edges], dtype=float)

    summary_rows: list[dict[str, str]] = []
    for graph_percentile, graph_tau in resolve_graph_thresholds(args, observed_scores):
        retained_edges = [(left, right) for left, right, score in scored_edges if score >= graph_tau]
        components = compute_components(n_nodes=len(nodes), retained_edges=retained_edges)
        component_sizes = np.asarray([len(component) for component in components], dtype=float)
        summary_rows.append(
            {
                "graph_percentile": graph_percentile,
                "graph_tau": format_float(graph_tau),
                "n_proteins": str(len(nodes)),
                "observed_pairs": str(len(scored_edges)),
                "retained_edges": str(len(retained_edges)),
                "n_components": str(len(components)),
                "n_singletons": str(sum(1 for size in component_sizes if size == 1)),
                "largest_component_size": str(int(np.max(component_sizes)) if component_sizes.size else 0),
                "mean_component_size": format_float(np.mean(component_sizes) if component_sizes.size else 0.0),
                "median_component_size": format_float(np.median(component_sizes) if component_sizes.size else 0.0),
                "non_singleton_component_count": str(sum(1 for size in component_sizes if size > 1)),
            }
        )

    write_tsv(args.output_dir / "graph_tau_component_summary.tsv", SUMMARY_COLUMNS, summary_rows)


if __name__ == "__main__":
    main()
