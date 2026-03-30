#!/usr/bin/env python3
"""
Example run:
PYTHONDONTWRITEBYTECODE=1 python3 \
   /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/validate_ood_splits.py \
   --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/valid_structures.tsv \
   --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/pairwise_foldseek.tsv \
   --split_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits \
   --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits \
   --percentiles 0 20 40 60 80 85 90 92 94 95 96 97 98 99
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


sys.dont_write_bytecode = True


DEFAULT_VALID_STRUCTURES = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/valid_structures.tsv"
)
DEFAULT_PAIRWISE_SCORES = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/pairwise_foldseek.tsv"
)
DEFAULT_SPLIT_DIR = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits"
)
SAMPLE_COLUMNS = [
    "effective_percentile",
    "pdbid",
    "chain",
    "component_id",
    "component_size",
    "is_singleton",
    "max_similarity_to_train",
    "nearest_train_pdbid",
    "nearest_train_chain",
]
SUMMARY_COLUMNS = [
    "percentile",
    "n_test",
    "min",
    "q1",
    "median",
    "mean",
    "q3",
    "p90",
    "p95",
    "max",
    "n_singletons_in_test",
]
EMBED_COLUMNS = ["pdbid", "chain", "umap1", "umap2"]


@dataclass(frozen=True)
class Node:
    pdbid: str
    chain: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate Foldseek-based OOD splits using test-to-train maximum similarity statistics and plots."
    )
    parser.add_argument("--valid_structures", type=Path, default=DEFAULT_VALID_STRUCTURES)
    parser.add_argument("--pairwise_scores", type=Path, default=DEFAULT_PAIRWISE_SCORES)
    parser.add_argument("--split_dir", type=Path, default=DEFAULT_SPLIT_DIR)
    parser.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--score_column", default="foldseek_score_max")
    parser.add_argument("--percentiles", type=float, nargs="+", default=[0, 20, 40, 60, 80])
    parser.add_argument("--random_seed", type=int, default=0)
    parser.add_argument("--skip_umap", action="store_true")
    return parser.parse_args()


def format_float(value: float) -> str:
    return format(float(value), ".6f")


def resolve_seed_split_dir(base_split_dir: Path, random_seed: int) -> Path:
    seed_name = f"random_seed_{random_seed}"
    if base_split_dir.name == seed_name:
        return base_split_dir
    return base_split_dir / seed_name


def resolve_seed_validation_dir(base_output_dir: Path, random_seed: int) -> Path:
    seed_name = f"random_seed_{random_seed}"
    if base_output_dir.name == "validation" and base_output_dir.parent.name == seed_name:
        return base_output_dir
    if base_output_dir.name == seed_name:
        return base_output_dir / "validation"
    return base_output_dir / seed_name / "validation"


def ensure_matplotlib(output_dir: Path):
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".mplconfig"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def load_nodes(path: Path) -> list[Node]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pdbid", "chain"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing valid_structures columns: {', '.join(sorted(missing))}")
        nodes = []
        seen = set()
        for row in reader:
            key = ((row.get("pdbid") or "").strip(), (row.get("chain") or "").strip())
            if not key[0] or not key[1]:
                continue
            if key in seen:
                raise ValueError(f"duplicate node: {key[0]}:{key[1]}")
            seen.add(key)
            nodes.append(Node(*key))
    return nodes


def load_similarity_matrix(path: Path, node_to_index: dict[tuple[str, str], int], score_column: str) -> np.ndarray:
    n_nodes = len(node_to_index)
    matrix = np.zeros((n_nodes, n_nodes), dtype=np.float32)
    np.fill_diagonal(matrix, 1.0)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pdbid1", "chain1", "pdbid2", "chain2", score_column}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing pairwise score columns: {', '.join(sorted(missing))}")
        for row in reader:
            left = ((row.get("pdbid1") or "").strip(), (row.get("chain1") or "").strip())
            right = ((row.get("pdbid2") or "").strip(), (row.get("chain2") or "").strip())
            if left not in node_to_index or right not in node_to_index:
                continue
            score = float((row.get(score_column) or "").strip())
            left_index = node_to_index[left]
            right_index = node_to_index[right]
            matrix[left_index, right_index] = score
            matrix[right_index, left_index] = score
    return matrix


def load_split_files(split_dir: Path, percentiles: list[float]) -> list[Path]:
    paths = []
    for percentile in percentiles:
        suffix = format(percentile, "g").replace(".0", "")
        path = split_dir / f"split_p{suffix}.tsv"
        if not path.exists():
            raise ValueError(f"missing split file: {path}")
        paths.append(path)
    if not paths:
        raise ValueError(f"no split_p*.tsv files found in {split_dir}")
    return paths


def load_split_rows(path: Path, node_to_index: dict[tuple[str, str], int]) -> tuple[str, list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pdbid", "chain", "split", "component_id", "component_size", "effective_percentile"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing split columns in {path}: {', '.join(sorted(missing))}")
        rows = list(reader)
    percentile = rows[0]["effective_percentile"] if rows else path.stem.removeprefix("split_p")
    for row in rows:
        key = ((row.get("pdbid") or "").strip(), (row.get("chain") or "").strip())
        if key not in node_to_index:
            raise ValueError(f"split row not found in valid structures: {key[0]}:{key[1]}")
    return percentile, rows


def compute_test_to_train(
    percentile: str,
    rows: list[dict[str, str]],
    node_to_index: dict[tuple[str, str], int],
    nodes: list[Node],
    similarity_matrix: np.ndarray,
) -> tuple[list[dict[str, str]], list[float]]:
    train_indices = []
    test_items = []
    for row in rows:
        key = ((row.get("pdbid") or "").strip(), (row.get("chain") or "").strip())
        index = node_to_index[key]
        split = (row.get("split") or "").strip()
        if split == "train":
            train_indices.append(index)
        elif split == "test":
            test_items.append((index, row))

    if not test_items:
        return [], []
    if not train_indices:
        raise ValueError(f"split_p{percentile} has no train samples")

    sample_rows: list[dict[str, str]] = []
    max_values: list[float] = []
    train_indices_array = np.asarray(train_indices, dtype=int)

    for test_index, row in test_items:
        sims = similarity_matrix[test_index, train_indices_array]
        best_position = int(np.argmax(sims))
        best_index = int(train_indices_array[best_position])
        best_score = float(sims[best_position])
        max_values.append(best_score)
        sample_rows.append(
            {
                "effective_percentile": percentile,
                "pdbid": nodes[test_index].pdbid,
                "chain": nodes[test_index].chain,
                "component_id": (row.get("component_id") or "").strip(),
                "component_size": (row.get("component_size") or "").strip(),
                "is_singleton": "1" if (row.get("component_size") or "").strip() == "1" else "0",
                "max_similarity_to_train": format_float(best_score),
                "nearest_train_pdbid": nodes[best_index].pdbid,
                "nearest_train_chain": nodes[best_index].chain,
            }
        )

    return sample_rows, max_values


def summarize_distribution(values: list[float], sample_rows: list[dict[str, str]], percentile: str) -> dict[str, str]:
    array = np.asarray(values, dtype=float)
    return {
        "percentile": percentile,
        "n_test": str(len(values)),
        "min": format_float(np.min(array)),
        "q1": format_float(np.percentile(array, 25)),
        "median": format_float(np.median(array)),
        "mean": format_float(np.mean(array)),
        "q3": format_float(np.percentile(array, 75)),
        "p90": format_float(np.percentile(array, 90)),
        "p95": format_float(np.percentile(array, 95)),
        "max": format_float(np.max(array)),
        "n_singletons_in_test": str(sum(1 for row in sample_rows if row["is_singleton"] == "1")),
    }


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary_text(path: Path, summary_rows: list[dict[str, str]]) -> None:
    lines = ["Test-to-Train Maximum Similarity Summary", ""]
    for row in summary_rows:
        lines.extend(
            [
                f"percentile: {row['percentile']}",
                f"n_test: {row['n_test']}",
                f"min: {row['min']}",
                f"q1: {row['q1']}",
                f"median: {row['median']}",
                f"mean: {row['mean']}",
                f"q3: {row['q3']}",
                f"p90: {row['p90']}",
                f"p95: {row['p95']}",
                f"max: {row['max']}",
                f"n_singletons_in_test: {row['n_singletons_in_test']}",
                "",
            ]
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def make_violin_plot(output_path: Path, values_by_percentile: list[tuple[str, list[float]]], plt) -> None:
    labels = [label for label, _ in values_by_percentile]
    values = [vals for _, vals in values_by_percentile]
    fig, ax = plt.subplots(figsize=(8, 5))
    parts = ax.violinplot(values, showmeans=True, showmedians=True)
    for body in parts["bodies"]:
        body.set_alpha(0.7)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlabel("Percentile")
    ax.set_ylabel("Max Similarity to Train")
    ax.set_title("Test-to-Train Maximum Similarity")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def make_ecdf_plot(output_path: Path, values_by_percentile: list[tuple[str, list[float]]], plt) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, values in values_by_percentile:
        array = np.sort(np.asarray(values, dtype=float))
        y = np.arange(1, len(array) + 1, dtype=float) / len(array)
        ax.step(array, y, where="post", label=f"p{label}")
    ax.set_xlabel("Max Similarity to Train")
    ax.set_ylabel("ECDF")
    ax.set_title("ECDF of Test Maximum Similarity")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def compute_umap_embedding(similarity_matrix: np.ndarray, random_seed: int) -> np.ndarray:
    try:
        import umap
    except Exception as exc:
        raise RuntimeError(
            "UMAP plotting requires umap-learn. Install 'umap-learn' to generate UMAP visualizations."
        ) from exc

    distance_matrix = np.asarray(1.0 - similarity_matrix, dtype=np.float32)
    np.fill_diagonal(distance_matrix, 0.0)
    reducer = umap.UMAP(metric="precomputed", random_state=random_seed)
    embedding = reducer.fit_transform(distance_matrix)
    return np.asarray(embedding, dtype=float)


def write_umap_embedding(path: Path, nodes: list[Node], embedding: np.ndarray) -> None:
    rows = [
        {
            "pdbid": node.pdbid,
            "chain": node.chain,
            "umap1": format_float(point[0]),
            "umap2": format_float(point[1]),
        }
        for node, point in zip(nodes, embedding, strict=True)
    ]
    write_tsv(path, EMBED_COLUMNS, rows)


def plot_umap_split(
    output_path: Path,
    embedding: np.ndarray,
    rows: list[dict[str, str]],
    node_to_index: dict[tuple[str, str], int],
    percentile: str,
    plt,
) -> None:
    colors = {"train": "#4C78A8", "test": "#E45756"}
    fig, ax = plt.subplots(figsize=(8, 6))
    for split_name in ["train", "test"]:
        indices = []
        for row in rows:
            if (row.get("split") or "").strip() != split_name:
                continue
            key = ((row.get("pdbid") or "").strip(), (row.get("chain") or "").strip())
            indices.append(node_to_index[key])
        if not indices:
            continue
        points = embedding[np.asarray(indices, dtype=int)]
        ax.scatter(points[:, 0], points[:, 1], s=10, alpha=0.8, label=split_name, color=colors[split_name])
    ax.set_title(f"UMAP Split Visualization (p{percentile})")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    split_dir = resolve_seed_split_dir(args.split_dir, args.random_seed)
    output_dir = resolve_seed_validation_dir(args.output_dir, args.random_seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_nodes(args.valid_structures)
    node_to_index = {(node.pdbid, node.chain): index for index, node in enumerate(nodes)}
    similarity_matrix = load_similarity_matrix(
        args.pairwise_scores,
        node_to_index=node_to_index,
        score_column=args.score_column,
    )

    split_paths = load_split_files(split_dir, percentiles=args.percentiles)
    all_summary_rows: list[dict[str, str]] = []
    values_by_percentile: list[tuple[str, list[float]]] = []
    split_rows_by_percentile: dict[str, list[dict[str, str]]] = {}

    for split_path in split_paths:
        percentile, split_rows = load_split_rows(split_path, node_to_index=node_to_index)
        sample_rows, values = compute_test_to_train(
            percentile=percentile,
            rows=split_rows,
            node_to_index=node_to_index,
            nodes=nodes,
            similarity_matrix=similarity_matrix,
        )
        write_tsv(output_dir / f"split_p{percentile}_test_vs_train.tsv", SAMPLE_COLUMNS, sample_rows)
        if values:
            all_summary_rows.append(summarize_distribution(values, sample_rows, percentile))
            values_by_percentile.append((percentile, values))
        split_rows_by_percentile[percentile] = split_rows

    write_tsv(output_dir / "test_max_similarity_summary.tsv", SUMMARY_COLUMNS, all_summary_rows)
    write_summary_text(output_dir / "test_max_similarity_summary.txt", all_summary_rows)

    plt = ensure_matplotlib(output_dir)
    make_violin_plot(output_dir / "test_max_similarity_violin.png", values_by_percentile, plt)
    make_ecdf_plot(output_dir / "test_max_similarity_ecdf.png", values_by_percentile, plt)

    if not args.skip_umap:
        embedding = compute_umap_embedding(similarity_matrix, random_seed=args.random_seed)
        write_umap_embedding(output_dir / "umap_embedding.tsv", nodes, embedding)
        for percentile, split_rows in split_rows_by_percentile.items():
            plot_umap_split(
                output_dir / f"umap_split_p{percentile}.png",
                embedding=embedding,
                rows=split_rows,
                node_to_index=node_to_index,
                percentile=percentile,
                plt=plt,
            )


if __name__ == "__main__":
    main()
