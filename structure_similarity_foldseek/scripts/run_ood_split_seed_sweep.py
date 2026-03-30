#!/usr/bin/env python3
"""
Example run:

PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/run_ood_split_seed_sweep.py \
    --valid_structures /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/valid_structures.tsv \
    --pairwise_scores /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/pairwise_foldseek.tsv \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits \
    --graph_percentile 99 \
    --candidate_test_fraction 0.30 \
    --percentiles 90 92 94 95 96 97 98 99 \
    --target_test_fraction 0.10 \
    --seeds 0 1 2 3 \
    --skip_umap

This will generate and validate:

* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_0/validation/
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_1/
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_1/validation/
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_2/
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_2/validation/
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_3/
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/random_seed_3/validation/

It will also write:

* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/seed_sweep_test_fraction.tsv
* /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits/seed_sweep_test_fraction.png
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path


sys.dont_write_bytecode = True


DEFAULT_VALID_STRUCTURES = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/valid_structures.tsv"
)
DEFAULT_PAIRWISE_SCORES = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek/pairwise_foldseek.tsv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/ood_splits"
)
SCRIPT_DIR = Path(__file__).resolve().parent
GENERATE_SCRIPT = SCRIPT_DIR / "generate_ood_splits.py"
VALIDATE_SCRIPT = SCRIPT_DIR / "validate_ood_splits.py"
AGGREGATE_COLUMNS = [
    "seed",
    "effective_percentile",
    "effective_tau",
    "final_test_size",
    "final_test_fraction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run generate_ood_splits.py and validate_ood_splits.py for multiple random seeds, "
            "writing each seed to its own random_seed_<n> directory."
        )
    )
    parser.add_argument("--valid_structures", type=Path, default=DEFAULT_VALID_STRUCTURES)
    parser.add_argument("--pairwise_scores", type=Path, default=DEFAULT_PAIRWISE_SCORES)
    parser.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--graph_percentile", type=float, default=99.0)
    parser.add_argument("--candidate_test_fraction", type=float, default=0.30)
    parser.add_argument("--percentiles", type=float, nargs="+", default=[90, 95, 97, 99])
    parser.add_argument("--target_test_fraction", type=float, default=0.10)
    parser.add_argument("--score_column", default="foldseek_score_max")
    parser.add_argument("--seeds", type=int, nargs="+", required=True)
    parser.add_argument("--skip_umap", action="store_true")
    return parser.parse_args()


def percentile_args(percentiles: list[float]) -> list[str]:
    return [format(percentile, "g") for percentile in percentiles]


def run_command(command: list[str]) -> None:
    print("$ " + " ".join(command), flush=True)
    env = os.environ.copy()
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    subprocess.run(command, check=True, env=env)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ensure_matplotlib(output_dir: Path):
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".mplconfig"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def collect_seed_summaries(base_output_dir: Path, seeds: list[int]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for seed in seeds:
        summary_path = base_output_dir / f"random_seed_{seed}" / "split_summary.tsv"
        with summary_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"effective_percentile", "effective_tau", "final_test_size", "final_test_fraction"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(f"missing columns in {summary_path}: {', '.join(sorted(missing))}")
            for row in reader:
                rows.append(
                    {
                        "seed": str(seed),
                        "effective_percentile": (row.get("effective_percentile") or "").strip(),
                        "effective_tau": (row.get("effective_tau") or "").strip(),
                        "final_test_size": (row.get("final_test_size") or "").strip(),
                        "final_test_fraction": (row.get("final_test_fraction") or "").strip(),
                    }
                )
    rows.sort(key=lambda item: (int(item["seed"]), float(item["effective_percentile"])))
    return rows


def make_seed_sweep_plot(output_path: Path, summary_rows: list[dict[str, str]], plt) -> None:
    by_seed: dict[str, list[tuple[float, float]]] = {}
    for row in summary_rows:
        seed = row["seed"]
        percentile = float(row["effective_percentile"])
        fraction = float(row["final_test_fraction"])
        by_seed.setdefault(seed, []).append((percentile, fraction))

    fig, ax = plt.subplots(figsize=(8, 5))
    for seed, points in sorted(by_seed.items(), key=lambda item: int(item[0])):
        points.sort(key=lambda item: item[0])
        ax.plot(
            [percentile for percentile, _ in points],
            [fraction for _, fraction in points],
            marker="o",
            linewidth=1.5,
            label=f"seed {seed}",
        )
    ax.set_xlabel("Effective Percentile")
    ax.set_ylabel("Final Test Fraction")
    ax.set_title("Final Test Fraction by Seed and Effective Percentile")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    percentiles = percentile_args(args.percentiles)

    for seed in args.seeds:
        print(f"[seed {seed}] generate", flush=True)
        generate_command = [
            sys.executable,
            str(GENERATE_SCRIPT),
            "--valid_structures",
            str(args.valid_structures),
            "--pairwise_scores",
            str(args.pairwise_scores),
            "--output_dir",
            str(args.output_dir),
            "--graph_percentile",
            format(args.graph_percentile, "g"),
            "--candidate_test_fraction",
            format(args.candidate_test_fraction, "g"),
            "--percentiles",
            *percentiles,
            "--target_test_fraction",
            format(args.target_test_fraction, "g"),
            "--random_seed",
            str(seed),
        ]
        run_command(generate_command)

        print(f"[seed {seed}] validate", flush=True)
        validate_command = [
            sys.executable,
            str(VALIDATE_SCRIPT),
            "--valid_structures",
            str(args.valid_structures),
            "--pairwise_scores",
            str(args.pairwise_scores),
            "--split_dir",
            str(args.output_dir),
            "--output_dir",
            str(args.output_dir),
            "--score_column",
            args.score_column,
            "--percentiles",
            *percentiles,
            "--random_seed",
            str(seed),
        ]
        if args.skip_umap:
            validate_command.append("--skip_umap")
        run_command(validate_command)

    aggregate_rows = collect_seed_summaries(args.output_dir, args.seeds)
    write_tsv(args.output_dir / "seed_sweep_test_fraction.tsv", AGGREGATE_COLUMNS, aggregate_rows)
    plt = ensure_matplotlib(args.output_dir)
    make_seed_sweep_plot(args.output_dir / "seed_sweep_test_fraction.png", aggregate_rows, plt)


if __name__ == "__main__":
    main()
