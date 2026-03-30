#!/usr/bin/env python3
"""
Example run:

PYTHONDONTWRITEBYTECODE=1 python3 \
    /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/scripts/run_pairwise_foldseek.py \
    --input_csv /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/csv_input/input_monomer_dataset.csv \
    --output_dir /Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek_based_TM_score \
    --foldseek_path /opt/homebrew/bin/foldseek \
    --tmalign_path /Users/byulee/PROJECT/E3NN/TMalign/TMalign
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np


sys.dont_write_bytecode = True


DEFAULT_INPUT_CSV = Path("/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/csv_input/input_monomer_dataset.csv")
DEFAULT_OUTPUT_DIR = Path(
    "/Users/byulee/PROJECT/E3NN/NMR_dataset_analysis/structure_similarity_foldseek/results_foldseek"
)
VALID_COLUMNS = ["pdbid", "chain", "pdb_path", "nresid"]
INVALID_COLUMNS = ["pdbid", "chain", "reason"]
RAW_HIT_COLUMNS = [
    "query_pdbid",
    "query_chain",
    "target_pdbid",
    "target_chain",
    "primary_score",
    "alntmscore",
    "qtmscore",
    "ttmscore",
    "evalue",
    "bits",
    "fident",
    "alnlen",
    "prob",
]
PAIR_COLUMNS = [
    "pdbid1",
    "chain1",
    "pdbid2",
    "chain2",
    "foldseek_score_max",
    "foldseek_score_min",
    "tm_score_max",
    "tm_score_min",
    "hit_direction_count",
]
FOLDSEEK_FIELDS = [
    "query",
    "target",
    "alntmscore",
    "qtmscore",
    "ttmscore",
    "evalue",
    "bits",
    "fident",
    "alnlen",
    "prob",
]
PRIMARY_SCORE_FIELD = "alntmscore"
TM_SCORE_PATTERN = re.compile(r"TM-score=\s*([0-9.]+)")


@dataclass(frozen=True)
class StructureRecord:
    pdbid: str
    chain: str
    nresid: int
    coordinates: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run all-vs-all Foldseek structural similarity on the NMR dataset."
    )
    parser.add_argument("--input_csv", type=Path, default=DEFAULT_INPUT_CSV)
    parser.add_argument("--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--foldseek_path", type=Path, required=True)
    parser.add_argument("--tmalign_path", type=Path, required=True)
    parser.add_argument(
        "--alignment_type",
        type=int,
        default=1,
        help="Foldseek alignment type. Default 1 uses TMalign realignment.",
    )
    parser.add_argument(
        "--sensitivity",
        type=float,
        default=None,
        help="Optional Foldseek -s sensitivity value.",
    )
    parser.add_argument(
        "--max_seqs",
        type=int,
        default=None,
        help="Optional Foldseek --max-seqs value.",
    )
    parser.add_argument(
        "--progress_interval",
        type=int,
        default=500,
        help="Report CSV parsing/materialization progress every N rows.",
    )
    return parser.parse_args()


def format_float(value: float) -> str:
    return format(float(value), ".6f")


def format_duration(seconds: float) -> str:
    total_seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds_part = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds_part:02d}"


def print_progress(message: str) -> None:
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def count_data_rows(input_csv: Path) -> int:
    with input_csv.open("r", encoding="utf-8", newline="") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def estimate_eta(start_time: float, completed: int, total: int) -> str:
    if completed <= 0 or total <= 0:
        return "unknown"
    elapsed = time.time() - start_time
    rate = elapsed / completed
    remaining = max(total - completed, 0) * rate
    return format_duration(remaining)


def parse_coordinates(raw_value: str, nresid: int) -> np.ndarray:
    parts = [part.strip() for part in raw_value.split(",") if part.strip()]
    expected = 3 * nresid
    if len(parts) != expected:
        raise ValueError(f"coordinate length mismatch: expected={expected}, found={len(parts)}")

    try:
        values = np.asarray([float(part) for part in parts], dtype=float)
    except ValueError as exc:
        raise ValueError(f"malformed coordinate data: {exc}") from exc

    if not np.isfinite(values).all():
        raise ValueError("coordinate contains non-finite values")
    return values.reshape(nresid, 3)


def choose_pdb_chain_id(chain: str) -> str:
    return chain if len(chain) == 1 and chain.isalnum() else "A"


def write_pseudo_pdb(record: StructureRecord, pdb_path: Path) -> None:
    chain_id = choose_pdb_chain_id(record.chain)
    with pdb_path.open("w", encoding="utf-8") as handle:
        for residue_index, (x, y, z) in enumerate(record.coordinates, start=1):
            handle.write(
                f"ATOM  {residue_index:5d}  CA  ALA {chain_id:1}{residue_index:4d}"
                f"    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )
        handle.write("END\n")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_and_materialize_structures(
    input_csv: Path,
    structures_dir: Path,
    progress_interval: int,
) -> tuple[list[StructureRecord], list[dict[str, str]], list[dict[str, str]], int]:
    valid_records: list[StructureRecord] = []
    valid_rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    total_rows = count_data_rows(input_csv)
    processed_rows = 0
    stage_start = time.time()

    print_progress(f"input_csv={input_csv}")
    print_progress(f"materializing CA-only structures total_rows={total_rows}")

    with input_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"pdbid", "chain", "nresid", "coordinate"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            processed_rows += 1
            pdbid = (row.get("pdbid") or "").strip()
            chain = (row.get("chain") or "").strip()
            try:
                nresid = int((row.get("nresid") or "").strip())
                if nresid <= 0:
                    raise ValueError("nresid must be positive")
                coordinates = parse_coordinates(row.get("coordinate") or "", nresid)
                record = StructureRecord(pdbid=pdbid, chain=chain, nresid=nresid, coordinates=coordinates)
                pdb_path = structures_dir / f"{pdbid}_{chain}.pdb"
                write_pseudo_pdb(record, pdb_path)
            except Exception as exc:
                invalid_rows.append({"pdbid": pdbid, "chain": chain, "reason": str(exc)})
                continue

            valid_records.append(record)
            valid_rows.append(
                {
                    "pdbid": pdbid,
                    "chain": chain,
                    "pdb_path": str(pdb_path.resolve()),
                    "nresid": str(nresid),
                }
            )

            if (
                processed_rows == total_rows
                or (progress_interval > 0 and processed_rows % progress_interval == 0)
            ):
                print_progress(
                    "structure_progress "
                    f"processed={processed_rows}/{total_rows} "
                    f"valid={len(valid_rows)} invalid={len(invalid_rows)} "
                    f"eta={estimate_eta(stage_start, processed_rows, total_rows)}"
                )

    return valid_records, valid_rows, invalid_rows, processed_rows


def run_command(command: list[str], log_lines: list[str], stage_name: str) -> None:
    print_progress(f"starting_stage={stage_name}")
    log_lines.append("COMMAND: " + " ".join(command))
    started = time.time()
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    captured_lines: list[str] = []
    assert process.stdout is not None
    for line in process.stdout:
        text = line.rstrip()
        captured_lines.append(text)
        print(f"[{stage_name}] {text}", flush=True)
    return_code = process.wait()
    if captured_lines:
        log_lines.append("OUTPUT:")
        log_lines.extend(captured_lines)
    print_progress(f"finished_stage={stage_name} elapsed={format_duration(time.time() - started)}")
    if return_code != 0:
        raise RuntimeError(
            f"command failed with exit code {return_code}: {' '.join(command)}"
        )


def build_foldseek_commands(
    foldseek_path: Path,
    structures_dir: Path,
    db_prefix: Path,
    result_db: Path,
    tmp_dir: Path,
    raw_hits_tmp: Path,
    alignment_type: int,
    sensitivity: float | None,
    max_seqs: int | None,
) -> list[list[str]]:
    createdb_command = [str(foldseek_path), "createdb", str(structures_dir), str(db_prefix)]
    search_command = [
        str(foldseek_path),
        "search",
        str(db_prefix),
        str(db_prefix),
        str(result_db),
        str(tmp_dir),
        "--alignment-type",
        str(alignment_type),
        "-a"
    ]
    if sensitivity is not None:
        search_command.extend(["-s", str(sensitivity)])
    if max_seqs is not None:
        search_command.extend(["--max-seqs", str(max_seqs)])

    convertalis_command = [
        str(foldseek_path),
        "convertalis",
        str(db_prefix),
        str(db_prefix),
        str(result_db),
        str(raw_hits_tmp),
        "--format-output",
        ",".join(FOLDSEEK_FIELDS),
    ]
    return [createdb_command, search_command, convertalis_command]


def parse_structure_id(raw_id: str) -> tuple[str, str]:
    value = Path(raw_id.strip()).stem
    if "_" not in value:
        raise ValueError(f"unable to split Foldseek structure ID: {raw_id!r}")
    pdbid, chain = value.split("_", 1)
    return pdbid, chain


def parse_foldseek_hits(raw_hits_tmp: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with raw_hits_tmp.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if not fields:
                continue
            if len(fields) != len(FOLDSEEK_FIELDS):
                raise ValueError(
                    f"unexpected Foldseek output column count: expected={len(FOLDSEEK_FIELDS)}, found={len(fields)}"
                )
            raw = dict(zip(FOLDSEEK_FIELDS, fields, strict=True))
            query_pdbid, query_chain = parse_structure_id(raw["query"])
            target_pdbid, target_chain = parse_structure_id(raw["target"])
            primary_score = float(raw[PRIMARY_SCORE_FIELD])
            rows.append(
                {
                    "query_pdbid": query_pdbid,
                    "query_chain": query_chain,
                    "target_pdbid": target_pdbid,
                    "target_chain": target_chain,
                    "primary_score": format_float(primary_score),
                    "alntmscore": raw["alntmscore"],
                    "qtmscore": raw["qtmscore"],
                    "ttmscore": raw["ttmscore"],
                    "evalue": raw["evalue"],
                    "bits": raw["bits"],
                    "fident": raw["fident"],
                    "alnlen": raw["alnlen"],
                    "prob": raw["prob"],
                }
            )
    return rows


def parse_tmalign_scores(stdout: str) -> tuple[float, float]:
    tm_scores = [float(match.group(1)) for match in TM_SCORE_PATTERN.finditer(stdout)]
    if len(tm_scores) < 2:
        raise ValueError("TM-align output parsing failure: TM-score not found twice")
    first, second = tm_scores[0], tm_scores[1]
    return max(first, second), min(first, second)


def canonical_pair(
    left: tuple[str, str],
    right: tuple[str, str],
) -> tuple[tuple[str, str], tuple[str, str]]:
    return (left, right) if left <= right else (right, left)


def build_pairwise_rows(raw_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[tuple[str, str], tuple[str, str]], dict[str, float]] = {}

    for row in raw_rows:
        left = (row["query_pdbid"], row["query_chain"])
        right = (row["target_pdbid"], row["target_chain"])
        if left == right:
            continue
        pair_left, pair_right = canonical_pair(left, right)
        direction_key = "forward" if (left, right) == (pair_left, pair_right) else "reverse"
        grouped.setdefault((pair_left, pair_right), {})
        score = float(row["primary_score"])
        previous = grouped[(pair_left, pair_right)].get(direction_key)
        if previous is None or score > previous:
            grouped[(pair_left, pair_right)][direction_key] = score

    pair_rows: list[dict[str, str]] = []
    for (left, right), direction_scores in sorted(grouped.items()):
        scores = list(direction_scores.values())
        pair_rows.append(
            {
                "pdbid1": left[0],
                "chain1": left[1],
                "pdbid2": right[0],
                "chain2": right[1],
                "foldseek_score_max": format_float(max(scores)),
                "foldseek_score_min": format_float(min(scores)),
                "tm_score_max": "",
                "tm_score_min": "",
                "hit_direction_count": str(len(direction_scores)),
            }
        )

    return pair_rows


def run_tmalign_rescoring(
    pair_rows: list[dict[str, str]],
    structures_dir: Path,
    tmalign_path: Path,
    progress_interval: int,
) -> list[dict[str, str]]:
    print_progress(f"starting_stage=tmalign_rescoring total_pairs={len(pair_rows)}")
    started = time.time()
    rescored_rows: list[dict[str, str]] = []

    for index, row in enumerate(pair_rows, start=1):
        left_path = structures_dir / f"{row['pdbid1']}_{row['chain1']}.pdb"
        right_path = structures_dir / f"{row['pdbid2']}_{row['chain2']}.pdb"
        if not left_path.exists():
            raise FileNotFoundError(f"TM-align input structure not found: {left_path}")
        if not right_path.exists():
            raise FileNotFoundError(f"TM-align input structure not found: {right_path}")

        completed = subprocess.run(
            [str(tmalign_path), str(left_path), str(right_path)],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            message = completed.stderr.strip() or completed.stdout.strip()
            raise RuntimeError(message or f"TM-align exited with code {completed.returncode}")

        tm_score_max, tm_score_min = parse_tmalign_scores(completed.stdout)
        updated_row = dict(row)
        updated_row["tm_score_max"] = format_float(tm_score_max)
        updated_row["tm_score_min"] = format_float(tm_score_min)
        rescored_rows.append(updated_row)

        if index == len(pair_rows) or (progress_interval > 0 and index % progress_interval == 0):
            print_progress(
                "tmalign_progress "
                f"processed={index}/{len(pair_rows)} "
                f"eta={estimate_eta(started, index, len(pair_rows))}"
            )

    print_progress(f"finished_stage=tmalign_rescoring elapsed={format_duration(time.time() - started)}")
    return rescored_rows


def compute_score_correlation_summary(pair_rows: list[dict[str, str]]) -> list[str]:
    if not pair_rows:
        return [
            "n_pairs: 0",
            "pearson_correlation_foldseek_vs_tm: NA",
        ]

    foldseek_scores = np.asarray([float(row["foldseek_score_max"]) for row in pair_rows], dtype=float)
    tm_scores = np.asarray([float(row["tm_score_max"]) for row in pair_rows], dtype=float)

    if foldseek_scores.size < 2:
        pearson = float("nan")
    else:
        pearson = float(np.corrcoef(foldseek_scores, tm_scores)[0, 1])

    return [
        f"n_pairs: {len(pair_rows)}",
        f"pearson_correlation_foldseek_vs_tm: {format_float(pearson) if math.isfinite(pearson) else 'NA'}",
        f"foldseek_score_max_mean: {format_float(np.mean(foldseek_scores))}",
        f"foldseek_score_max_median: {format_float(np.median(foldseek_scores))}",
        f"tm_score_max_mean: {format_float(np.mean(tm_scores))}",
        f"tm_score_max_median: {format_float(np.median(tm_scores))}",
    ]


def write_score_correlation_summary(path: Path, pair_rows: list[dict[str, str]]) -> None:
    path.write_text("\n".join(compute_score_correlation_summary(pair_rows)) + "\n", encoding="utf-8")


def write_run_log(
    path: Path,
    *,
    input_csv: Path,
    total_rows: int,
    valid_count: int,
    invalid_count: int,
    structures_dir: Path,
    foldseek_path: Path,
    tmalign_path: Path,
    db_prefix: Path,
    commands: list[list[str]],
    raw_hit_count: int,
    pair_count: int,
    runtime_seconds: float,
) -> None:
    lines = [
        f"input_csv: {input_csv}",
        f"total_rows: {total_rows}",
        f"valid_structures: {valid_count}",
        f"invalid_structures: {invalid_count}",
        f"structure_output_directory: {structures_dir}",
        f"foldseek_executable_path: {foldseek_path}",
        f"tmalign_executable_path: {tmalign_path}",
        f"foldseek_database_path: {db_prefix}",
        f"foldseek_primary_score_field: {PRIMARY_SCORE_FIELD}",
        "pair_aggregation_rule: foldseek_score_max=max(score_ij,score_ji), foldseek_score_min=min(score_ij,score_ji), tm_score_max=max(tm_norm_query,tm_norm_target), tm_score_min=min(tm_norm_query,tm_norm_target), primary_for_downstream=tm_score_max",
        f"total_raw_hits: {raw_hit_count}",
        f"total_unique_pairs: {pair_count}",
        f"runtime: {format_duration(runtime_seconds)}",
        "",
        "foldseek_commands:",
    ]
    lines.extend("  " + " ".join(command) for command in commands)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    start_time = time.time()

    if shutil.which(str(args.foldseek_path)) is None and not args.foldseek_path.exists():
        raise FileNotFoundError(f"Foldseek executable not found: {args.foldseek_path}")
    if shutil.which(str(args.tmalign_path)) is None and not args.tmalign_path.exists():
        raise FileNotFoundError(f"TM-align executable not found: {args.tmalign_path}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    structures_dir = args.output_dir / "structures"
    foldseek_db_dir = args.output_dir / "foldseek_db"
    tmp_dir = args.output_dir / "tmp"
    structures_dir.mkdir(parents=True, exist_ok=True)
    foldseek_db_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    valid_records, valid_rows, invalid_rows, total_rows = load_and_materialize_structures(
        input_csv=args.input_csv,
        structures_dir=structures_dir,
        progress_interval=args.progress_interval,
    )

    write_tsv(args.output_dir / "valid_structures.tsv", VALID_COLUMNS, valid_rows)
    write_tsv(args.output_dir / "invalid_structures.tsv", INVALID_COLUMNS, invalid_rows)

    if not valid_records:
        raise ValueError("no valid structures available for Foldseek processing")

    db_prefix = foldseek_db_dir / "structures_db"
    result_db = foldseek_db_dir / "search_result_db"
    raw_hits_tmp = foldseek_db_dir / "foldseek_raw_hits.tmp.tsv"

    commands = build_foldseek_commands(
        foldseek_path=args.foldseek_path,
        structures_dir=structures_dir,
        db_prefix=db_prefix,
        result_db=result_db,
        tmp_dir=tmp_dir,
        raw_hits_tmp=raw_hits_tmp,
        alignment_type=args.alignment_type,
        sensitivity=args.sensitivity,
        max_seqs=args.max_seqs,
    )

    command_logs: list[str] = []
    for stage_name, command in zip(["createdb", "search", "convertalis"], commands, strict=True):
        run_command(command, command_logs, stage_name=stage_name)

    raw_hit_rows = parse_foldseek_hits(raw_hits_tmp)
    write_tsv(args.output_dir / "foldseek_raw_hits.tsv", RAW_HIT_COLUMNS, raw_hit_rows)
    pair_rows = build_pairwise_rows(raw_hit_rows)
    pair_rows = run_tmalign_rescoring(
        pair_rows=pair_rows,
        structures_dir=structures_dir,
        tmalign_path=args.tmalign_path,
        progress_interval=args.progress_interval,
    )
    write_tsv(args.output_dir / "pairwise_foldseek.tsv", PAIR_COLUMNS, pair_rows)
    write_score_correlation_summary(args.output_dir / "score_correlation_summary.txt", pair_rows)

    runtime_seconds = time.time() - start_time
    write_run_log(
        path=args.output_dir / "run.log",
        input_csv=args.input_csv,
        total_rows=total_rows,
        valid_count=len(valid_rows),
        invalid_count=len(invalid_rows),
        structures_dir=structures_dir,
        foldseek_path=args.foldseek_path,
        tmalign_path=args.tmalign_path,
        db_prefix=db_prefix,
        commands=commands,
        raw_hit_count=len(raw_hit_rows),
        pair_count=len(pair_rows),
        runtime_seconds=runtime_seconds,
    )

    if raw_hits_tmp.exists():
        raw_hits_tmp.unlink()

    print(f"valid_structures={len(valid_rows)} invalid_structures={len(invalid_rows)}")
    print(f"raw_hits={len(raw_hit_rows)} unique_pairs={len(pair_rows)}")
    print(f"runtime={format_duration(runtime_seconds)}")


if __name__ == "__main__":
    main()
