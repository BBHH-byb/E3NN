from pathlib import Path

import numpy as np
import pandas as pd


BASE_DIR = Path("../../../domenico/27_NMR_screening/")
MASTER_LIST = BASE_DIR / "final_list_RMSD.txt"
PCA_DIR = BASE_DIR / "pca_results"
CENTROID_DIR = BASE_DIR / "pdb_centroids"
OUTPUT_DIR = Path("csv_input")
OUTPUT_COLUMNS = [
    "pdbid",
    "chain",
    "nresid",
    "resid_type",
    "resid_chain",
    "coordinate",
    "PC",
    "var",
]


def warn(message):
    print(f"Warning: {message}")


def format_number(value):
    return np.format_float_positional(float(value), trim="-")


def format_vector(values):
    return ",".join(format_number(value) for value in values)


def build_entry_prefix(pdbid, chain):
    return f"{pdbid}_{chain}"


def get_entry_paths(pdbid, chain):
    prefix = build_entry_prefix(pdbid, chain)
    return prefix, {
        "centroid": CENTROID_DIR / f"{prefix}_centroid.pdb",
        "pc": PCA_DIR / f"{prefix}_PCs.txt",
        "variance": PCA_DIR / f"{prefix}_pc_variances.txt",
    }


def get_missing_paths(paths):
    return [str(path) for path in paths.values() if not path.exists()]


def is_ca_atom_record(line):
    if not line.startswith(("ATOM", "HETATM")):
        return False

    atom_name = line[12:16].strip()
    alt_loc = line[16].strip()
    return atom_name == "CA" and alt_loc in {"", "A"}


def parse_ca_record(line):
    residue_name = line[17:20].strip()
    chain_id = line[21].strip()
    residue_seq = line[22:26].strip()
    insertion_code = line[26].strip()
    coordinates = [
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54]),
    ]
    residue_key = (chain_id, residue_seq, insertion_code)
    return residue_key, residue_name, chain_id, coordinates


def parse_centroid_pdb(pdb_path):
    residue_types = []
    residue_chains = []
    coordinates = []
    seen_residues = set()

    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not is_ca_atom_record(line):
                continue

            residue_key, residue_name, chain_id, atom_coordinates = parse_ca_record(line)
            if residue_key in seen_residues:
                continue

            seen_residues.add(residue_key)
            residue_types.append(residue_name)
            residue_chains.append(chain_id)
            coordinates.extend(atom_coordinates)

    return {
        "residue_types": residue_types,
        "residue_chains": residue_chains,
        "coordinates": coordinates,
    }


def load_pc_vector(pc_path):
    matrix = np.loadtxt(pc_path)
    array = np.asarray(matrix, dtype=float)

    if array.ndim == 0:
        return np.array([float(array)], dtype=float)
    if array.ndim == 1:
        return array[:1]
    return array[:, 0]


def load_variance(variance_path):
    values = np.loadtxt(variance_path)
    return float(np.ravel(np.asarray(values, dtype=float))[0])


def validate_entry_data(prefix, nresid, centroid_data, pc_vector):
    residue_count = len(centroid_data["residue_types"])
    coordinate_count = len(centroid_data["coordinates"])
    expected_vector_length = 3 * nresid

    if residue_count != nresid:
        warn(
            f"residue count mismatch for {prefix}: master list={nresid}, centroid={residue_count}"
        )
        return False

    if coordinate_count != expected_vector_length:
        warn(
            f"coordinate length mismatch for {prefix}: expected={expected_vector_length}, found={coordinate_count}"
        )
        return False

    if len(pc_vector) != expected_vector_length:
        warn(
            f"PC1 length mismatch for {prefix}: expected={expected_vector_length}, found={len(pc_vector)}"
        )
        return False

    return True


def make_dataset_row(pdbid, chain, nresid, centroid_data, pc_vector, variance):
    return {
        "pdbid": pdbid,
        "chain": chain,
        "nresid": nresid,
        "resid_type": ",".join(centroid_data["residue_types"]),
        "resid_chain": "".join(centroid_data["residue_chains"]),
        "coordinate": format_vector(centroid_data["coordinates"]),
        "PC": format_vector(pc_vector),
        "var": variance,
    }


def build_row(pdbid, chain, nresid):
    prefix, paths = get_entry_paths(pdbid, chain)
    missing_paths = get_missing_paths(paths)
    if missing_paths:
        warn(f"missing file for {prefix}: {', '.join(missing_paths)}")
        return None

    try:
        centroid_data = parse_centroid_pdb(paths["centroid"])
        pc_vector = load_pc_vector(paths["pc"])
        variance = load_variance(paths["variance"])
    except Exception as exc:
        warn(f"failed to parse {prefix}: {exc}")
        return None

    if not validate_entry_data(prefix, nresid, centroid_data, pc_vector):
        return None

    return make_dataset_row(pdbid, chain, nresid, centroid_data, pc_vector, variance)


def parse_master_list_line(raw_line, line_number):
    line = raw_line.strip()
    if not line:
        return None

    parts = line.split()
    if len(parts) < 4:
        warn(f"invalid master list line {line_number}: {raw_line.rstrip()}")
        return None

    pdbid, chain, _, length = parts[:4]
    try:
        nresid = int(length)
    except ValueError:
        warn(f"invalid residue count on line {line_number}: {length}")
        return None

    return pdbid, chain, nresid


def load_master_entries(master_list_path):
    entries = []

    with master_list_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            entry = parse_master_list_line(raw_line, line_number)
            if entry is not None:
                entries.append(entry)

    return entries


def build_dataset():
    if not MASTER_LIST.exists():
        warn(f"missing file: {MASTER_LIST}")
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    rows = []
    for pdbid, chain, nresid in load_master_entries(MASTER_LIST):
        row = build_row(pdbid, chain, nresid)
        if row is not None:
            rows.append(row)

    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def write_outputs(dataset):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    csv_path = OUTPUT_DIR / "input_dataset.csv"
    pkl_path = OUTPUT_DIR / "input_dataset.pkl"

    dataset.to_csv(csv_path, index=False)
    dataset.to_pickle(pkl_path)

    print(f"Saved {len(dataset)} rows to {csv_path}")
    print(f"Saved {len(dataset)} rows to {pkl_path}")


def main():
    dataset = build_dataset()
    write_outputs(dataset)


if __name__ == "__main__":
    main()
