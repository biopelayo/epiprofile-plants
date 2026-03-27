#!/usr/bin/env python3
"""
predict_rt_deeplc.py
====================
Predict retention times for EpiProfile-PLANTS H2A/H2B/H1 peptides using DeepLC.

This script:
  1. Reads the plant_module_manifest.tsv to get peptide sequences and modifications
  2. Parses the corresponding MATLAB modules to extract mod_type strings
  3. Converts EpiProfile modification notation to DeepLC format
  4. Optionally uses known H3/H4 RTs as calibration anchors
  5. Runs DeepLC prediction (or generates the input files for manual execution)
  6. Outputs a TSV with predicted RTs for all peptides

Requirements:
  pip install deeplc pandas numpy

If DeepLC is not installed, the script generates input files and prints the
command to run DeepLC manually.

Author: Pelayo Gonzalez de Lena (@pelamovic)
Date: 2026-03-26

Usage:
  python predict_rt_deeplc.py                    # Generate inputs only
  python predict_rt_deeplc.py --run               # Run DeepLC prediction
  python predict_rt_deeplc.py --run --calibrate   # Run with H3/H4 calibration
"""

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = Path(r"D:/Antigravity/revision_tesis/epiprofile_plants_expanded")
MANIFEST_PATH = BASE_DIR / "metadata" / "plant_module_manifest.tsv"
RT_COMPARISON_PATH = BASE_DIR / "metadata" / "RT_COMPARISON.tsv"
OUTPUT_DIR = BASE_DIR / "metadata"
SCRIPTS_DIR = BASE_DIR / "scripts"

# Module directories (for reading mod_type from MATLAB files)
MODULE_DIRS = [
    Path(r"D:/Antigravity/repos/epiprofile_yuan_ecosistem/audit/generated_modules/AT"),
    Path(r"E:/EpiProfile_Proyecto/EpiProfile_20_AT/EpiProfile_PLANTS_resources/src"),
]

# ── Modification Mass Deltas (from GetMods.m) ───────────────────────────────
# These are monoisotopic mass deltas used in EpiProfile's propionylation workflow.
# IMPORTANT: In EpiProfile, me1 on K includes propionyl (70.041865 = 56.026215 + 14.01565)
MOD_MASSES = {
    "pr":  56.026215,   # Propionylation (derivatization)
    "ac":  42.010565,   # Acetylation
    "me1": 70.041865,   # Monomethyl (propionyl + methyl in EpiProfile workflow)
    "me2": 28.031300,   # Dimethylation
    "me3": 42.046950,   # Trimethylation
    "ph":  79.966331,   # Phosphorylation
    "cr":  68.026215,   # Crotonylation
    "ub":  114.042927,  # Ubiquitylation (diglycine remnant)
    "ox":  15.994915,   # Oxidation
    "pic": 56.026215,   # Alias for propionylation used in some modules
}

# ── DeepLC Modification Mapping ─────────────────────────────────────────────
# Map EpiProfile mod codes to DeepLC/Unimod names.
# DeepLC uses Unimod nomenclature for modifications.
DEEPLC_MOD_MAP = {
    "pr":  "Propionyl",
    "pic": "Propionyl",
    "ac":  "Acetyl",
    "me1": "Methyl",        # NOTE: DeepLC uses raw methyl, not propionyl+methyl
    "me2": "Dimethyl",
    "me3": "Trimethyl",
    "ph":  "Phospho",
    "cr":  "Crotonyl",
    "ox":  "Oxidation",
}

# ── Known H3/H4 Calibration Anchors ────────────────────────────────────────
# From H3/H4 MATLAB modules with real rt_ref values.
# Format: (sequence, mod_type_string, observed_RT_minutes)
# These are for the CoA-single propionylation chemistry on ZenoTOF / Evosep 30SPD.
CALIBRATION_ANCHORS = [
    # H3 3-8: TKQTAR
    ("TKQTAR", "0,pr;2,pr;", 18.34),          # unmod
    ("TKQTAR", "0,pr;2,me1;", 21.58),         # K4me1
    ("TKQTAR", "0,pr;2,me2;", 10.82),         # K4me2
    ("TKQTAR", "0,pr;2,me3;", 10.80),         # K4me3
    ("TKQTAR", "0,pr;2,ac;", 16.64),          # K4ac
    # H3 18-26: KQLATKAAR
    ("KQLATKAAR", "0,pr;1,pr;6,pr;", 37.41),  # unmod
    # H3 27-40: KSAPATGGVKKPHR (AT motif)
    ("KSAPATGGVKKPHR", "0,pr;1,pr;10,pr;11,pr;", 26.08),  # unmod
    # H3 117-128: VTIMPKDIQLAR
    ("VTIMPKDIQLAR", "0,pr;6,pr;", 48.69),    # unmod
]


# ═══════════════════════════════════════════════════════════════════════════
# PARSING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def parse_manifest(path: Path) -> List[Dict[str, Any]]:
    """Parse plant_module_manifest.tsv into list of module dicts."""
    modules = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Only AT species modules
            if row.get("species", "") != "AT":
                continue
            modules.append({
                "file_name": row["file_name"],
                "histone_name": row["histone_name"],
                "histone_type": row["histone_type"],
                "pep_seq": row["pep_seq"],
                "pep_len": int(row["pep_len"]),
                "pos_start": int(row["pos_start"]),
                "pos_end": int(row["pos_end"]),
                "n_mods": int(row["n_mods"]),
                "mod_shorts": row["mod_shorts"],
                "charges": row["charges"],
                "mz_unmod_ch1": float(row["mz_unmod_ch1"]),
                "variant_code": row["variant_code"],
                "path": row["path"],
            })
    return modules


def parse_mod_type_from_matlab(filepath: Path) -> List[Tuple[str, str, float]]:
    """
    Parse a MATLAB module file to extract (mod_short, mod_type, rt_ref) tuples.

    Returns list of (mod_label, mod_type_string, rt_ref_value) for each peptidoform
    in the module's init_histone() function.
    """
    results = []
    if not filepath.exists():
        return results

    content = filepath.read_text(encoding="utf-8", errors="replace")

    # Extract pep_seq
    seq_match = re.search(r"His\.pep_seq\s*=\s*'(\w+)'", content)
    if not seq_match:
        return results
    pep_seq = seq_match.group(1)

    # Extract mod_short (cell array)
    short_match = re.search(r"His\.mod_short\s*=\s*\{([^}]+)\}", content)
    if not short_match:
        return results
    mod_shorts = re.findall(r"'([^']+)'", short_match.group(1))

    # Extract mod_type (cell array)
    type_match = re.search(r"His\.mod_type\s*=\s*\{([^}]+)\}", content)
    if not type_match:
        return results
    mod_types = re.findall(r"'([^']+)'", type_match.group(1))

    # Extract rt_ref (numeric array)
    rt_match = re.search(r"His\.rt_ref\s*=\s*\[([^\]]+)\]", content)
    if not rt_match:
        return results
    rt_values = [float(x.strip()) for x in rt_match.group(1).replace(";", " ").split()
                 if x.strip()]

    # Align arrays
    n = min(len(mod_shorts), len(mod_types), len(rt_values))
    for i in range(n):
        results.append((mod_shorts[i], mod_types[i], rt_values[i]))

    return results


def find_module_file(file_name: str) -> Optional[Path]:
    """Locate a MATLAB module file in known directories."""
    for base_dir in MODULE_DIRS:
        # Direct search
        candidate = base_dir / file_name
        if candidate.exists():
            return candidate
        # Search in subdirectories (H2A, H2B, H1, H3, H4)
        for subdir in ["H2A", "H2B", "H1", "H3", "H4"]:
            candidate = base_dir / subdir / file_name
            if candidate.exists():
                return candidate
    return None


# ═══════════════════════════════════════════════════════════════════════════
# DEEPLC FORMAT CONVERSION
# ═══════════════════════════════════════════════════════════════════════════

def epiprofile_mod_to_deeplc(pep_seq: str, mod_type_str: str) -> str:
    """
    Convert EpiProfile mod_type string to DeepLC modification format.

    EpiProfile format: "0,pr;2,me1;6,pr;" (0-indexed, 0=N-term)
    DeepLC format: "1|Propionyl|3|Methyl|7|Propionyl" (1-indexed, position in sequence)

    Note on propionylation handling:
    - EpiProfile treats propionylation as a derivatization applied to all unmodified K
      and the N-terminus. In DeepLC, we encode these as "Propionyl" modifications.
    - For me1 in EpiProfile's propionyl workflow, the mass delta (70.04) includes
      propionyl + methyl. In DeepLC, we encode this as just "Methyl" since DeepLC's
      model was NOT trained on propionylated peptides. Instead, we add a separate
      "Propionyl" on the same position if needed. However, DeepLC may not handle
      stacked modifications on the same residue.
    - PRAGMATIC APPROACH: We encode each EpiProfile modification as-is using the
      DeepLC mod map, and accept that predictions for propionylated peptides will
      have higher uncertainty.
    """
    if not mod_type_str or mod_type_str.strip() == "":
        return ""

    deeplc_mods = []

    # Parse EpiProfile mod_type: "pos,type;" pairs
    pairs = [p.strip() for p in mod_type_str.split(";") if p.strip()]
    for pair in pairs:
        parts = pair.split(",")
        if len(parts) != 2:
            continue
        pos = int(parts[0])
        mod_code = parts[1].strip()

        # Map to DeepLC name
        deeplc_name = DEEPLC_MOD_MAP.get(mod_code)
        if deeplc_name is None:
            print(f"  WARNING: Unknown mod '{mod_code}' at pos {pos}, skipping",
                  file=sys.stderr)
            continue

        # DeepLC uses 1-indexed positions.
        # EpiProfile pos 0 = N-terminal modification.
        # In DeepLC, N-terminal mods go on position 0 or the first residue (1).
        if pos == 0:
            deeplc_pos = 0  # N-terminal
        else:
            deeplc_pos = pos  # 1-indexed in EpiProfile = 1-indexed in sequence

        deeplc_mods.append(f"{deeplc_pos}|{deeplc_name}")

    return "|".join(deeplc_mods)


def build_deeplc_input(modules: List[Dict], include_calibration: bool = True) -> Tuple[list, list]:
    """
    Build DeepLC input tables.

    Returns:
        prediction_rows: List of dicts for peptides needing RT prediction
        calibration_rows: List of dicts for calibration peptides (with known RT)
    """
    prediction_rows = []
    calibration_rows = []
    seen_peptides = set()  # avoid duplicates (same seq + same mods)

    # --- Process modules from manifest ---
    for mod in modules:
        file_name = mod["file_name"]
        pep_seq = mod["pep_seq"]

        # Find and parse the MATLAB module
        mod_file = find_module_file(file_name)
        if mod_file is None:
            # Fallback: create a minimal mod_type for unmodified peptide
            # with N-terminal propionylation
            mod_type_str = "0,pr;"
            mod_short = "unmod"
            rt_ref = 0.0
            peptidoforms = [(mod_short, mod_type_str, rt_ref)]
        else:
            peptidoforms = parse_mod_type_from_matlab(mod_file)
            if not peptidoforms:
                mod_type_str = "0,pr;"
                mod_short = "unmod"
                rt_ref = 0.0
                peptidoforms = [(mod_short, mod_type_str, rt_ref)]

        for mod_short, mod_type_str, rt_ref in peptidoforms:
            key = f"{pep_seq}|{mod_type_str}"
            if key in seen_peptides:
                continue
            seen_peptides.add(key)

            deeplc_mods = epiprofile_mod_to_deeplc(pep_seq, mod_type_str)

            row = {
                "seq": pep_seq,
                "modifications": deeplc_mods,
                "mod_type_epiprofile": mod_type_str,
                "mod_short": mod_short,
                "histone_type": mod["histone_type"],
                "file_name": file_name,
                "rt_ref_original": rt_ref,
            }

            if rt_ref > 0:
                row["tr"] = rt_ref  # known RT -> can serve as calibration
            else:
                row["tr"] = ""  # needs prediction

            prediction_rows.append(row)

    # --- Add calibration anchors ---
    if include_calibration:
        for seq, mod_type, rt in CALIBRATION_ANCHORS:
            key = f"{seq}|{mod_type}"
            if key in seen_peptides:
                # Update existing entry with known RT
                for row in prediction_rows:
                    if row["seq"] == seq and row["mod_type_epiprofile"] == mod_type:
                        row["tr"] = rt
                continue

            deeplc_mods = epiprofile_mod_to_deeplc(seq, mod_type)
            calibration_rows.append({
                "seq": seq,
                "modifications": deeplc_mods,
                "mod_type_epiprofile": mod_type,
                "mod_short": "calibration",
                "histone_type": "H3/H4",
                "file_name": "calibration_anchor",
                "rt_ref_original": rt,
                "tr": rt,
            })

    return prediction_rows, calibration_rows


# ═══════════════════════════════════════════════════════════════════════════
# DEEPLC EXECUTION
# ═══════════════════════════════════════════════════════════════════════════

def run_deeplc(prediction_file: Path, calibration_file: Optional[Path] = None,
               output_file: Optional[Path] = None) -> Optional[Path]:
    """
    Run DeepLC prediction. Returns path to results file.

    If DeepLC is not installed, prints instructions and returns None.
    """
    try:
        from deeplc import DeepLC
    except ImportError:
        print("\n" + "=" * 70)
        print("DeepLC is not installed in this environment.")
        print("To install: pip install deeplc")
        print("\nAlternatively, run DeepLC manually:")
        print(f"  deeplc --file_pred {prediction_file}", end="")
        if calibration_file:
            print(f" \\\n    --file_cal {calibration_file}", end="")
        if output_file:
            print(f" \\\n    --file_pred_out {output_file}", end="")
        print("\n" + "=" * 70)
        return None

    import pandas as pd

    # Load prediction peptides
    df_pred = pd.read_csv(prediction_file, sep=",")

    # Load calibration peptides if provided
    df_cal = None
    if calibration_file and calibration_file.exists():
        df_cal = pd.read_csv(calibration_file, sep=",")
        df_cal = df_cal.rename(columns={"tr": "tr"})

    # Initialize DeepLC
    dlc = DeepLC(
        split_cal=10,    # number of calibration splits
        n_jobs=4,        # parallel jobs
        verbose=True,
    )

    # Calibrate if calibration data is provided
    if df_cal is not None and len(df_cal) > 0:
        print(f"Calibrating DeepLC with {len(df_cal)} anchor peptides...")
        dlc.calibrate_preds(seq_df=df_cal)
        print("Calibration complete.")

    # Predict
    print(f"Predicting RTs for {len(df_pred)} peptides...")
    preds = dlc.make_preds(seq_df=df_pred)

    # Add predictions to dataframe
    df_pred["rt_predicted_deeplc"] = preds

    # Save results
    if output_file is None:
        output_file = OUTPUT_DIR / "rt_predictions_deeplc.tsv"
    df_pred.to_csv(output_file, sep="\t", index=False)
    print(f"Predictions saved to: {output_file}")

    return output_file


# ═══════════════════════════════════════════════════════════════════════════
# iRT CONVERSION
# ═══════════════════════════════════════════════════════════════════════════

def compute_irt_coefficients(anchor_rts_instrument: Dict[str, float],
                             anchor_irt_reference: Dict[str, float]) -> Tuple[float, float]:
    """
    Compute iRT linear transformation coefficients.

    iRT = a * RT_instrument + b

    Args:
        anchor_rts_instrument: {peptide_key: RT_on_this_instrument}
        anchor_irt_reference: {peptide_key: iRT_reference_value}

    Returns:
        (a, b) coefficients for the linear transformation
    """
    common_keys = set(anchor_rts_instrument.keys()) & set(anchor_irt_reference.keys())
    if len(common_keys) < 3:
        raise ValueError(f"Need >= 3 common anchors, got {len(common_keys)}")

    x = np.array([anchor_rts_instrument[k] for k in sorted(common_keys)])
    y = np.array([anchor_irt_reference[k] for k in sorted(common_keys)])

    # Linear regression: iRT = a * RT + b
    A = np.vstack([x, np.ones(len(x))]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]

    return float(a), float(b)


def rt_to_irt(rt: float, a: float, b: float) -> float:
    """Convert absolute RT to iRT using linear coefficients."""
    return a * rt + b


def irt_to_rt(irt: float, a: float, b: float) -> float:
    """Convert iRT back to absolute RT for a specific instrument."""
    return (irt - b) / a


# ═══════════════════════════════════════════════════════════════════════════
# OUTPUT GENERATION
# ═══════════════════════════════════════════════════════════════════════════

def write_deeplc_input_csv(rows: list, filepath: Path, include_tr: bool = False):
    """Write DeepLC-compatible CSV input file."""
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if include_tr:
            writer.writerow(["seq", "modifications", "tr"])
            for row in rows:
                writer.writerow([row["seq"], row["modifications"], row.get("tr", "")])
        else:
            writer.writerow(["seq", "modifications"])
            for row in rows:
                writer.writerow([row["seq"], row["modifications"]])


def write_full_output(prediction_rows: list, calibration_rows: list,
                      filepath: Path):
    """Write complete output TSV with all metadata."""
    fieldnames = [
        "file_name", "histone_type", "seq", "mod_short",
        "mod_type_epiprofile", "modifications_deeplc",
        "rt_ref_original", "needs_prediction",
    ]

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        all_rows = prediction_rows + calibration_rows
        for row in all_rows:
            writer.writerow({
                "file_name": row["file_name"],
                "histone_type": row["histone_type"],
                "seq": row["seq"],
                "mod_short": row["mod_short"],
                "mod_type_epiprofile": row["mod_type_epiprofile"],
                "modifications_deeplc": row["modifications"],
                "rt_ref_original": row["rt_ref_original"],
                "needs_prediction": "YES" if row.get("tr", "") == "" else "NO",
            })


def write_deeplc_custom_mods(filepath: Path):
    """Write the DeepLC custom modifications JSON file."""
    custom_mods = {
        "Propionyl": {
            "unimod_name": "Propionyl",
            "unimod_accession": 58,
            "mass_shift": 56.026215,
            "composition": {"H": 4, "C": 3, "O": 1},
        },
        "Acetyl": {
            "unimod_name": "Acetyl",
            "unimod_accession": 1,
            "mass_shift": 42.010565,
            "composition": {"H": 2, "C": 2, "O": 1},
        },
        "Methyl": {
            "unimod_name": "Methyl",
            "unimod_accession": 34,
            "mass_shift": 14.01565,
            "composition": {"H": 2, "C": 1},
        },
        "Dimethyl": {
            "unimod_name": "Dimethyl",
            "unimod_accession": 36,
            "mass_shift": 28.031300,
            "composition": {"H": 4, "C": 2},
        },
        "Trimethyl": {
            "unimod_name": "Trimethyl",
            "unimod_accession": 37,
            "mass_shift": 42.046950,
            "composition": {"H": 6, "C": 3},
        },
        "Phospho": {
            "unimod_name": "Phospho",
            "unimod_accession": 21,
            "mass_shift": 79.966331,
            "composition": {"H": 1, "O": 3, "P": 1},
        },
        "Crotonyl": {
            "unimod_name": "Crotonyl",
            "unimod_accession": 1363,
            "mass_shift": 68.026215,
            "composition": {"H": 4, "C": 4, "O": 1},
        },
        "Oxidation": {
            "unimod_name": "Oxidation",
            "unimod_accession": 35,
            "mass_shift": 15.994915,
            "composition": {"O": 1},
        },
    }

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(custom_mods, f, indent=2)


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Predict retention times for EpiProfile-PLANTS peptides using DeepLC"
    )
    parser.add_argument("--run", action="store_true",
                        help="Actually run DeepLC (requires deeplc package)")
    parser.add_argument("--calibrate", action="store_true",
                        help="Use H3/H4 calibration anchors for DeepLC")
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH,
                        help="Path to plant_module_manifest.tsv")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR,
                        help="Output directory for generated files")
    args = parser.parse_args()

    print("=" * 70)
    print("DeepLC RT Prediction for EpiProfile-PLANTS")
    print("=" * 70)

    # --- 1. Read manifest ---
    print(f"\n[1/5] Reading manifest: {args.manifest}")
    modules = parse_manifest(args.manifest)
    print(f"  Found {len(modules)} AT modules")

    # Count by histone type
    type_counts = {}
    for m in modules:
        ht = m["histone_type"]
        type_counts[ht] = type_counts.get(ht, 0) + 1
    for ht, count in sorted(type_counts.items()):
        print(f"    {ht}: {count}")

    # --- 2. Build DeepLC input ---
    print(f"\n[2/5] Building DeepLC input tables...")
    prediction_rows, calibration_rows = build_deeplc_input(
        modules, include_calibration=args.calibrate
    )

    n_need_pred = sum(1 for r in prediction_rows if r.get("tr", "") == "")
    n_have_rt = sum(1 for r in prediction_rows if r.get("tr", "") != "")
    print(f"  Peptidoforms needing prediction: {n_need_pred}")
    print(f"  Peptidoforms with known RT: {n_have_rt}")
    print(f"  Calibration anchors: {len(calibration_rows)}")

    # --- 3. Write output files ---
    print(f"\n[3/5] Writing output files to: {args.output_dir}")

    # DeepLC prediction input
    pred_input_file = args.output_dir / "deeplc_input_prediction.csv"
    pred_rows_no_rt = [r for r in prediction_rows if r.get("tr", "") == ""]
    write_deeplc_input_csv(pred_rows_no_rt, pred_input_file, include_tr=False)
    print(f"  Prediction input: {pred_input_file} ({len(pred_rows_no_rt)} peptides)")

    # DeepLC calibration input
    if args.calibrate:
        cal_input_file = args.output_dir / "deeplc_input_calibration.csv"
        cal_rows = [r for r in prediction_rows if r.get("tr", "") != ""]
        cal_rows.extend(calibration_rows)
        write_deeplc_input_csv(cal_rows, cal_input_file, include_tr=True)
        print(f"  Calibration input: {cal_input_file} ({len(cal_rows)} peptides)")
    else:
        cal_input_file = None

    # Custom modifications JSON
    mods_file = args.output_dir / "deeplc_custom_mods.json"
    write_deeplc_custom_mods(mods_file)
    print(f"  Custom mods file: {mods_file}")

    # Full metadata table
    full_output_file = args.output_dir / "rt_prediction_input_table.tsv"
    write_full_output(prediction_rows, calibration_rows, full_output_file)
    print(f"  Full metadata table: {full_output_file}")

    # --- 4. Run DeepLC (optional) ---
    if args.run:
        print(f"\n[4/5] Running DeepLC prediction...")
        result_file = run_deeplc(
            prediction_file=pred_input_file,
            calibration_file=cal_input_file,
            output_file=args.output_dir / "rt_predictions_deeplc.tsv",
        )
        if result_file:
            print(f"  Results: {result_file}")
    else:
        print(f"\n[4/5] Skipping DeepLC execution (use --run to enable)")
        print(f"\n  To run manually:")
        print(f"    pip install deeplc")
        print(f"    python predict_rt_deeplc.py --run --calibrate")

    # --- 5. Summary ---
    print(f"\n[5/5] Summary")
    print(f"  Generated files:")
    print(f"    - {pred_input_file}")
    if args.calibrate:
        print(f"    - {cal_input_file}")
    print(f"    - {mods_file}")
    print(f"    - {full_output_file}")
    print()
    print("  IMPORTANT CAVEATS:")
    print("  1. DeepLC was NOT trained on propionylated histone peptides.")
    print("     Predictions for derivatized peptides will have higher uncertainty.")
    print("  2. The me1 mass in EpiProfile (70.04 Da) includes propionyl + methyl.")
    print("     DeepLC sees 'Methyl' only (14.02 Da). RT prediction may differ.")
    print("  3. Always validate predictions against empirical data.")
    print("  4. Recommended: run extract_empirical_rt.py FIRST, use DeepLC only")
    print("     for peptides not detected experimentally.")
    print()
    print("=" * 70)


if __name__ == "__main__":
    main()
