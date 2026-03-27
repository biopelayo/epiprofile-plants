#!/usr/bin/env python3
"""
extract_empirical_rt.py
=======================
Extract empirical retention times for EpiProfile-PLANTS peptides from .ms1 files.

This script:
  1. Reads peptide definitions from plant_module_manifest.tsv
  2. Parses MATLAB module files to get exact m/z values at each charge state
  3. For each peptide, scans all .ms1 files to extract XICs (m/z +/- 10 ppm)
  4. Finds peak apex RT by smoothed XIC analysis
  5. Computes consensus RT across all samples (median, IQR)
  6. Outputs a TSV with empirical RTs ready for module update

Requirements:
  pip install numpy scipy pandas

Data format:
  The .ms1 files use pFind/xtract centroided format:
    H  header lines
    S  scan_number  scan_number
    I  RetTime  <seconds>
    I  NumberOfPeaks  <N>
    <m/z> <intensity>
    ...

Author: Pelayo Gonzalez de Lena (@pelamovic)
Date: 2026-03-26

Usage:
  python extract_empirical_rt.py                          # All peptides, all ms1 files
  python extract_empirical_rt.py --max-files 5            # Limit to 5 ms1 files
  python extract_empirical_rt.py --histone-type H2A       # Only H2A peptides
  python extract_empirical_rt.py --use-mat                # Use .mat files instead of .ms1
"""

import argparse
import csv
import os
import re
import struct
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Optional imports (checked at runtime)
try:
    from scipy.signal import savgol_filter
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import scipy.io as sio
    HAS_SCIPY_IO = True
except ImportError:
    HAS_SCIPY_IO = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = Path(r"D:/Antigravity/revision_tesis/epiprofile_plants_expanded")
MANIFEST_PATH = BASE_DIR / "metadata" / "plant_module_manifest.tsv"
OUTPUT_DIR = BASE_DIR / "metadata"

MS1_DIR = Path(r"E:/EpiProfile_Proyecto/EpiProfile_20_AT/RawData/MS1")

# MATLAB module locations
MODULE_DIRS = [
    Path(r"D:/Antigravity/repos/epiprofile_yuan_ecosistem/audit/generated_modules/AT"),
    Path(r"E:/EpiProfile_Proyecto/EpiProfile_20_AT/EpiProfile_PLANTS_resources/src"),
]

# ── Constants ────────────────────────────────────────────────────────────────
PPM_TOLERANCE = 10.0     # Mass tolerance in ppm (ZenoTOF centroided)
MIN_PEAK_INTENSITY = 50  # Minimum XIC peak intensity to consider
MIN_PEAK_WIDTH = 0.3     # Minimum peak width in minutes
MAX_PEAK_WIDTH = 5.0     # Maximum peak width in minutes
MIN_SAMPLES_DETECTED = 3 # Minimum samples where peptide must be detected
RT_OUTLIER_SD = 2.0      # SD threshold for RT outlier removal
SMOOTHING_WINDOW = 7     # Savitzky-Golay window (must be odd, >= 3)
SMOOTHING_ORDER = 2      # Savitzky-Golay polynomial order


# ── Amino acid masses (monoisotopic) ─────────────────────────────────────────
# From Getaamass.m in EpiProfile
AA_MASSES = {
    "G": 57.021464, "A": 71.037114, "V": 99.068414, "L": 113.084064,
    "I": 113.084064, "P": 97.052764, "F": 147.068414, "W": 186.079313,
    "M": 131.040485, "S": 87.032028, "T": 101.047679, "C": 103.009185,
    "Y": 163.063329, "H": 137.058912, "D": 115.026943, "E": 129.042593,
    "N": 114.042927, "Q": 128.058578, "K": 128.094963, "R": 156.101111,
}

# Water mass (lost during peptide bond formation)
H2O = 18.010565
# Proton mass
PROTON = 1.007276


# ── Modification masses (from GetMods.m) ─────────────────────────────────────
MOD_MASSES = {
    "pr":  56.026215,
    "ac":  42.010565,
    "me1": 70.041865,   # propionyl + methyl in EpiProfile workflow
    "me2": 28.031300,
    "me3": 42.046950,
    "ph":  79.966331,
    "cr":  68.026215,
    "ub":  114.042927,
    "ox":  15.994915,
    "pic": 56.026215,   # alias for propionylation
}


# ═══════════════════════════════════════════════════════════════════════════
# PEPTIDE M/Z CALCULATION
# ═══════════════════════════════════════════════════════════════════════════

def calculate_peptide_mz(sequence: str, mod_type_str: str, charge: int) -> float:
    """
    Calculate theoretical monoisotopic m/z for a peptide with modifications.

    Mirrors EpiProfile's calculate_pepmz logic:
    - Sum amino acid residue masses
    - Add H2O (water for free termini)
    - Add modification mass deltas
    - Divide by charge and add proton mass

    Args:
        sequence: Peptide sequence (e.g., "TKQTAR")
        mod_type_str: EpiProfile mod_type (e.g., "0,pr;2,pr;")
        charge: Charge state (1, 2, 3, etc.)

    Returns:
        Monoisotopic m/z value
    """
    # Sum residue masses
    neutral_mass = sum(AA_MASSES.get(aa, 0.0) for aa in sequence) + H2O

    # Add modification deltas
    if mod_type_str and mod_type_str.strip():
        pairs = [p.strip() for p in mod_type_str.split(";") if p.strip()]
        for pair in pairs:
            parts = pair.split(",")
            if len(parts) == 2:
                mod_code = parts[1].strip()
                delta = MOD_MASSES.get(mod_code, 0.0)
                neutral_mass += delta

    # Calculate m/z
    mz = (neutral_mass + charge * PROTON) / charge

    return mz


# ═══════════════════════════════════════════════════════════════════════════
# MS1 FILE PARSING
# ═══════════════════════════════════════════════════════════════════════════

def parse_ms1_file(filepath: Path) -> List[Dict]:
    """
    Parse a .ms1 file (pFind centroided format) into a list of scans.

    Each scan is a dict with:
      - scan_num: int
      - rt_seconds: float
      - rt_minutes: float
      - peaks: list of (mz, intensity) tuples

    Memory-efficient: yields scans one at a time via generator.
    For XIC extraction, we process scan-by-scan rather than loading all into memory.
    """
    scans = []
    current_scan = None

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("H\t"):
                continue  # Header line

            if line.startswith("S\t"):
                # Save previous scan
                if current_scan is not None:
                    scans.append(current_scan)

                parts = line.split("\t")
                scan_num = int(parts[1])
                current_scan = {
                    "scan_num": scan_num,
                    "rt_seconds": 0.0,
                    "rt_minutes": 0.0,
                    "peaks": [],
                }

            elif line.startswith("I\t"):
                if current_scan is not None:
                    parts = line.split("\t")
                    if len(parts) >= 3 and parts[1] == "RetTime":
                        rt_sec = float(parts[2])
                        current_scan["rt_seconds"] = rt_sec
                        current_scan["rt_minutes"] = rt_sec / 60.0

            else:
                # Peak line: m/z intensity
                if current_scan is not None:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            current_scan["peaks"].append((mz, intensity))
                        except ValueError:
                            pass

        # Don't forget the last scan
        if current_scan is not None:
            scans.append(current_scan)

    return scans


def extract_xic_from_scans(scans: List[Dict], target_mz: float,
                            ppm_tol: float = PPM_TOLERANCE) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract XIC (extracted ion chromatogram) for a target m/z from parsed scans.

    Args:
        scans: List of parsed scan dicts
        target_mz: Target m/z value
        ppm_tol: Mass tolerance in ppm

    Returns:
        (rt_array, intensity_array) in minutes
    """
    mz_tol = target_mz * ppm_tol * 1e-6

    rts = []
    intensities = []

    for scan in scans:
        rt = scan["rt_minutes"]
        # Sum intensities of all peaks within tolerance
        total_int = 0.0
        for mz, intensity in scan["peaks"]:
            if abs(mz - target_mz) <= mz_tol:
                total_int += intensity

        rts.append(rt)
        intensities.append(total_int)

    return np.array(rts), np.array(intensities)


def extract_xic_streaming(filepath: Path, target_mz: float,
                           ppm_tol: float = PPM_TOLERANCE) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract XIC from .ms1 file in streaming mode (lower memory than full parse).

    Reads the file line-by-line, only keeping RT and matching intensities.
    Preferred for large files.
    """
    mz_tol = target_mz * ppm_tol * 1e-6

    rts = []
    intensities = []
    current_rt = 0.0
    current_intensity = 0.0
    in_scan = False

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("S\t"):
                # Save previous scan's XIC point
                if in_scan:
                    rts.append(current_rt)
                    intensities.append(current_intensity)
                in_scan = True
                current_intensity = 0.0

            elif line.startswith("I\tRetTime\t"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    current_rt = float(parts[2]) / 60.0  # seconds -> minutes

            elif line.startswith("H\t") or line.startswith("I\t"):
                continue

            elif in_scan:
                # Peak line
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        if abs(mz - target_mz) <= mz_tol:
                            current_intensity += float(parts[1])
                    except ValueError:
                        pass

    # Last scan
    if in_scan:
        rts.append(current_rt)
        intensities.append(current_intensity)

    return np.array(rts), np.array(intensities)


# ═══════════════════════════════════════════════════════════════════════════
# MAT FILE SUPPORT (ALTERNATIVE)
# ═══════════════════════════════════════════════════════════════════════════

def extract_xic_from_mat(peaks_file: Path, scans_file: Path,
                          target_mz: float, ppm_tol: float = PPM_TOLERANCE
                          ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract XIC from MATLAB .mat files (MS1peaks + MS1scans format).

    EpiProfile preprocesses .ms1 into .mat files with:
      - MS1_peaks: sparse matrix of (m/z, intensity) per scan
      - MS1_index: [scan_num, RT, ...] per scan

    This is faster than reparsing .ms1 but requires scipy.io.
    """
    if not HAS_SCIPY_IO:
        raise ImportError("scipy.io required for .mat file reading. pip install scipy")

    # Load MAT files
    peaks_data = sio.loadmat(str(peaks_file))
    scans_data = sio.loadmat(str(scans_file))

    # MS1_index columns: [scan_number, retention_time, ...]
    # The exact structure depends on EpiProfile version
    ms1_index = None
    ms1_peaks = None

    for key in scans_data:
        if not key.startswith("_"):
            arr = scans_data[key]
            if hasattr(arr, "shape") and len(arr.shape) == 2 and arr.shape[1] >= 2:
                ms1_index = arr
                break

    for key in peaks_data:
        if not key.startswith("_"):
            arr = peaks_data[key]
            if hasattr(arr, "shape"):
                ms1_peaks = arr
                break

    if ms1_index is None or ms1_peaks is None:
        raise ValueError(f"Could not parse MAT file structure: {peaks_file}")

    mz_tol = target_mz * ppm_tol * 1e-6
    rts = []
    intensities = []

    n_scans = ms1_index.shape[0]
    for i in range(n_scans):
        rt = ms1_index[i, 1] / 60.0  # Assuming column 1 is RT in seconds

        # Extract peaks for this scan
        # Structure varies -- this is a best-effort implementation
        # The actual peaks matrix format depends on EpiProfile's Raw2MS output
        total_int = 0.0
        # ms1_peaks is typically a cell array or sparse matrix
        # Fallback: skip MAT parsing if structure is unclear
        # TODO: Adapt to actual EpiProfile MAT structure after inspection

        rts.append(rt)
        intensities.append(total_int)

    return np.array(rts), np.array(intensities)


# ═══════════════════════════════════════════════════════════════════════════
# PEAK DETECTION
# ═══════════════════════════════════════════════════════════════════════════

def smooth_xic(intensities: np.ndarray, window: int = SMOOTHING_WINDOW,
               order: int = SMOOTHING_ORDER) -> np.ndarray:
    """Apply Savitzky-Golay smoothing to XIC. Falls back to moving average."""
    if len(intensities) < window:
        return intensities

    if HAS_SCIPY:
        return savgol_filter(intensities, window, order)
    else:
        # Simple moving average fallback
        kernel = np.ones(min(window, len(intensities))) / min(window, len(intensities))
        return np.convolve(intensities, kernel, mode="same")


def find_xic_peaks(rts: np.ndarray, intensities: np.ndarray,
                    min_intensity: float = MIN_PEAK_INTENSITY,
                    min_width: float = MIN_PEAK_WIDTH,
                    max_width: float = MAX_PEAK_WIDTH,
                    ) -> List[Dict[str, float]]:
    """
    Find peaks in an XIC.

    Returns list of peak dicts with:
      - rt_apex: retention time at peak maximum
      - intensity: peak maximum intensity
      - rt_start: left boundary of peak
      - rt_end: right boundary of peak
      - width: peak width (FWHM estimate)
    """
    if len(rts) < 5:
        return []

    # Smooth
    smoothed = smooth_xic(intensities)

    # Find local maxima (points higher than both neighbors)
    peaks = []
    for i in range(1, len(smoothed) - 1):
        if smoothed[i] > smoothed[i-1] and smoothed[i] > smoothed[i+1]:
            if smoothed[i] >= min_intensity:
                # Estimate peak boundaries (descend to half-max on each side)
                half_max = smoothed[i] / 2.0

                # Left boundary
                left = i
                while left > 0 and smoothed[left] > half_max:
                    left -= 1
                rt_start = rts[left]

                # Right boundary
                right = i
                while right < len(smoothed) - 1 and smoothed[right] > half_max:
                    right += 1
                rt_end = rts[right]

                width = rt_end - rt_start

                if min_width <= width <= max_width:
                    peaks.append({
                        "rt_apex": rts[i],
                        "intensity": float(smoothed[i]),
                        "rt_start": rt_start,
                        "rt_end": rt_end,
                        "width": width,
                    })

    # Sort by intensity (highest first)
    peaks.sort(key=lambda p: p["intensity"], reverse=True)

    return peaks


# ═══════════════════════════════════════════════════════════════════════════
# MODULE PARSING
# ═══════════════════════════════════════════════════════════════════════════

def parse_manifest(path: Path) -> List[Dict[str, Any]]:
    """Parse plant_module_manifest.tsv."""
    modules = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
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


def find_module_file(file_name: str) -> Optional[Path]:
    """Locate a MATLAB module file."""
    for base_dir in MODULE_DIRS:
        candidate = base_dir / file_name
        if candidate.exists():
            return candidate
        for subdir in ["H2A", "H2B", "H1", "H3", "H4"]:
            candidate = base_dir / subdir / file_name
            if candidate.exists():
                return candidate
    return None


def parse_module_peptidoforms(filepath: Path) -> List[Dict]:
    """
    Parse MATLAB module to extract all peptidoforms with their m/z values.

    Returns list of dicts with: pep_seq, mod_short, mod_type, charges, mz_values, rt_ref
    """
    if not filepath.exists():
        return []

    content = filepath.read_text(encoding="utf-8", errors="replace")

    # Extract pep_seq
    seq_match = re.search(r"His\.pep_seq\s*=\s*'(\w+)'", content)
    if not seq_match:
        return []
    pep_seq = seq_match.group(1)

    # Extract mod_short
    short_match = re.search(r"His\.mod_short\s*=\s*\{([^}]+)\}", content)
    if not short_match:
        return []
    mod_shorts = re.findall(r"'([^']+)'", short_match.group(1))

    # Extract mod_type
    type_match = re.search(r"His\.mod_type\s*=\s*\{([^}]+)\}", content)
    if not type_match:
        return []
    mod_types = re.findall(r"'([^']+)'", type_match.group(1))

    # Extract rt_ref
    rt_match = re.search(r"His\.rt_ref\s*=\s*\[([^\]]+)\]", content)
    rt_values = []
    if rt_match:
        rt_values = [float(x.strip()) for x in rt_match.group(1).replace(";", " ").split()
                     if x.strip()]

    # Extract charge states: His.pep_ch = repmat([1 2 3], ...)
    ch_match = re.search(r"His\.pep_ch\s*=\s*repmat\(\[([^\]]+)\]", content)
    if ch_match:
        charges = [int(x) for x in ch_match.group(1).split()]
    else:
        charges = [1, 2, 3]  # default

    # Build peptidoforms
    n = min(len(mod_shorts), len(mod_types))
    peptidoforms = []
    for i in range(n):
        rt = rt_values[i] if i < len(rt_values) else 0.0
        # Calculate m/z for each charge state
        mz_values = {}
        for ch in charges:
            mz = calculate_peptide_mz(pep_seq, mod_types[i], ch)
            mz_values[ch] = mz

        peptidoforms.append({
            "pep_seq": pep_seq,
            "mod_short": mod_shorts[i],
            "mod_type": mod_types[i],
            "charges": charges,
            "mz_values": mz_values,
            "rt_ref": rt,
        })

    return peptidoforms


# ═══════════════════════════════════════════════════════════════════════════
# MAIN EXTRACTION PIPELINE
# ═══════════════════════════════════════════════════════════════════════════

def build_peptide_targets(modules: List[Dict],
                           histone_filter: Optional[str] = None,
                           only_missing_rt: bool = True,
                           ) -> List[Dict]:
    """
    Build list of (peptide, charge, target_mz) tuples for XIC extraction.

    Args:
        modules: Parsed manifest modules
        histone_filter: If set, only include this histone type (H2A, H2B, H1)
        only_missing_rt: If True, only include peptides with rt_ref = 0

    Returns:
        List of target dicts with peptide info and m/z values
    """
    targets = []
    seen = set()

    for mod in modules:
        if histone_filter and mod["histone_type"] != histone_filter:
            continue

        mod_file = find_module_file(mod["file_name"])
        if mod_file is None:
            # Use manifest info to build a minimal target
            pep_seq = mod["pep_seq"]
            mz_ch1 = mod["mz_unmod_ch1"]
            charges = [int(x) for x in mod["charges"].split()]

            key = f"{pep_seq}|unmod"
            if key not in seen:
                seen.add(key)
                # Estimate m/z for charge 2 from ch1 (most common for XIC)
                # mz_ch2 = (mz_ch1 * 1 - PROTON + 2*PROTON) / 2
                # Simplification: (mz_ch1 + PROTON) / 2
                mz_ch2 = (mz_ch1 + PROTON) / 2.0 if 2 in charges else None
                targets.append({
                    "file_name": mod["file_name"],
                    "pep_seq": pep_seq,
                    "mod_short": "unmod",
                    "mod_type": "0,pr;",
                    "histone_type": mod["histone_type"],
                    "target_charge": charges[0] if len(charges) == 1 else 2,
                    "target_mz": mz_ch2 if mz_ch2 else mz_ch1,
                    "all_mz": {1: mz_ch1},
                    "rt_ref": 0.0,
                })
            continue

        peptidoforms = parse_module_peptidoforms(mod_file)
        for pf in peptidoforms:
            if only_missing_rt and pf["rt_ref"] > 0:
                continue

            key = f"{pf['pep_seq']}|{pf['mod_type']}"
            if key in seen:
                continue
            seen.add(key)

            # Choose best charge for XIC (prefer charge 2 for typical tryptic peptides)
            preferred_charge = 2 if 2 in pf["charges"] else pf["charges"][0]

            targets.append({
                "file_name": mod["file_name"],
                "pep_seq": pf["pep_seq"],
                "mod_short": pf["mod_short"],
                "mod_type": pf["mod_type"],
                "histone_type": mod["histone_type"],
                "target_charge": preferred_charge,
                "target_mz": pf["mz_values"].get(preferred_charge, 0.0),
                "all_mz": pf["mz_values"],
                "rt_ref": pf["rt_ref"],
            })

    return targets


def extract_rts_from_ms1_files(targets: List[Dict], ms1_files: List[Path],
                                 verbose: bool = True) -> Dict[str, List[Dict]]:
    """
    Extract RTs for all target peptides across all ms1 files.

    Returns:
        Dict mapping peptide_key -> list of per-sample detection dicts
    """
    results = defaultdict(list)  # peptide_key -> [{sample, rt, intensity, ...}]

    n_files = len(ms1_files)
    n_targets = len(targets)
    print(f"\nExtracting XICs: {n_targets} targets x {n_files} files")
    print(f"{'=' * 60}")

    for fi, ms1_file in enumerate(ms1_files):
        sample_name = ms1_file.stem
        if verbose:
            print(f"\n  [{fi+1}/{n_files}] {sample_name}", end="", flush=True)

        t0 = time.time()

        # Process each target peptide
        for ti, target in enumerate(targets):
            key = f"{target['pep_seq']}|{target['mod_type']}"
            mz = target["target_mz"]

            if mz <= 0:
                continue

            # Extract XIC (streaming mode for memory efficiency)
            rts, intensities = extract_xic_streaming(ms1_file, mz, PPM_TOLERANCE)

            if len(rts) == 0:
                continue

            # Find peaks
            peaks = find_xic_peaks(rts, intensities)

            if peaks:
                # Take the highest-intensity peak
                best_peak = peaks[0]
                results[key].append({
                    "sample": sample_name,
                    "rt_apex": best_peak["rt_apex"],
                    "intensity": best_peak["intensity"],
                    "width": best_peak["width"],
                    "n_peaks_found": len(peaks),
                })

            if verbose and (ti + 1) % 20 == 0:
                print(".", end="", flush=True)

        elapsed = time.time() - t0
        if verbose:
            print(f"  ({elapsed:.1f}s)")

    return dict(results)


def compute_consensus_rts(results: Dict[str, List[Dict]],
                           targets: List[Dict]) -> List[Dict]:
    """
    Compute consensus RT for each peptide across samples.

    Applies outlier removal and computes statistics.
    """
    target_by_key = {}
    for t in targets:
        key = f"{t['pep_seq']}|{t['mod_type']}"
        target_by_key[key] = t

    consensus = []

    for key, detections in results.items():
        target = target_by_key.get(key, {})
        rts = np.array([d["rt_apex"] for d in detections])
        intensities = np.array([d["intensity"] for d in detections])

        n_detected = len(rts)

        if n_detected >= 3:
            # Remove outliers (> 2 SD from median)
            median_rt = np.median(rts)
            sd_rt = np.std(rts)
            if sd_rt > 0:
                mask = np.abs(rts - median_rt) <= RT_OUTLIER_SD * sd_rt
                rts_clean = rts[mask]
                intensities_clean = intensities[mask]
            else:
                rts_clean = rts
                intensities_clean = intensities

            consensus_rt = float(np.median(rts_clean))
            iqr = float(np.percentile(rts_clean, 75) - np.percentile(rts_clean, 25))
            mean_intensity = float(np.mean(intensities_clean))
        elif n_detected > 0:
            consensus_rt = float(np.median(rts))
            iqr = float(np.max(rts) - np.min(rts)) if n_detected > 1 else 0.0
            mean_intensity = float(np.mean(intensities))
        else:
            consensus_rt = 0.0
            iqr = 0.0
            mean_intensity = 0.0

        # Quality flag
        if n_detected >= MIN_SAMPLES_DETECTED and iqr < 1.0:
            quality = "HIGH"
        elif n_detected >= 2 and iqr < 2.0:
            quality = "MEDIUM"
        elif n_detected >= 1:
            quality = "LOW"
        else:
            quality = "NOT_DETECTED"

        consensus.append({
            "pep_seq": target.get("pep_seq", key.split("|")[0]),
            "mod_short": target.get("mod_short", ""),
            "mod_type": target.get("mod_type", key.split("|")[1] if "|" in key else ""),
            "histone_type": target.get("histone_type", ""),
            "file_name": target.get("file_name", ""),
            "target_mz": target.get("target_mz", 0.0),
            "target_charge": target.get("target_charge", 0),
            "rt_consensus": round(consensus_rt, 2),
            "rt_iqr": round(iqr, 2),
            "n_samples_detected": n_detected,
            "n_samples_total": 34,  # hardcoded for this dataset
            "mean_intensity": round(mean_intensity, 1),
            "quality": quality,
            "rt_ref_original": target.get("rt_ref", 0.0),
        })

    # Add targets that were never detected
    detected_keys = set(results.keys())
    for target in targets:
        key = f"{target['pep_seq']}|{target['mod_type']}"
        if key not in detected_keys:
            consensus.append({
                "pep_seq": target["pep_seq"],
                "mod_short": target["mod_short"],
                "mod_type": target["mod_type"],
                "histone_type": target["histone_type"],
                "file_name": target["file_name"],
                "target_mz": target["target_mz"],
                "target_charge": target["target_charge"],
                "rt_consensus": 0.0,
                "rt_iqr": 0.0,
                "n_samples_detected": 0,
                "n_samples_total": 34,
                "mean_intensity": 0.0,
                "quality": "NOT_DETECTED",
                "rt_ref_original": target.get("rt_ref", 0.0),
            })

    # Sort by histone type, then sequence
    consensus.sort(key=lambda x: (x["histone_type"], x["pep_seq"], x["mod_short"]))

    return consensus


# ═══════════════════════════════════════════════════════════════════════════
# OUTPUT
# ═══════════════════════════════════════════════════════════════════════════

def write_results_tsv(consensus: List[Dict], filepath: Path):
    """Write consensus RT results to TSV."""
    fieldnames = [
        "file_name", "histone_type", "pep_seq", "mod_short", "mod_type",
        "target_mz", "target_charge",
        "rt_consensus", "rt_iqr", "n_samples_detected", "n_samples_total",
        "mean_intensity", "quality", "rt_ref_original",
    ]

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in consensus:
            writer.writerow(row)


def write_matlab_update_snippet(consensus: List[Dict], filepath: Path):
    """
    Generate a MATLAB-compatible snippet for updating rt_ref values.

    Creates a lookup table that can be used by a MATLAB script to patch modules.
    """
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("%% RT Reference Updates from Empirical Extraction\n")
        f.write("%% Generated by extract_empirical_rt.py\n")
        f.write("%% Date: 2026-03-26\n")
        f.write("%% Quality: HIGH = detected in >=3 samples, IQR < 1 min\n")
        f.write("%%          MEDIUM = detected in >=2 samples, IQR < 2 min\n")
        f.write("%%          LOW = detected in 1 sample only\n")
        f.write("%%          NOT_DETECTED = not found in any sample\n\n")
        f.write("rt_updates = struct();\n\n")

        for row in consensus:
            if row["quality"] == "NOT_DETECTED":
                continue
            safe_key = row["pep_seq"] + "_" + row["mod_short"].replace("|", "_")
            safe_key = re.sub(r"[^a-zA-Z0-9_]", "_", safe_key)
            f.write(f"rt_updates.{safe_key}.seq = '{row['pep_seq']}';\n")
            f.write(f"rt_updates.{safe_key}.mod = '{row['mod_type']}';\n")
            f.write(f"rt_updates.{safe_key}.rt = {row['rt_consensus']:.2f};  "
                    f"% quality={row['quality']}, n={row['n_samples_detected']}, "
                    f"IQR={row['rt_iqr']:.2f}\n\n")


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Extract empirical RTs from .ms1 files for EpiProfile-PLANTS peptides"
    )
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH,
                        help="Path to plant_module_manifest.tsv")
    parser.add_argument("--ms1-dir", type=Path, default=MS1_DIR,
                        help="Directory containing .ms1 files")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR,
                        help="Output directory")
    parser.add_argument("--max-files", type=int, default=0,
                        help="Limit number of .ms1 files processed (0 = all)")
    parser.add_argument("--histone-type", type=str, default=None,
                        help="Filter by histone type (H2A, H2B, H1)")
    parser.add_argument("--include-known-rt", action="store_true",
                        help="Also extract for peptides with known rt_ref (for validation)")
    parser.add_argument("--use-mat", action="store_true",
                        help="Use .mat files instead of .ms1 (requires scipy)")
    parser.add_argument("--ppm", type=float, default=PPM_TOLERANCE,
                        help=f"Mass tolerance in ppm (default: {PPM_TOLERANCE})")
    parser.add_argument("--verbose", action="store_true", default=True,
                        help="Print progress information")
    args = parser.parse_args()

    global PPM_TOLERANCE
    PPM_TOLERANCE = args.ppm

    print("=" * 70)
    print("Empirical RT Extraction for EpiProfile-PLANTS")
    print("=" * 70)

    # --- 1. Read manifest ---
    print(f"\n[1/5] Reading manifest: {args.manifest}")
    modules = parse_manifest(args.manifest)
    print(f"  Found {len(modules)} AT modules")

    # --- 2. Build targets ---
    print(f"\n[2/5] Building peptide targets...")
    targets = build_peptide_targets(
        modules,
        histone_filter=args.histone_type,
        only_missing_rt=not args.include_known_rt,
    )
    print(f"  Targets to extract: {len(targets)}")

    if not targets:
        print("  No targets found. Check filters.")
        return

    # Print target summary
    by_type = defaultdict(int)
    for t in targets:
        by_type[t["histone_type"]] += 1
    for ht, count in sorted(by_type.items()):
        print(f"    {ht}: {count}")

    # --- 3. Find ms1 files ---
    print(f"\n[3/5] Finding .ms1 files in: {args.ms1_dir}")
    if args.use_mat:
        mat_files = sorted(args.ms1_dir.glob("*_MS1peaks.mat"))
        print(f"  Found {len(mat_files)} MAT file pairs")
        if args.max_files > 0:
            mat_files = mat_files[:args.max_files]
            print(f"  Limited to {args.max_files} files")
        ms1_files = mat_files  # Will need special handling
    else:
        ms1_files = sorted(args.ms1_dir.glob("*.ms1"))
        print(f"  Found {len(ms1_files)} .ms1 files")
        if args.max_files > 0:
            ms1_files = ms1_files[:args.max_files]
            print(f"  Limited to {args.max_files} files")

    if not ms1_files:
        print("  No data files found. Check path.")
        return

    # --- 4. Extract XICs ---
    print(f"\n[4/5] Extracting XICs...")
    print(f"  Mass tolerance: {PPM_TOLERANCE} ppm")
    print(f"  Min peak intensity: {MIN_PEAK_INTENSITY}")
    print(f"  Smoothing: Savitzky-Golay (window={SMOOTHING_WINDOW}, order={SMOOTHING_ORDER})")

    t_start = time.time()
    results = extract_rts_from_ms1_files(targets, ms1_files, verbose=args.verbose)
    t_elapsed = time.time() - t_start

    n_detected = sum(1 for v in results.values() if len(v) > 0)
    print(f"\n  Extraction complete in {t_elapsed:.1f}s")
    print(f"  Peptides detected: {n_detected}/{len(targets)}")

    # --- 5. Compute consensus and save ---
    print(f"\n[5/5] Computing consensus RTs and saving results...")
    consensus = compute_consensus_rts(results, targets)

    # Output files
    tsv_path = args.output_dir / "rt_empirical_extraction.tsv"
    write_results_tsv(consensus, tsv_path)
    print(f"  Results TSV: {tsv_path}")

    matlab_path = args.output_dir / "rt_empirical_matlab_updates.m"
    write_matlab_update_snippet(consensus, matlab_path)
    print(f"  MATLAB snippet: {matlab_path}")

    # --- Summary ---
    quality_counts = defaultdict(int)
    for row in consensus:
        quality_counts[row["quality"]] += 1

    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Total targets:    {len(targets)}")
    print(f"  HIGH quality:     {quality_counts.get('HIGH', 0)}")
    print(f"  MEDIUM quality:   {quality_counts.get('MEDIUM', 0)}")
    print(f"  LOW quality:      {quality_counts.get('LOW', 0)}")
    print(f"  NOT DETECTED:     {quality_counts.get('NOT_DETECTED', 0)}")
    print(f"\n  Output files:")
    print(f"    {tsv_path}")
    print(f"    {matlab_path}")
    print()

    # Validation hint
    if not args.include_known_rt:
        print("  TIP: Run with --include-known-rt to also extract H3/H4 peptides")
        print("       and validate against their known rt_ref values.")

    print("=" * 70)


if __name__ == "__main__":
    main()
