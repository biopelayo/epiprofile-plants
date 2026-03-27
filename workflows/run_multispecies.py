#!/usr/bin/env python3
"""
EpiProfile-PLANTS multi-species launcher.
Runs AT, CR, MP bundles on shared MS data with separate output dirs.
Supports parallel (if MATLAB license allows) or sequential fallback.

Usage:
    python run_multispecies.py
    python run_multispecies.py --species AT CR
    python run_multispecies.py --sequential
    python run_multispecies.py --raw-path "C:\\path\\to\\MS1_MS2"
"""
import argparse
import subprocess
import os
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# === CONFIGURATION ===
MATLAB_EXE = r"C:\Program Files\MATLAB\R2023a\bin\matlab.exe"
SRC_ROOT = Path(r"D:\Antigravity\repos\epiprofile-plants\SRC_update")
RAW_PATH = r"C:\Users\geope\Desktop\Tesis_Epigenomica\EpiProfile_PLANTS\test_chapter2\MS1_MS2"
RUNS_ROOT = SRC_ROOT / "RUNS"
VALID_SPECIES = ["AT", "CR", "MP"]


def generate_paras(species: str, raw_path: str, output_path: str) -> str:
    """Generate paras.txt for a species run. Returns path to generated file."""
    run_dir = Path(output_path)
    run_dir.mkdir(parents=True, exist_ok=True)

    paras_path = run_dir / f"paras_{species}.txt"
    content = (
        f"[EpiProfile]\n"
        f"raw_path={raw_path}\n"
        f"norganism=1\n"
        f"nsource=1\n"
        f"nsubtype=0\n"
        f"output_path={output_path}\n"
    )
    paras_path.write_text(content)
    print(f"  Generated {paras_path}")
    return str(paras_path)


def build_matlab_command(species: str, paras_path: str) -> list:
    """Build MATLAB batch command for a species."""
    bundle_path = str(SRC_ROOT / species).replace("\\", "/")
    paras_path_m = paras_path.replace("\\", "/")

    matlab_cmd = (
        f"addpath('{bundle_path}'); "
        f"EpiProfile('{paras_path_m}')"
    )
    return [MATLAB_EXE, "-batch", matlab_cmd]


def run_species(species: str, paras_path: str) -> dict:
    """Run EpiProfile for one species. Returns result dict."""
    cmd = build_matlab_command(species, paras_path)

    output_dir = Path(paras_path).parent
    stdout_log = output_dir / f"matlab_stdout_{species}.log"

    print(f"\n  [{species}] Starting MATLAB...")
    print(f"  [{species}] Command: matlab -batch \"...EpiProfile('paras_{species}.txt')\"")

    t0 = time.time()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour max
        )
        elapsed = time.time() - t0

        # Save stdout/stderr
        with open(stdout_log, "w") as f:
            f.write(f"=== STDOUT ===\n{result.stdout}\n\n=== STDERR ===\n{result.stderr}\n")

        # Check for license error
        license_error = "license" in result.stderr.lower() if result.stderr else False

        print(f"  [{species}] Finished in {elapsed:.1f}s (exit code {result.returncode})")

        return {
            "species": species,
            "exit_code": result.returncode,
            "elapsed_seconds": elapsed,
            "success": result.returncode == 0,
            "license_error": license_error,
            "stdout_log": str(stdout_log),
        }
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        print(f"  [{species}] TIMEOUT after {elapsed:.1f}s")
        return {
            "species": species,
            "exit_code": -1,
            "elapsed_seconds": elapsed,
            "success": False,
            "license_error": False,
            "stdout_log": str(stdout_log),
        }


def verify_run(species: str, output_path: str, n_expected: int = 5) -> dict:
    """Verify a completed run produced expected output."""
    layout_dir = Path(output_path) / "histone_layouts"
    checks = {}

    # Check layout dir exists
    checks["layout_dir_exists"] = layout_dir.is_dir()

    # Count sample subdirs (NN_rawname pattern)
    sample_dirs = [d for d in layout_dir.iterdir() if d.is_dir() and d.name[:2].isdigit()] if layout_dir.is_dir() else []
    checks["n_samples"] = len(sample_dirs)
    checks["samples_match"] = len(sample_dirs) == n_expected

    # Check ratios file
    ratios_file = layout_dir / "histone_ratios_single_PTMs.xls"
    checks["ratios_exists"] = ratios_file.is_file()
    checks["ratios_size"] = ratios_file.stat().st_size if ratios_file.is_file() else 0

    # Check log
    log_file = Path(output_path) / "histone_logs.txt"
    checks["log_exists"] = log_file.is_file()

    # Check detail .mat files in first sample
    if sample_dirs:
        detail_dir = sample_dirs[0] / "detail"
        mat_files = list(detail_dir.glob("*.mat")) if detail_dir.is_dir() else []
        checks["n_mat_files"] = len(mat_files)
    else:
        checks["n_mat_files"] = 0

    checks["PASS"] = (
        checks["layout_dir_exists"] and
        checks["samples_match"] and
        checks["ratios_exists"] and
        checks["ratios_size"] > 0
    )

    return checks


def try_parallel(tasks: list) -> list:
    """Try running all species in parallel. Falls back to sequential on license error."""
    results = []

    with ProcessPoolExecutor(max_workers=len(tasks)) as executor:
        futures = {
            executor.submit(run_species, sp, pp): sp
            for sp, pp in tasks
        }

        for future in as_completed(futures):
            result = future.result()
            results.append(result)

            if result.get("license_error"):
                print(f"\n  LICENSE ERROR detected for {result['species']}!")
                print("  Falling back to sequential execution...")
                # Cancel remaining futures
                for f in futures:
                    f.cancel()
                # Run remaining species sequentially
                done_species = {r["species"] for r in results}
                remaining = [(sp, pp) for sp, pp in tasks if sp not in done_species]
                for sp, pp in remaining:
                    results.append(run_species(sp, pp))
                break

    return results


def run_sequential(tasks: list) -> list:
    """Run all species one at a time."""
    results = []
    for sp, pp in tasks:
        results.append(run_species(sp, pp))
    return results


def print_report(results: list, verifications: dict, total_elapsed: float):
    """Print summary report."""
    print("\n" + "=" * 70)
    print("  EpiProfile-PLANTS Multi-Species Run Report")
    print("=" * 70)
    print(f"\n{'Species':<8} {'Exit':<6} {'Time(min)':<10} {'Samples':<9} {'Ratios':<8} {'Status'}")
    print("-" * 60)

    for r in sorted(results, key=lambda x: x["species"]):
        sp = r["species"]
        v = verifications.get(sp, {})

        exit_str = str(r["exit_code"])
        time_str = f"{r['elapsed_seconds']/60:.1f}"
        samples_str = f"{v.get('n_samples', '?')}/5"
        ratios_str = "OK" if v.get("ratios_exists") else "MISSING"
        status = "PASS" if v.get("PASS") else "FAIL"

        print(f"{sp:<8} {exit_str:<6} {time_str:<10} {samples_str:<9} {ratios_str:<8} {status}")

    print(f"\nTotal elapsed: {total_elapsed/60:.1f} min")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description="EpiProfile-PLANTS multi-species launcher")
    parser.add_argument("--species", nargs="+", default=VALID_SPECIES,
                        help=f"Species to run (default: {VALID_SPECIES})")
    parser.add_argument("--raw-path", default=RAW_PATH,
                        help="Path to MS1/MS2 data")
    parser.add_argument("--runs-root", default=str(RUNS_ROOT),
                        help="Root directory for run outputs")
    parser.add_argument("--parallel", action="store_true", default=False,
                        help="Try parallel execution (default: sequential)")
    parser.add_argument("--n-samples", type=int, default=5,
                        help="Expected number of samples")
    args = parser.parse_args()

    # Validate
    for sp in args.species:
        if sp not in VALID_SPECIES:
            print(f"ERROR: Unknown species '{sp}'. Valid: {VALID_SPECIES}")
            sys.exit(1)

    if not Path(MATLAB_EXE).is_file():
        print(f"ERROR: MATLAB not found at {MATLAB_EXE}")
        sys.exit(1)

    print("=" * 70)
    print("  EpiProfile-PLANTS Multi-Species Launcher")
    print("=" * 70)
    print(f"  MATLAB:    {MATLAB_EXE}")
    print(f"  Raw data:  {args.raw_path}")
    print(f"  Output:    {args.runs_root}")
    print(f"  Species:   {', '.join(args.species)}")
    print(f"  Mode:      {'parallel' if args.parallel else 'sequential'}")
    print()

    # Generate configs
    print("Generating paras.txt files...")
    tasks = []
    for sp in args.species:
        out_dir = str(Path(args.runs_root) / sp)
        paras_path = generate_paras(sp, args.raw_path, out_dir)
        tasks.append((sp, paras_path))

    # Run
    print(f"\nLaunching {len(tasks)} species runs...")
    t0 = time.time()

    if args.parallel:
        results = try_parallel(tasks)
    else:
        results = run_sequential(tasks)

    total_elapsed = time.time() - t0

    # Verify
    print("\nVerifying outputs...")
    verifications = {}
    for r in results:
        sp = r["species"]
        out_dir = str(Path(args.runs_root) / sp)
        v = verify_run(sp, out_dir, args.n_samples)
        verifications[sp] = v
        status = "PASS" if v["PASS"] else "FAIL"
        print(f"  [{sp}] {status} — {v['n_samples']} samples, "
              f"{v['n_mat_files']} .mat files, "
              f"ratios={'OK' if v['ratios_exists'] else 'MISSING'}")

    # Report
    print_report(results, verifications, total_elapsed)


if __name__ == "__main__":
    main()
