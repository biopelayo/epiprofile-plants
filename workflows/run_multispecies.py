#!/usr/bin/env python3
"""
EpiProfile-PLANTS multi-species launcher.

Runs the AT, CR and MP bundles of this repository on one MS1/MS2 dataset and
keeps the outputs of each species apart.

How it works
------------
For every species the script:
  1. writes a `paras_<SP>.txt` (raw_path / norganism / nsource / nsubtype),
  2. launches `matlab -batch` with ONLY that species bundle on the path
     (`addpath(genpath('bundles/<SP>/src'))`, exactly as the README quickstart),
  3. moves `<raw_path>/histone_layouts/`, `<raw_path>/histone_ratios.xls` and
     `<raw_path>/histone_logs.txt` to `<runs_root>/<SP>/` once MATLAB exits.

Step 3 is needed because EpiProfile always writes next to the data
(`raw_path/histone_layouts`); it has no output-path parameter. That is also
why species run one after another: two bundles writing into the same
`histone_layouts/` at the same time would corrupt each other's outputs and
share the RT reference (`0_ref_info.mat`).

Requirements
------------
- MATLAB (R2019b+ for `-batch`). Location taken from, in order:
  `--matlab`, the `MATLAB_EXE` environment variable, `matlab` on PATH,
  `C:\\Program Files\\MATLAB\\<newest>\\bin\\matlab.exe`.
- A dataset folder laid out as EpiProfile expects:
      <raw_path>/<run>.raw   (may be 0-byte placeholders; only the names are read)
      <raw_path>/MS1/<run>.MS1
      <raw_path>/MS2/<run>.ms2

Usage
-----
    python workflows/run_multispecies.py --raw-path D:/data/MS1_MS2
    python workflows/run_multispecies.py --raw-path D:/data/MS1_MS2 --species AT CR
    python workflows/run_multispecies.py --raw-path D:/data/MS1_MS2 --runs-root D:/runs
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional, List

REPO_ROOT = Path(__file__).resolve().parents[1]
BUNDLES_ROOT = REPO_ROOT / "bundles"
VALID_SPECIES = ["AT", "CR", "MP"]


# ---------------------------------------------------------------------------
def find_matlab(explicit: Optional[str]) -> Optional[str]:
    """Locate matlab executable. Returns a path or None."""
    if explicit:
        return explicit if Path(explicit).is_file() else None
    env = os.environ.get("MATLAB_EXE")
    if env and Path(env).is_file():
        return env
    on_path = shutil.which("matlab")
    if on_path:
        return on_path
    if sys.platform.startswith("win"):
        cands = sorted(glob.glob(r"C:\Program Files\MATLAB\R*\bin\matlab.exe"))
        if cands:
            return cands[-1]
    return None


def check_dataset(raw_path: Path) -> List[str]:
    """Return the run basenames EpiProfile will see (from *.raw), with sanity checks."""
    raws = sorted(p.stem for p in raw_path.glob("*.raw"))
    if not raws:
        print(f"ERROR: no *.raw files in {raw_path}. EpiProfile derives the run list from "
              f"them (0-byte placeholders are enough when MS1/MS2 are already extracted).")
        sys.exit(1)
    missing = []
    for r in raws:
        ms1 = list((raw_path / "MS1").glob(f"{r}.[Mm][Ss]1"))
        ms2 = list((raw_path / "MS2").glob(f"{r}.[Mm][Ss]2"))
        if not ms1 or not ms2:
            missing.append(r)
    if missing:
        print("WARNING: these runs lack MS1/ or MS2/ files and will need RawToMS1.exe/xtract.exe "
              "next to the bundle to be converted: " + ", ".join(missing))
    return raws


def generate_paras(species: str, raw_path: Path, run_dir: Path) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    paras_path = run_dir / f"paras_{species}.txt"
    paras_path.write_text(
        "[EpiProfile]\n"
        "% the datapath of raw files (must contain <run>.raw, MS1/, MS2/)\n"
        f"raw_path={raw_path}\n\n"
        "% organism index of the bundle (all PLANTS bundles use 1)\n"
        "norganism=1\n\n"
        "% 1: histone_normal (the only mode ported to PLANTS)\n"
        "nsource=1\n\n"
        "% only used by nsource=4 (N15)\n"
        "nsubtype=0\n",
        encoding="ascii",
    )
    return paras_path


def build_matlab_command(matlab_exe: str, species: str, paras_path: Path) -> List[str]:
    bundle_src = (BUNDLES_ROOT / species / "src").as_posix()
    paras = paras_path.as_posix()
    matlab_cmd = (
        "restoredefaultpath; "
        f"addpath(genpath('{bundle_src}')); "
        "w = which('EpiProfile','-all'); "
        "assert(numel(w)==1, 'MATLAB path contaminated: %d EpiProfile.m on path', numel(w)); "
        f"EpiProfile('{paras}');"
    )
    return [matlab_exe, "-batch", matlab_cmd]


def run_species(matlab_exe: str, species: str, paras_path: Path, raw_path: Path,
                run_dir: Path, timeout_s: int) -> dict:
    cmd = build_matlab_command(matlab_exe, species, paras_path)
    stdout_log = run_dir / f"matlab_stdout_{species}.log"
    print(f"\n  [{species}] MATLAB -batch ... EpiProfile('{paras_path.name}')")
    t0 = time.time()
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
        rc, out, err = res.returncode, res.stdout, res.stderr
    except subprocess.TimeoutExpired as e:
        rc, out, err = -1, (e.stdout or ""), f"TIMEOUT after {timeout_s}s"
    elapsed = time.time() - t0
    stdout_log.write_text(f"=== STDOUT ===\n{out}\n\n=== STDERR ===\n{err}\n", encoding="utf-8",
                          errors="replace")

    # Move outputs out of raw_path so the next species starts clean.
    # (histone_ratios.xls is written next to the data, not inside histone_layouts/)
    moved = {}
    for name in ("histone_layouts", "histone_logs.txt", "histone_ratios.xls"):
        src = raw_path / name
        if src.exists():
            dst = run_dir / name
            if dst.exists():
                shutil.rmtree(dst) if dst.is_dir() else dst.unlink()
            shutil.move(str(src), str(dst))
            moved[name] = str(dst)

    print(f"  [{species}] exit {rc} in {elapsed/60:.1f} min; outputs -> {run_dir}")
    return {"species": species, "exit_code": rc, "elapsed_seconds": elapsed,
            "success": rc == 0, "stdout_log": str(stdout_log), "moved": moved,
            "license_error": "license" in (err or "").lower()}


def verify_run(run_dir: Path, n_expected: int) -> dict:
    layout_dir = run_dir / "histone_layouts"
    checks = {"layout_dir_exists": layout_dir.is_dir()}
    sample_dirs = [d for d in layout_dir.iterdir() if d.is_dir() and d.name[:2].isdigit()] \
        if layout_dir.is_dir() else []
    checks["n_samples"] = len(sample_dirs)
    checks["samples_match"] = len(sample_dirs) == n_expected
    ratios = run_dir / "histone_ratios.xls"
    checks["ratios_exists"] = ratios.is_file()
    checks["ratios_size"] = ratios.stat().st_size if ratios.is_file() else 0
    checks["single_ptms_exists"] = (layout_dir / "histone_ratios_single_PTMs.xls").is_file()
    checks["log_exists"] = (run_dir / "histone_logs.txt").is_file()
    if sample_dirs:
        detail = sample_dirs[0] / "detail"
        checks["n_mat_files"] = len(list(detail.glob("*.mat"))) if detail.is_dir() else 0
    else:
        checks["n_mat_files"] = 0
    checks["PASS"] = (checks["layout_dir_exists"] and checks["samples_match"]
                      and checks["ratios_exists"] and checks["ratios_size"] > 0)
    return checks


def print_report(results: list, verifications: dict, total_elapsed: float, n_expected: int):
    print("\n" + "=" * 70)
    print("  EpiProfile-PLANTS Multi-Species Run Report")
    print("=" * 70)
    print(f"\n{'Species':<8} {'Exit':<6} {'Time(min)':<10} {'Samples':<9} {'Ratios':<8} {'Status'}")
    print("-" * 60)
    for r in sorted(results, key=lambda x: x["species"]):
        v = verifications.get(r["species"], {})
        print(f"{r['species']:<8} {r['exit_code']:<6} {r['elapsed_seconds']/60:<10.1f} "
              f"{str(v.get('n_samples','?'))+'/'+str(n_expected):<9} "
              f"{'OK' if v.get('ratios_exists') else 'MISSING':<8} "
              f"{'PASS' if v.get('PASS') else 'FAIL'}")
    print(f"\nTotal elapsed: {total_elapsed/60:.1f} min")
    print("=" * 70)


def main():
    p = argparse.ArgumentParser(description="EpiProfile-PLANTS multi-species launcher",
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__)
    p.add_argument("--raw-path", required=True, type=Path,
                   help="dataset folder with <run>.raw placeholders, MS1/ and MS2/")
    p.add_argument("--species", nargs="+", default=VALID_SPECIES, choices=VALID_SPECIES,
                   help="species bundles to run (default: all)")
    p.add_argument("--runs-root", type=Path, default=None,
                   help="where per-species outputs go (default: <raw_path>/RUNS)")
    p.add_argument("--matlab", default=None, help="path to matlab executable")
    p.add_argument("--timeout", type=int, default=6 * 3600, help="seconds per species (default 6 h)")
    args = p.parse_args()

    raw_path = args.raw_path.resolve()
    if not raw_path.is_dir():
        print(f"ERROR: raw path not found: {raw_path}")
        sys.exit(1)
    runs_root = (args.runs_root or raw_path / "RUNS").resolve()

    matlab_exe = find_matlab(args.matlab)
    if not matlab_exe:
        print("ERROR: MATLAB not found. Use --matlab, set MATLAB_EXE, or add matlab to PATH.")
        sys.exit(1)

    for sp in args.species:
        if not (BUNDLES_ROOT / sp / "src" / "TIER2" / "EpiProfile.m").is_file():
            print(f"ERROR: bundle {sp} not found under {BUNDLES_ROOT}")
            sys.exit(1)

    runs = check_dataset(raw_path)
    if (raw_path / "histone_layouts").exists():
        print(f"WARNING: {raw_path/'histone_layouts'} already exists. It will be picked up (and "
              f"moved) by the first species run; delete it first for a clean run.")

    print("=" * 70)
    print("  EpiProfile-PLANTS Multi-Species Launcher")
    print("=" * 70)
    print(f"  MATLAB:    {matlab_exe}")
    print(f"  Raw data:  {raw_path}  ({len(runs)} runs)")
    print(f"  Output:    {runs_root}")
    print(f"  Species:   {', '.join(args.species)}  (sequential)")

    t0 = time.time()
    results, verifications = [], {}
    for sp in args.species:
        run_dir = runs_root / sp
        paras_path = generate_paras(sp, raw_path, run_dir)
        r = run_species(matlab_exe, sp, paras_path, raw_path, run_dir, args.timeout)
        results.append(r)
        verifications[sp] = verify_run(run_dir, len(runs))
        if r["license_error"]:
            print(f"  [{sp}] MATLAB reported a license problem, see {r['stdout_log']}")

    print_report(results, verifications, time.time() - t0, len(runs))
    sys.exit(0 if all(v["PASS"] for v in verifications.values()) else 1)


if __name__ == "__main__":
    main()
