# Changelog

All notable changes to EpiProfile_PLANTS. Dates are ISO (YYYY-MM-DD).

## [Unreleased] — branch `update-operability-2026-08`

Operability audit of the public repository (2026-08-18, MATLAB R2023a, Windows 10). A clean
clone passed the static checks but crashed on real *Arabidopsis* data inside the first H3
module (`get_histone12.m` line 12, empty RT window), because the published core lacked the
guards and fixes that the thesis bundle carries. See `_qc_local/` in the working copy for the
full report (not versioned).

### Fixed
- `EpiProfile.m` (all bundles): refuses `nsource` other than 1 with a message instead of
  failing with "Undefined function DrawISOProfile3/4/5" after 20 min of parsing; calls
  `Raw2MS` only when some MS1/MS2 file is missing; accepts `.MS1/.ms1` and `.ms2/.MS2`;
  clear message when `paras.txt` or `raw_path` is missing. Carries the 2026-05-24 fixes:
  no `clc` (hung `-batch`), startup timestamp, abort when `0_ref_info.mat` was not written.
- `check_otherparas.m` (all bundles): run list falls back to `MS1/*.MS1` when there is no
  `*.raw`; sorted, de-duplicated run list (deterministic `01_`, `02_` numbering);
  `nfigure = 0` by default (figures need two toolboxes and a display).
- CR and MP bundles: the generic core files gained the empty-window / `rt_ref<=0` guards
  (`get_histone1/2/10/12/13/22`, `get_rts0/2/22`, `check_ref`), stable `dir()` ordering
  (`GetBenchmark`, `OutputTogether`, `OutputFigures`), NaN below the 1e3 intensity threshold
  (`output_histone`, `output_histone2`), tolerant MS1/MS2 parsers (`GetMS1ScanNo` reads
  ProteoWizard headers, `GetMS2ScanNo` skips malformed peak lines) and the `H3_Snapshot`
  `out_filename` fix (`H3_06_54_63` → `H3_06_53_63`). Those 18 files moved from `TIER1/` to
  `TIER2/` because they now differ from upstream. Byte-identical duplicates
  `TIER3/H3_06_53_63.m` and `TIER3/H3_Snapshot.m` removed (they produced two hits in `which`).
- `bundles/AT/scripts/extract_empirical_rt.py` did not compile (`global` after use).
- `bundles/AT/scripts/*.py` and `workflows/run_multispecies.py` used absolute paths of the
  author's workstation (`D:/Antigravity/...`, `E:/EpiProfile_Proyecto/...`, `SRC_update/`);
  they are now repository-relative. `run_multispecies.py` also ran the three species on the
  same `histone_layouts/`; it now runs them one after another and moves each species' outputs
  to `<runs_root>/<SP>/`.
- `.gitignore`: `*.xls` no longer hides the reference tables under `workflows/`.
- README requirements (2026-09-17): the Curve Fitting Toolbox is required for every run
  (`smooth`); the Statistics and Machine Learning and Bioinformatics toolboxes only for the QC
  figures. The earlier text said no toolbox was needed without figures.
- Workstation paths removed from `tier_audit_<XX>.tsv` (`tools/tier_audit.py` now writes
  `upstream_file` relative to `--upstream`), `plant_module_manifest.tsv`,
  `build_AT_nucleosome_master.py`, MANUAL §55.2 and `workflows/PXD014739/README`.

### Added
- `docs/INSTALL.md` and `docs/RUNNING.md` (2026-09-17): requirements checked with
  `matlab.codetools.requiredFilesAndProducts`, MS1/MS2 format as parsed, Thermo and SCIEX
  conversion routes, commented `paras.txt`, run-time switches, outputs, and which code
  reproduces the thesis matrices.
- `paras.example.txt` at the root; `paras.txt` was undocumented.
- `tools/tier_audit.py`: provenance (T1/T2/T3 vs upstream) and invocation (T4) audit of a
  bundle, implementing `docs/tiers.md`.
- README: what the MATLAB core actually reads (`MS1/`, `MS2/`, pFind text layout), where the
  external extractors come from, required toolboxes, tested quickstart with `paras.txt`,
  measured timings.

### Changed
- MANUAL §10.3 smoke test (called `DrawISOProfile1` with `special = 0`, which cannot work) and
  §13 naming conventions (`<run>_MS1.txt` in one folder is not what the code reads) rewritten
  to match the code.
- `metadata/species.tsv`: CR and MP are assembled bundles (not "planned"), not yet validated
  on data.
- Consolidated the 2026-07/08 sequence-reference work (real UniProt sequences in
  `assets/sequences/`, H2A/H2B/H1 catalogues, `tiers.md` T4 definition, PXD014739 provenance).

### Known gaps (not fixed here)
- `EpiProfile.m` aborts when `0_ref_info.mat` is missing, but `DrawISOProfile0.m` does not write
  it with `ndebug` other than 0 or with 100 runs or more in one folder, so those cases stop after
  parsing (`docs/RUNNING.md`, section 5).
- No converter from mzML to the MS1/MS2 text layout ships in the repository yet; the
  README no longer points at the non-existent `workflows/00_convert_*` / `01_extract_*`.
- `metadata/audit_master.tsv` predates the AT bundle expansion; regenerate it with
  `tools/tier_audit.py` once the AT core decision below is taken.

## [0.1.0] — 2026-02-20 … 2026-08-03
Initial public bundles for AT, CR and MP; H2A/H2B/H1 catalogue modules for AT (2026-03-27);
auxiliary functions catalogue and documentation updates.
