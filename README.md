# EpiProfile-PLANTS `v1.1-tesis`: the AT bundle as run for the thesis

This branch holds the MATLAB code that quantified the histone peptidoforms of the doctoral thesis
behind EpiProfile-PLANTS, at the version the thesis cites, `v1.1-tesis`. It is an orphan branch
with no shared history with `main`, and it is frozen. The maintained bundles, their documentation
and later fixes live on `main`.

- `v1.0-tesis-asrun`: the 119 `.m` files exactly as they ran on 2026-04-16 (not published).
- `v1.1-tesis` (this branch): `v1.0` plus a guard in `check_ref.m` against an `rt_ref` of 0 returned
  by the calibration. When the calibration finds every reference peptide, the output is the same as
  with `v1.0`.

The code derives from EpiProfile 2.0 (Yuan et al., *J. Proteome Res.* 2018, 17, 2533-2541) and is
distributed under GPL-3.0 (`LICENSE`).

## Contents

| Path | What |
|---|---|
| `src/*.m` | the 119 MATLAB files of `v1.1-tesis`, byte for byte |
| `src/paras.txt` | the `paras.txt` of the thesis runs with a generic `raw_path` (only that value changed) |
| `MANIFEST.sha256` | SHA-256 of every file under `src/` |
| `.gitattributes` | `* -text`, so Git never rewrites line endings and the manifest keeps matching |

Check a checkout with `sha256sum -c MANIFEST.sha256` (Git Bash, Linux, macOS) or compare
`Get-FileHash -Algorithm SHA256` in PowerShell.

No data, logs, results or executables are included. The data of the thesis will be published with
its deposit.

## How to run it

Requirements, installation, input layout, `paras.txt`, outputs and known pitfalls are documented in
`docs/INSTALL.md` and `docs/RUNNING.md` of the main line of the repository
([main](https://github.com/biopelayo/epiprofile-plants/tree/main/docs); until pull request #1 is
merged,
[update-operability-2026-08](https://github.com/biopelayo/epiprofile-plants/tree/update-operability-2026-08/docs)).
What is specific to this version:

- Tested environment: Windows 10, MATLAB R2023a. This version writes the QC figures by default
  (`nfigure = 1` in `src/check_otherparas.m`), so it needs the Curve Fitting Toolbox, the
  Statistics and Machine Learning Toolbox and the Bioinformatics Toolbox.
- Every run needs a `<raw_path>/<run>.raw` file (a 0-byte placeholder is enough when
  `MS1/<run>.MS1` and `MS2/<run>.ms2` already exist). The run list comes only from those files.
- Start from a `raw_path` without a `histone_layouts/` folder, and run from `src/`, where the
  extractors must sit if `.raw` files still need converting:

  ```
  cd <checkout>\src
  "C:\Program Files\MATLAB\R2023a\bin\matlab.exe" -batch "restoredefaultpath; addpath(pwd); EpiProfile('paras.txt')"
  ```

  after setting `raw_path` in `src/paras.txt`.
- With 100 runs or more in one folder, or with `ndebug` other than 0, `DrawISOProfile0.m` builds no
  RT reference and the run continues without it.

## Extractors (not redistributed)

`RawToMS1.exe` and `xtract.exe` convert Thermo `.raw` files into the MS1/MS2 text files that the
code reads. Their licence does not allow redistribution here. Both ship inside the upstream
packages at <https://github.com/zfyuan/EpiProfile2.0_Family>. Unpacked copies of
`EpiProfile2.1_1Basic.zip` and `EpiProfile2.2.zip` carry the exact files used for the thesis:

| File | Version | SHA-256 | MD5 |
|---|---|---|---|
| `RawToMS1.exe` | 1.0.0.1 | `27b162ddc6a8fcdfe461efad5e6112b4ad1d1bbe34002f0f0c63adb016474be8` | `a0bd03d9302643f3437463be5fbcfffe` |
| `xtract.exe` | pXtract 2.1.0 | `5f9d157a9a4de3f531ef367e8f91c174a72e639d0b608d14211240424b75188b` | `d229889469c16db4bad48343a0d88baf` |

`EpiProfile2.0_1Basic.zip` carries an older `xtract.exe` with another hash. Given a single `.raw`
file, this `xtract.exe` writes no MS2 and still returns 0. Extract MS2 by passing the folder
instead (see `docs/RUNNING.md`, section 3). Data from other vendors (SCIEX `.wiff`) do not use these
executables; the route is described in the same section.
