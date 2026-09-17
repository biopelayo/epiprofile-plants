# Installing EpiProfile-PLANTS

EpiProfile-PLANTS is a set of MATLAB bundles (`bundles/AT`, `bundles/CR`, `bundles/MP`) derived
from EpiProfile 2.0 basic (Yuan et al., *J. Proteome Res.* 2018). There is nothing to compile: a
run needs MATLAB with one toolbox, one bundle on the MATLAB path, and the spectra already written
as MS1/MS2 text files (or Thermo `.raw` files plus two upstream extractors).

What to do next: [`RUNNING.md`](RUNNING.md) (data layout, `paras.txt`, a minimal run, outputs and
which code reproduces the thesis).

## How the requirements below were checked

- **Toolboxes:** `matlab.codetools.requiredFilesAndProducts` (MATLAB R2023a Update 2,
  2026-09-17) over every `.m` file of `bundles/AT/src`, `bundles/CR/src`, `bundles/MP/src`, the
  thesis as-run bundle and upstream EpiProfile 2.0 basic. All five trees give the same answer.
- **Parameters and file layout:** read from the code (`EpiProfile.m`, `ReadInput.m`,
  `check_otherparas.m`, `Raw2MS.m`, `GetMS1ScanNo.m`, `GetMS2ScanNo.m`, `Output*.m`) and checked
  in MATLAB where noted in `RUNNING.md`.
- **End-to-end runs:** Windows 10, MATLAB R2023a, *Arabidopsis* SCIEX ZenoTOF 7600 DDA data
  (5 runs from this branch on 2026-08-18; the thesis runs with the as-run bundle in 2026).

## Requirements

| Component | Version | Needed for | Status |
|---|---|---|---|
| MATLAB | R2019b or later | everything | Tested only with **R2023a** (9.14, Update 2). R2019b is the declared floor, not a tested one (see note 1) |
| **Curve Fitting Toolbox** | matching MATLAB | **every run** | Required. `smooth` is called by `get_area.m` (reached from every quantification module), `get_histone2/4/6.m` and `get_rts2/22.m` |
| Statistics and Machine Learning Toolbox | matching MATLAB | QC figures only (`nfigure = 1`) | `OutputFigures.m`: `zscore`, `boxplot`, `pca`, `kmeans` |
| Bioinformatics Toolbox | matching MATLAB | QC figures only (`nfigure = 1`) | `OutputFigures.m`: `HeatMap`, `clustergram` |
| Operating system | Windows 10 | RAW conversion; all tested runs | The MATLAB code makes no OS-specific call except `Raw2MS.m`, which runs Windows `.exe` files. Linux and macOS with pre-extracted MS1/MS2 are untested |
| `RawToMS1.exe` and `xtract.exe` | from the EpiProfile 2.0 distribution | converting Thermo `.raw` only | Not redistributed here (see note 2) |
| ProteoWizard `msconvert` | any recent release | non-Thermo vendors (SCIEX `.wiff`) | Produces mzML; a second extractor writes MS1/MS2 (see `RUNNING.md`, section 3) |
| Python 3 | tested with 3.13 | optional helpers: `tools/*.py`, `workflows/run_multispecies.py` | Not needed by the MATLAB pipeline |

Disk: plan for roughly 200 MB per run for the MS1/MS2 text files and the `.mat` caches that the
first run writes next to them (SCIEX ZenoTOF 7600 DDA), plus a few MB per run under
`histone_layouts/`.

**Note 1, MATLAB version.** A static scan of the TIER1-TIER3 code of the AT bundle found no function
introduced after R2016b. The command lines in these documents use double-quoted strings (R2017a or
later) and `matlab -batch` (R2019a or later). Nothing older than R2023a has been run.

**Note 2, extractors.** Both executables come with upstream EpiProfile 2.0. The thesis used
`RawToMS1.exe` 1.0.0.1 (MD5 `a0bd03d9302643f3437463be5fbcfffe`) and `xtract.exe` pXtract 2.1.0
(MD5 `d229889469c16db4bad48343a0d88baf`). Different files with the same name and version string
circulate, so compare the MD5, not the name. The upstream manual adds two warnings that were not
re-tested here: Thermo MSFileReader may be needed for `.raw` access, and `xtract.exe` stops working
at the end of each year until it is replaced by a newer copy.

## Install

1. Clone the repository into a path without spaces:

   ```
   git clone https://github.com/biopelayo/epiprofile-plants.git C:\tools\epiprofile-plants
   ```

2. Install MATLAB with the Curve Fitting Toolbox (add the Statistics and Machine Learning Toolbox
   and the Bioinformatics Toolbox if you want the QC figures).

3. Only for Thermo `.raw` input: copy `RawToMS1.exe` and `xtract.exe` into the folder MATLAB will
   run from (`pwd`). `Raw2MS.m` looks for them there and builds the command without quotes, so that
   folder must not contain spaces.

4. Check the installation in MATLAB, from the repository root:

   ```matlab
   restoredefaultpath;
   addpath(genpath('bundles/AT/src'));          % ONE bundle: AT, CR or MP
   which EpiProfile -all                        % must print exactly one path
   which smooth                                 % must point into toolbox\curvefit
   license('test','Curve_Fitting_Toolbox')      % 1
   license('test','Statistics_Toolbox')         % 1 if you want figures
   license('test','Bioinformatics_Toolbox')     % 1 if you want figures
   ```

   If `which EpiProfile -all` prints more than one path, another EpiProfile copy is on the path and
   the run will pick modules from both. Start again from `restoredefaultpath`.

Never add two bundles (or a bundle and an upstream EpiProfile folder) to the path at the same time:
module names are shared, and MATLAB silently uses the first one it finds.
