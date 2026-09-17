# Running EpiProfile-PLANTS

This page covers one run from spectra to ratio tables with a species bundle of this repository.
Installation and requirements: [`INSTALL.md`](INSTALL.md). The run loop, `paras.txt` and the output
names are those of upstream EpiProfile 2.0; the differences are listed where they matter.

Statements marked *checked in MATLAB* were executed with MATLAB R2023a on 2026-09-17; the rest
come from reading the code of `bundles/AT/src` (the CR and MP bundles share the same entry point and
parsers).

## 1. Data layout

Everything for one dataset lives under one folder, `raw_path`, and all outputs are written there:

```
<raw_path>/
    <run>.raw          optional. Only the NAME is read when MS1/MS2 exist (0-byte files are fine).
                       Without any *.raw, the run list is taken from MS1/*.MS1 (checked in MATLAB).
    MS1/<run>.MS1      centroided MS1 scans (.ms1 also accepted)
    MS2/<run>.ms2      MS2 scans (.MS2 also accepted)
    paras.txt          may live here or anywhere else
```

- Run names are sorted, and each run gets a two-digit prefix in that order (`01_<run>`, `02_<run>`).
- Keep one acquisition regime per `raw_path` (same instrument, column and gradient). The retention
  time reference built from the first runs is applied to every run in the folder.
- Start every run from a folder **without** `histone_layouts/`. `0_ref_info.mat` is reused if it
  exists, and `OutputTogether.m` builds the table from the module files it finds under
  `histone_layouts/01_<run>/detail/`, so leftovers from an earlier run end up in the new table.
- The first run writes `<run>_MS1scans.mat`, `<run>_MS1peaks.mat` (in `MS1/`) and
  `<run>_MS2scans.mat`, `<run>_MS2peaks.mat` (in `MS2/`) and reuses them afterwards without
  looking at the text files again. Delete them whenever the MS1/MS2 text files change.

## 2. MS1/MS2 text format

`GetMS1ScanNo.m` and `GetMS2ScanNo.m` read the pFind text layout written by `RawToMS1.exe` and
`xtract.exe`:

One MS2 scan (fields separated by tabs; the peak lines by a space):

```
H	DataType	Centroided
S	000361	000361	311.8268403
I	NumberOfPeaks	804
I	RetTime	133.683000
I	IonInjectionTime	0.000000
I	ActivationType	HCD
I	InstrumentType	FTMS
I	PrecursorScan	360
I	ActivationCenter	311.83
I	MonoiosotopicMz	311.8268403
Z	2	622.646404
100.05055 8.0
100.06298 3.7
```

- `H` lines form the file header. `H DataType Profile` stops the run.
- `S` starts a scan (MS2: scan number twice, then the precursor m/z). `Z` gives charge and MH+.
- `I RetTime` is read in seconds and converted to minutes when the last MS1 scan is above 1000.
  `I InstrumentType ITMS` in the MS1 file switches the tolerance to 1000 ppm.
- MS2 is parsed **by position**: after each `S` line the parser expects those eight `I` lines in
  that order, then the `Z` line(s), then the peaks. MS2 files need the `H DataType` line.
- MS1 is more tolerant: it accepts `I RTime` (ProteoWizard) instead of `I RetTime`, skips unknown
  `I` lines and assumes centroid data when `H DataType` is missing.
- ProteoWizard's own `--ms1`/`--ms2` output does not follow the MS2 layout above.

## 3. Producing MS1/MS2

### Thermo `.raw`

Put the `.raw` files in `raw_path` and `RawToMS1.exe` plus `xtract.exe` in the folder MATLAB runs
from. When a run lacks its MS1 or MS2 file, `EpiProfile` calls `Raw2MS.m`, which runs:

```
RawToMS1.exe "<raw_path>\<run>.RAW"                      -> MS1\<run>.MS1
xtract.exe -a -i 1 -m 5 -ms2 -o "<raw_path>\MS2" "<raw_path>\<run>.RAW"
```

Known problem (pXtract 2.1.0, seen on 2026-08-18): given a single `.raw` file, `xtract.exe` prints
`Find no raws in <file>`, returns 0 and writes nothing, so the run stops later in `GetMS2ScanNo`.
Given the **folder**, it extracts every file. Extract MS2 beforehand from the folder that holds
`xtract.exe`, then start `EpiProfile`, which skips the conversion because the `.ms2` files exist:

```
xtract.exe -a -i 1 -m 5 -ms2 -o "<raw_path>\MS2" "<raw_path>"
```

`RawToMS1.exe` works file by file.

### SCIEX `.wiff` (the route of the thesis cohort)

1. Convert to mzML with vendor centroiding, for example
   `msconvert <run>.wiff --mzML --64 --zlib --filter "peakPicking vendor"`. Without vendor peak
   picking the MS1 are profile spectra and the run stops with
   `MS1 is profile mode, convert to centroid mode first!`.
2. Write MS1/MS2 with pFind's `xtract_xml.exe` (version 3.0.0 in the thesis):
   `xtract_xml.exe -ms -o <out_folder> <run>.mzML`. It writes `<run>.ms1` and
   `<run>.HCD.FTMS.ms2`.
3. Rename the MS2 file to `<run>.ms2`, then move the files into `MS1/` and `MS2/`.

Two cautions. The files written this way declare `Vendor Thermo` in their header whatever the
instrument; the parser ignores that line. The same extractor applied to Thermo Orbitrap Elite mzML
produced MS1 intensities that overflow (see [`../workflows/PXD014739/README`](../workflows/PXD014739/README)),
so use the `.raw` route above for Thermo data and look at the intensity range of the first MS1 file
before a long run.

## 4. `paras.txt`

Template: [`../paras.example.txt`](../paras.example.txt).

```
[EpiProfile]
% folder with MS1/, MS2/ (and optionally <run>.raw); every output is written here
raw_path=C:\data\MS1_MS2

% organism index inside the bundle; every PLANTS bundle uses 1
norganism=1

% 1 = histone_normal, the only mode ported to PLANTS (2 SILAC, 3 C13, 4 N15, 5 13CD3 are refused)
nsource=1

% sub-mode, only read by nsource=4; keep 0
nsubtype=0
```

How `ReadInput.m` reads it (the first three points *checked in MATLAB*):

- `raw_path=` must start in column 1. An indented line leaves `raw_path` empty, and the run stops
  with `raw_path does not exist`.
- The other keys are searched in the order `norganism`, `nsource`, `nsubtype`, anywhere in a line.
  A comment line **below** `raw_path` that contains one of those words is taken as the key line:
  the value becomes `NaN` and no error is raised. Keep key names out of comments.
- The value is everything after `=`, spaces included. Do not leave trailing spaces after the path.
- The data path may contain spaces; the folder MATLAB runs from should not when `.raw` conversion
  is needed (the `.exe` path is not quoted).

## 5. Run-time switches (`TIER2/check_otherparas.m`)

These are set in code, not in `paras.txt`. Bundle defaults *checked in MATLAB*.

| Switch | Bundle default | Values |
|---|---|---|
| `def_ptol` | `10` | MS1 tolerance in ppm (set to 1000 automatically for ITMS data) |
| `soutput` | `'31'` | first digit: `1` H3/H4 basic, `2` basic + H3S10ph, `3` all H3/H4 modules; second digit: `1` adds H1/H2A/H2B, `0` leaves them out |
| `nfigure` | `0` | `1` writes the QC figures (needs the Statistics and Bioinformatics toolboxes; only with 3 or more runs) |
| `ndebug` | `0` | `0` normal run with an RT reference; `1` manual XIC correction; `2` no RT reference |

Limitation of this bundle: `EpiProfile.m` stops with
`ERROR: calibration failed (0_ref_info.mat was not created)` whenever `DrawISOProfile0.m` does not
write the reference, and `DrawISOProfile0.m` (upstream code) skips it when `ndebug` is not 0 or the
folder holds 100 runs or more. Until that is changed, keep `ndebug = 0` and split large datasets into
folders of fewer than 100 runs of the same regime. (Found by reading the code, not by a failing run.)

## 6. Minimal run

Interactive, from the repository root:

```matlab
restoredefaultpath;
addpath(genpath('bundles/AT/src'));
which EpiProfile -all                      % exactly one hit
EpiProfile('C:\data\MS1_MS2\paras.txt')    % or EpiProfile with paras.txt in the current folder
```

Unattended (Windows command prompt). Use an absolute bundle path, because the current folder must
be the one with the `.exe` files when conversion is needed:

```
cd C:\tools\extractors
"C:\Program Files\MATLAB\R2023a\bin\matlab.exe" -batch "restoredefaultpath; addpath(genpath('C:\tools\epiprofile-plants\bundles\AT\src')); EpiProfile('C:\data\MS1_MS2\paras.txt')"
```

Three bundles on the same data, one after another, each with its own output folder:
`python workflows/run_multispecies.py --raw-path C:\data\MS1_MS2`.

What happens, in order: `Raw2MS` (only if an MS1/MS2 file is missing) → `GetMS1ScanNo` and
`GetMS2ScanNo` (text to `.mat` caches) → `DrawISOProfile0` (RT reference, opens the diary) →
`DrawISOProfile1` (every module for every run; DDA or DIA detected per run) → `OutputTogether` →
`OutputSinglePTMs` → `OutputFigures` (if `nfigure = 1`).

Reference timing (laptop, R2023a, 2026-08-18): 5 *Arabidopsis* runs from a SCIEX ZenoTOF 7600 (DDA,
100-180 MB per MS1 or MS2 file) took about 20 min to parse the first time and 14.5 min to quantify
with the full AT panel.

## 7. Outputs

All under `raw_path` (names checked against a complete run):

| Path | Content |
|---|---|
| `histone_ratios.xls` | The cohort table. Tab-separated text with CRLF line ends despite the extension. Row 1: `i,<run>` headers, repeated for three blocks. Row 2: `Peptide`, then `Ratio`, `Area` and `RT(min)` blocks separated by an empty column. Then, per module, one line `<sequence>(<region>)` followed by one line per peptidoform, `<region> <modification>`, with its ratio, area and RT in every run |
| `histone_logs.txt` | MATLAB diary from the RT reference onwards; one block per run and the elapsed time at the end |
| `histone_layouts/0_ref_info.mat` | RT reference of this dataset |
| `histone_layouts/NN_<run>/<module>.pdf` | extracted-ion chromatogram layout per module |
| `histone_layouts/NN_<run>/detail/<module>.xls`, `.mat` | per-module quantification (the `.mat` files feed `histone_ratios.xls`) |
| `histone_layouts/NN_<run>/detail/psm/` | `<module>.mat`, `<module>.plabel` and `identification_list.xls` (DDA runs) |
| `histone_layouts/histone_ratios_single_PTMs.xls`, `.mat` | one row per single modification, one ratio column per run |
| `histone_layouts/heatmap_clustering/` | with `nfigure = 1`: `01_bar_peptide_number.pdf`, `02_boxplot_peptide_intensity.pdf`, `03_pca.pdf`, `10_heatmap_ratio_single_PTMs.pdf`, `10_heatmap_zscore_single_PTMs.pdf`, `11_heatmap_ratio_<H1,H2A,H2B,H3,H4>.pdf`, `Figure Legends.txt` |
| `MS1/*_MS1*.mat`, `MS2/*_MS2*.mat` | parsing caches (section 1) |

A failed run usually shows in `histone_logs.txt` (a run block without its module list, or no
`elapsed time` line) before it shows in the table.

## 8. Which code reproduces the thesis

The ratio matrices of the doctoral thesis were **not** produced with `bundles/AT/src` of this
repository. They come from the as-run AT bundle (119 `.m` files, `RawToMS1.exe`, `xtract.exe` and
the `paras.txt` of the run), which the author keeps with a SHA-256 manifest and three tags:

| Tag | Change | Effect |
|---|---|---|
| `v1.0-tesis-asrun` | code as it ran on 2026-04-16 | reference matrix |
| `v1.1-tesis` | guard in `check_ref.m` against `rt_ref = 0` from the calibration | same matrix for the reference cohort; recovers `H4_20_23` in datasets where the calibration returned 0 |
| `v1.2-tesis` | `DrawISOProfile0.m` threshold raised from 100 to 10000 runs | identical to `v1.1-tesis` below 100 runs; calibrates datasets with 100 runs or more |

Those tags live in a separate archival repository that is not published yet. Until it is, results
from this repository should not be expected to match the thesis tables.

How this tree differs from `v1.2-tesis` (file contents compared after normalising line endings,
2026-09-17):

- 117 files have the same name in both; 102 are identical and 15 differ: `DrawISOProfile0`,
  `DrawISOProfile1`, `EpiProfile`, `GetBenchmark`, `H3_07_73_83`, `HH2A_AT_Snapshot`,
  `HH2B_AT_Snapshot`, `OutputFigures`, `OutputSinglePTMs`, `OutputTogether`, `check_otherparas`,
  `draw_layout`, `init_histone0`, `output_histone`, `output_histone2`.
- `v1.2-tesis` also runs `HH2A_AT09_can_36_42` and `HH2B_AT08_shared_80_86`, which this tree leaves
  out (retired on 2026-05-24).
- This tree adds 113 files: 112 in `TIER4/` (the 108 generated H1/H2A/H2B catalogue modules and
  snapshots, plus `H3_04v3a_27_40`, `H3_05b_41_49`, `H4_03_24_35` and `HH2B_01u_104_145`, none of
  them called) and `TIER3/run_calibrate.m`.
- Defaults: the as-run bundle writes the QC figures (`nfigure = 1`); this tree does not.
