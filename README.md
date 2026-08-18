# EpiProfile_PLANTS — Preprocessing, Data Model, Outputs, Provenance, RT References, and Troubleshooting

This document summarizes the recommended preprocessing pipeline (“gold-standard artifact generation”), deterministic MATLAB execution, the EpiProfile_PLANTS quantification data model (hDP/hPF/hPTM), expected outputs, provenance/audit conventions (T1–T4), RT reference hygiene, common failure modes, practical chemistry notes, and citation/license guidance.

---
<img width="1024" height="572" alt="image" src="https://github.com/user-attachments/assets/8b96045b-2299-4b97-a733-eba23f62410c" />

## ⚙️ Preprocessing pipeline (what the MATLAB core reads)

The MATLAB core does not read RAW, WIFF or mzML. It reads two text files per run,
placed in two subfolders of the data folder (`raw_path`):

```
<raw_path>/<run>.raw        run list only (0-byte placeholders are fine, or omit them
                            and the run list is taken from MS1/)
<raw_path>/MS1/<run>.MS1    centroided MS1 scans, pFind/RawToMS1 text layout
<raw_path>/MS2/<run>.ms2    MS2 scans, xtract text layout
```

`GetMS1ScanNo.m` / `GetMS2ScanNo.m` parse that layout literally
(`H\tDataType`, `S`, `I\tRetTime`, then `I\tIonInjectionTime` and `I\tInstrumentType`
lines, then `m/z intensity` pairs). Files written by other tools must reproduce it;
ProteoWizard's `--ms1/--ms2` writers do not (different `I` keys), so a converter is
needed for that route.

### Stage 1 — RAW → MS1/MS2 with the upstream extractors (validated route)

`Raw2MS.m` calls **`RawToMS1.exe`** (MS1) and **`xtract.exe`** (MS2) on Thermo `.RAW`
files. Both executables ship with the upstream EpiProfile 2.0 distribution and are
**not redistributed here** (their licence does not allow it): copy them into the folder
you run MATLAB from (`pwd`), and `EpiProfile` converts any run whose MS1/MS2 files are
missing. If the MS1/MS2 files already exist, the executables are not needed at all.

### Stage 2 — Other vendors (WIFF, .d) and mzML

Convert with ProteoWizard `msconvert` (vendor peak-picking ON, 64-bit, zlib) and then
write MS1/MS2 files in the layout above. There is no converter in this repository yet
(the earlier README pointed at `workflows/00_convert_*` and `workflows/01_extract_*`,
which were never added). Until one exists, use the extractor that produced your MS1/MS2
files and check the header keys against `GetMS1ScanNo.m`.

Common pitfall: mzML produced without vendor peak-picking yields profile spectra, and
`GetMS1ScanNo` stops with `MS1 is profile mode, convert to centroid mode first!`.

---

## ⚡ Quickstart (MATLAB, deterministic)

Requirements: MATLAB R2019b or newer (tested with R2023a on Windows 10, 2026-08-18),
Statistics and Machine Learning Toolbox (`boxplot`, `zscore`, `pca`), and Bioinformatics
Toolbox (`HeatMap`, `clustergram`) for the QC figures. Without the toolboxes set
`nfigure = 0` in `check_otherparas.m`; the ratio tables do not need them.

1. Lay out the data as shown above and write a `paras.txt` (template:
   [`paras.example.txt`](paras.example.txt) at the repository root):

   ```
   [EpiProfile]
   raw_path=C:\data\MS1_MS2
   norganism=1
   nsource=1
   nsubtype=0
   ```

2. In MATLAB, from the repository root, with **exactly one** species bundle on the path:

   ```matlab
   restoredefaultpath;                       % clean environment
   addpath(genpath("bundles/AT/src"));       % ONE bundle (AT, CR or MP)
   which EpiProfile -all                     % MUST return exactly one path
   EpiProfile("C:\data\MS1_MS2\paras.txt");  % or plain EpiProfile with paras.txt in pwd
   ```

3. Outputs land next to the data:
   - `<raw_path>/histone_ratios.xls`: cohort table (TSV inside), one block per hDP with
     Ratio / Area / RT(min) columns per run;
   - `<raw_path>/histone_layouts/histone_ratios_single_PTMs.xls`: site-level marginals;
   - `<raw_path>/histone_layouts/NN_<run>/detail/`: per-module `.mat` and XIC PDFs;
   - `<raw_path>/histone_layouts/0_ref_info.mat`: the RT reference (dataset-scoped);
   - `<raw_path>/histone_logs.txt`: the run diary.

Reference timing (laptop, R2023a, 2026-08-18): 5 *Arabidopsis* runs (SCIEX ZenoTOF 7600,
DDA, 100-130 MB MS1 and 100-180 MB MS2 text files) take about 20 min of MS1/MS2 parsing on
the first run (cached as `.mat` afterwards) plus 14.5 min of quantification for the full AT
panel (H3, H4, H2A, H2B, H1 modules), 306-line `histone_ratios.xls`.

To run the three bundles on the same dataset one after another, with outputs kept apart,
use `python workflows/run_multispecies.py --raw-path <raw_path>`.

# Operational rule
If `which EpiProfile -all` returns more than one entry, your MATLAB path is contaminated → results become non-deterministic. `nsource` other than 1 (SILAC, C13, N15, 13CD3) is refused: those upstream runners were not ported.

# 🧬 Data model (what you actually quantify)
EpiProfile_PLANTS uses three explicit layers:

1. **hDP (Histone-Derived Peptides)**  
   Backbone regions (e.g., H3_27_40, H4_4_17) and sequence/variant definitions.  
   Provides the “container” for families of modified peptideforms.

2. **hPF (Histone Peptideforms)**  
   Combinatorial PTM states on a given hDP.  
   Quantified by RT-aware extracted ion chromatograms (XIC/EIC) and AUC integration.

3. **hPTM (Site-level PTM summaries)**  
   Marginalization over hPFs to yield site-wise marks (e.g., H3K27me3, H4K16ac).  
   Computed via explicit include/exclude pattern logic (to prevent mis-assignment).

**Recommended:** Use validators to ensure sums do not accidentally include incompatible forms.

# 🧾 Outputs (what you should expect)
EpiProfile_PLANTS produces:

**A) Cohort-level matrices (analysis-ready)**  
- hDP matrix: peptide-level ratios/intensities.  
- hPF matrix: peptideform ratios/intensities.  
- hPTM matrix: site-level marginal ratios (derived from hPFs).

**B) QC/summary plots (cohort)**  
Typical set includes:  
- Identified peptide counts vs theoretical catalog size.  
- Intensity distributions (e.g., log10 peptide intensities per sample).  
- PCA of sample ratios.  
- Heatmaps of single-PTM ratios / z-scores.  
- Separate H3/H4 ratio heatmaps (when relevant).

**C) Per-sample “evidence” artifacts (audit trail)**  
- XIC plots per region/peptide.  
- PSM tables/lists and intermediate exports.  
- Identification lists and detail folders.

**Important file‑format note:** `histone_ratios.xls` is often a TSV-formatted text file (despite the extension). Downstream parsing should treat it as TSV and validate header/layout blocks.

# 🧾 Provenance & audit model (T1–T4)
All MATLAB functions are tracked under a four‑tier provenance model:

- **T1 (Reused):** Upstream functions used unchanged.  
- **T2 (Modified):** Upstream code adapted for plant workflows.  
- **T3 (New):** Newly implemented modules (e.g., species initializers such as `init_histone0_AT.m`).  
- **T4 (Not invoked):** Modules that ship in the bundle but no other file calls, kept for reference or future work.

Upstream functions that were never ported (SILAC, C13/N15 runners) are out of the bundle, not a tier.
T4 takes precedence over T1–T3, so the four tiers add up to the bundle. See `docs/tiers.md`.

**Source of truth:** `metadata/audit_master.tsv`

**Recommended practice:** Keep dataset-level manifests under `metadata/` so every figure/table can be mapped back to:  
- Bundle version  
- Layouts used  
- Conversion/extraction parameters  
- Code provenance tiering (T1–T4)

# 🕒 RT reference system (critical operational detail)
EpiProfile-style pipelines can use RT references to guide peak finding. This is powerful but easy to misuse.

**Key operational facts (high impact):**  
- RT references are effectively “global per folder” if you reuse the same `raw_path` across heterogeneous runs.  
- `check_ref` will typically look for: `raw_path/histone_layouts/0_ref_info.mat` and may apply it automatically depending on debug settings.  

**Failure mode: Cross-run RT contamination**  
If you mix instruments/columns/batches in the same folder, an RT reference learned from one run can be imposed on another → mis-quantification or missing peptides.

**Conservative operational recommendations (no code changes):**  
- One folder = one homogeneous group (same instrument/column/method family).  
- When reusing a folder for a new dataset, delete `0_ref_info.mat` (or regenerate cleanly).  
- If you want to run each RAW independently without applying a stored reference, use a “no-reference” mode (commonly `ndebug=2`, depending on your bundle conventions).

**Additional subtle contamination layer:**  
Some regions (notably `H3_27_40`-style layouts) may read RT anchors from region-specific files in the output path (e.g., `H3_04_27_40.xls` inside `His.outpath`). Reusing an outpath across runs can propagate RT assumptions.

**Bottom line:** Treat RT reference artifacts as **dataset-scoped**, not “project-scoped”.

# ⚠️ Known failure modes & troubleshooting

1. **MATLAB path contamination**  
   *Symptom:* `which EpiProfile -all` shows >1 hit.  
   *Fix:*  
   ```matlab
   restoredefaultpath;
   addpath(genpath("bundles/<ONE_BUNDLE>"));
   which EpiProfile -all
   ```

   ## Empty MS1/MS2

**Symptom:** near-zero identifications; missing chromatograms.  
**Usual cause:** mzML generated without vendor peak-picking.  
**Fix:**
- rerun `msconvert` with vendor peak-picking enabled
- validate mzML integrity/size and confirm non-empty scans

---

## RT window drift / missing peaks

**Symptom:** peptides not found despite MS evidence.  
**Fix:**
- widen curated RT windows (where appropriate)
- ensure RT references match the dataset
- avoid cross-run references (folder hygiene)

---

## Crashes due to empty RT windows

**Root cause class:** `rt_ref = 0` / out-of-range → empty scan windows → indexing failures.  
**Note:** In EpiProfile_PLANTS the intent is to harden these cases so pipelines degrade gracefully (zeros + logs), but RT reference hygiene remains first-line prevention.

---

## Path encoding / special characters

Avoid spaces/Unicode in dataset paths (MATLAB + legacy tooling are brittle).  
Prefer `C:\data\datasetA\` over `C:\Users\Mi Usuario\Datos\`.

---

## 🔬 Chemistry notes (practical)

EpiProfile_PLANTS is commonly used with derivatization-based histone workflows (often propionylation):
- Propionylation blocks free lysines, reducing charge-state variability and simplifying tryptic behavior.
- Trypsin digestion behaves effectively Arg-C–like (since Lys is blocked).
- Some protocols include a high-pH step to reverse over-derivatization on Ser/Thr/Tyr to avoid signal dispersion.
- Isobaric ambiguities exist and must be handled conservatively (validation rules, cutoffs, and careful PTM assignment logic).

---

## 📝 Citation & license

**License:**  
- GPL-3.0 (inherits from upstream EpiProfile 2.0 licensing constraints; check repository files for details)

**How to cite:**
- Cite EpiProfile_PLANTS via `CITATION.cff`
- Cite upstream: Zuo-Fei Yuan et al., EpiProfile 2.0

---

## 🔧 Optional repo helpers (templates)
<img width="1024" height="572" alt="image" src="https://github.com/user-attachments/assets/5536709d-12b9-4e60-86dd-a0ed2d843f2e" />



This repository can include (or you can generate) ready-to-use templates:
- `CITATION.cff` (GitHub citation UI compatible)
- `metadata/audit_master.tsv` scaffold (T1–T4 audit baseline)
- dataset manifest template (conversion parameters, extraction params, bundle hash/version, layouts used, RT reference mode/settings, output paths)

---

## Acknowledgements

EpiProfile_PLANTS builds upon the foundational work of **Dr. Zuo-Fei Yuan** and the **Garcia Lab** (University of Pennsylvania), who created [EpiProfile 2.0](https://github.com/zfyuan/EpiProfile2.0_Family) (Yuan et al., *J. Proteome Res.* 2018, 17, 2533--2541; [DOI: 10.1021/acs.jproteome.8b00133](https://doi.org/10.1021/acs.jproteome.8b00133)). We gratefully acknowledge Dr. Yuan's guidance during the development of the plant-specific extensions.
