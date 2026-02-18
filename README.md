# EpiProfile_PLANTS — Preprocessing, Data Model, Outputs, Provenance, RT References, and Troubleshooting

This document summarizes the recommended preprocessing pipeline (“gold-standard artifact generation”), deterministic MATLAB execution, the EpiProfile_PLANTS quantification data model (hDP/hPF/hPTM), expected outputs, provenance/audit conventions (T1–T4), RT reference hygiene, common failure modes, practical chemistry notes, and citation/license guidance.

---

## ⚙️ Preprocessing pipeline (gold-standard artifact generation)

### 🧱 Stage 1 — Conversion (RAW/WIFF → mzML)

Tool: **ProteoWizard `msconvert`** (recommended via **Docker**)

Required settings (typical, recommended):
- **Vendor peak-picking: ON**
- **64-bit**
- **zlib compression** (stable + portable)

Script:
- `workflows/00_convert_raw_or_wiff_to_mzml.*`

Common pitfall:
- Forgetting **vendor peak-picking** can produce mzML that later yields **empty / near-empty MS1/MS2**, leading to near-zero identifications downstream.

---

### 🧩 Stage 2 — Extraction (mzML → MS1/MS2)

Tool: **`xtract_xml.exe`** (external dependency)

Script:
- `workflows/01_extract_ms1_ms2_from_mzml.*`

Licensing note:
- `xtract_xml.exe` is treated as an **external dependency** and is **not redistributed** unless explicitly permitted.
- Treat it like a system tool:
  - available in `PATH`, **or**
  - referenced via configuration.

---

## ⚡ Quickstart (MATLAB, deterministic)

Deterministic execution requires the MATLAB path to include **exactly one** species bundle.

```matlab
% 1) Clean environment
restoredefaultpath;

% 2) Load a single bundle (example: Arabidopsis)
addpath(genpath("bundles/AT/src"));

% 3) Sanity check: MUST return EXACTLY one path
which EpiProfile -all

% 4) Launch
EpiProfile;
```

# Operational rule
If `which EpiProfile -all` returns more than one entry, your MATLAB path is contaminated → results become non-deterministic.

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
- **T4 (Excluded):** Intentionally removed / unused modules (e.g., SILAC/C13 modes, if not supported).

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

This repository can include (or you can generate) ready-to-use templates:
- `CITATION.cff` (GitHub citation UI compatible)
- `metadata/audit_master.tsv` scaffold (T1–T4 audit baseline)
- dataset manifest template (conversion parameters, extraction params, bundle hash/version, layouts used, RT reference mode/settings, output paths)
