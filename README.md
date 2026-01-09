# 🌿 EpiProfile_PLANTS

**EpiProfile_PLANTS** is a plant-oriented, auditable, and reproducible extension of **EpiProfile 2.0** (Zuo-Fei Yuan et al., GPL-3.0) designed for precise and transparent quantification of histone modifications in botanical models.

It supports the quantification of:

- 🧬 **hDP** — *Histone-derived peptides* (peptide backbones)
- 🧩 **hPF** — *Histone peptideforms* (PTM combinatorial states)
- 🧷 **hPTM** — *Site-level histone PTMs* (marginal summaries derived from hPFs)

---

## 🚀 Overview

### ✅ Gold-standard input (the “computational truth”)

Unlike legacy workflows, EpiProfile_PLANTS defines the **true computational input** using three mandatory artifacts:

- 📄 `*.mzML` — high-fidelity conversion output
- 📄 `*.MS1` — extracted ion peaks
- 📄 `*.MS2` — fragment evidence

Optional legacy compatibility:
- 🧯 `*.raw` — optional dummy placeholder for upstream/legacy expectations

> **Key principle:** quantification is driven by **mzML + MS1 + MS2**.  
> `.raw` is **not required** for quantification.

---

## 🔁 Canonical workflow

```text
RAW/WIFF  ->  (Docker + msconvert)  ->  mzML  ->  (xtract_xml.exe)  ->  MS1/MS2  ->  EpiProfile_PLANTS
```

### Mermaid diagram (optional; works in GitHub)

```mermaid
graph LR
    A[RAW/WIFF] --> B(Docker + msconvert)
    B --> C[mzML]
    C --> D(xtract_xml.exe)
    D --> E[MS1/MS2]
    E --> F[EpiProfile_PLANTS]
```

---

## ✨ Key differentiators

| Feature | EpiProfile 2.0 (Original) | EpiProfile_PLANTS |
|---|---|---|
| Target species | General / human-centric | Plant-specific (AT, MP, CR, PP) |
| Input handling | Implicit / variable | **Gold-standard** (`mzML + MS1/MS2`) |
| Runtime model | Shared MATLAB path | **Isolated species bundles** |
| Auditability | Manual / ad hoc | **T1–T4 provenance model** |
| Chemistry focus | Multiple chemistries | **Anhydride derivatization (BottomMap-like)** |

---

## ⚙️ Preprocessing pipeline (gold-standard input generation)

Reproducibility depends on generating the gold-standard artifacts correctly. Follow these two stages:

### 🧱 Stage 1 — Conversion (RAW/WIFF → mzML)

- **Tool:** ProteoWizard `msconvert` (recommended via Docker)
- **Settings:** 64-bit, vendor peak-picking, zlib compression
- **Script:** `workflows/00_convert_raw_or_wiff_to_mzml.*`

---

### 🧩 Stage 2 — Extraction (mzML → MS1/MS2)

- **Tool:** `xtract_xml.exe` *(external dependency)*
- **Script:** `workflows/01_extract_ms1_ms2_from_mzml.*`

> ⚠️ **Licensing note:** `xtract_xml.exe` is treated as an **external dependency** and is not redistributed unless explicitly permitted.

---

## 🗂️ Project structure

```text
epiprofile-plants/
├── bundles/                # Isolated species-specific runtimes (AT, MP, CR, PP)
├── docs/                   # Architecture, manuals, provenance tiers
├── metadata/               # Audit snapshots and dataset manifests
├── reference/              # Upstream EpiProfile 2.0 baseline
└── workflows/              # Reproducible conversion scripts
```

---

## ⚡ Quickstart (MATLAB)

To ensure deterministic execution, you must isolate the MATLAB path to **exactly one** species bundle.

```matlab
% 1) Clean environment
restoredefaultpath;

% 2) Load a single bundle (example: Arabidopsis)
addpath(genpath("bundles/EpiProfile2.0_AT"));

% 3) Sanity check: MUST return EXACTLY one path
which EpiProfile -all

% 4) Launch GUI
EpiProfile;
```

✅ If `which EpiProfile -all` returns more than one entry, your MATLAB path is contaminated and results become non-deterministic.

---

## 🧾 Provenance & audit model (T1–T4)

All MATLAB functions are tracked using a four-tier provenance model:

- **T1 (Reused):** upstream functions used without changes  
- **T2 (Modified):** upstream code adapted for plant-specific workflows  
- **T3 (New):** new implementations (e.g., `init_histone0_AT.m`)  
- **T4 (Excluded):** intentionally removed / unused modules (e.g., SILAC/C13 modes)

> **Detailed audit records:** `metadata/audit_master.tsv`

---

## ⚠️ Known failure modes

- **Path contamination:** if `which EpiProfile -all` returns more than one entry, execution is non-deterministic.
- **Empty MS1/MS2:** often caused by missing **vendor peak-picking** during `msconvert`.
- **RT window drift:** if quantification fails, widen RT windows in curated layouts.
- **Special characters in paths:** avoid spaces or Unicode in dataset paths  
  (use `C:\data\` instead of `C:\Users\Mi Usuario\Data`).

---

## 📝 Citation & license

### License
Distributed under **GPL-3.0**.

### How to cite
- Cite **EpiProfile_PLANTS** using the repository `CITATION.cff`
- Cite the upstream publication: **Zuo-Fei Yuan et al., EpiProfile 2.0**

---

## 🔧 Optional repo setup helpers

If you want, I can generate:
- `CITATION.cff` (ready for GitHub citation UI)
- a template `metadata/audit_master.tsv` (T1–T4 audit scaffold)
