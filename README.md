# EpiProfile_PLANTS

MATLAB extension of EpiProfile 2.0 for **plant histone proteomics**.

EpiProfile_PLANTS provides:

- Histone peptide catalogs and layouts for:
  - *Arabidopsis thaliana* (AT)
  - *Marchantia polymorpha* (MP)
  - *Chlamydomonas reinhardtii* (CR)
- Adapted MATLAB functions for plant-specific histone peptides and variants.
- Utilities for QC and basic reporting of histone PTM quantification.

> This repository focuses on the **MATLAB code and histone catalogs**.
> For the end-to-end workflow (Docker + Snakemake from WIFF to EpiProfile_PLANTS),
> see the companion repository `epiprofile-plants-workflow`.

---

## 1. Installation

EpiProfile_PLANTS currently runs in **MATLAB** (tested versions to be documented).

1. Clone or download this repository.
2. Add the `src/` folder to the MATLAB path.
3. Ensure that EpiProfile 2.0 binaries and required executables (e.g. `RawToMS1.exe`, `xtract_xml`) are available in your system.

Instructions for using MATLAB Runtime / standalone packages will be added in future releases.

---

## 2. Repository structure (planned)

- `src/` – MATLAB source code.
  - `core/` – functions inherited or adapted from EpiProfile 2.0.
  - `plants/` – new plant-specific functions and layouts.
  - `qc/` – QC and auditing utilities.
- `layouts/` – text/MAT files describing histone panels (H3/H4, etc.).
- `catalogs/` – `init_histone0_*` and histone peptide catalogs for AT, MP, CR.
- `examples/` – small example datasets and scripts.
- `docs/` – user and developer documentation.

The exact structure may evolve as the project is cleaned up and documented.

---

## 3. Usage (high level)

Typical usage within a WIFF→mzML→MS1/MS2 pipeline:

1. Convert vendor files to `mzML` (see `epiprofile-plants-workflow`).
2. Generate `MS1`/`MS2` text files with `Raw2MS` / `xtract_xml`.
3. Configure the appropriate **layout** and **init_histone0_*`** file for the species.
4. Run `EpiProfile_PLANTS` main script (to be documented) to obtain histone PTM areas.
5. Export peptide/PTM matrices for downstream analysis (R, Python).

Detailed step-by-step examples will be provided in the `examples/` folder.

---

## 4. Relation to EpiProfile 2.0

EpiProfile_PLANTS is a **derivative work** based on EpiProfile 2.0:

- Some MATLAB files are reused without changes.
- Others are adapted for plant histone peptides.
- New functions and catalogs are added specifically for plants.

A full provenance map (copied / modified / new / not used) is being prepared
as part of the associated PhD thesis and software paper.

---

## 5. Citation and license

Citation (preliminary):

> González de Lena P. *EpiProfile_PLANTS: reproducible histone PTM quantification in plant systems* (manuscript in preparation).

License: to be finalised (GPL-3.0 is the planned default for the combined work).

---

## 6. Related repositories

- Workflow (Docker + Snakemake):  
  <https://github.com/biopelayo/epiprofile-plants-workflow>

- Thesis repository (H3K79 in Arabidopsis):  
  <https://github.com/biopelayo/thesis-h3k79-arabidopsis>
