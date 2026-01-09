🌿 EpiProfile_PLANTS

EpiProfile_PLANTS is a plant-oriented, auditable, and reproducible extension of EpiProfile 2.0 (Zuo-Fei Yuan et al., GPL-3.0) for the quantification of:

🧬 hDP — histone-derived peptides (peptide backbones)

🧩 hPF — histone peptideforms (PTM combinatorial states)

🧷 hPTM — site-level histone PTMs (marginal summaries derived from hPFs)

…from bottom-up LC–MS/MS histone datasets generated under anhydride derivatization workflows (BottomMap-like).

🚀 TL;DR (one screen overview)
✅ Gold-standard input

📄 *.mzML

📄 *.MS1

📄 *.MS2

🧯 *.raw (optional dummy placeholder; legacy compatibility only)

✅ Canonical workflow
RAW/WIFF  ->  (Docker + msconvert)  ->  mzML  ->  (xtract_xml.exe)  ->  MS1/MS2  ->  EpiProfile_PLANTS

✅ Deterministic execution

📦 One MATLAB bundle at a time. Always.

✨ What makes this repo different
✅ 1) Gold-standard computational input is NOT RAW/WIFF

EpiProfile_PLANTS standardizes the true computational input as:

📦 Gold-standard input

*.mzML

*.MS1

*.MS2

🧯 Optional legacy compatibility

*.raw → may be a dummy empty placeholder if required by legacy assumptions

Key idea: quantification is driven by mzML + MS1 + MS2.
.raw is not required for quantification.

✅ 2) Fully explicit preprocessing pipeline (Docker + xtract_xml)

The gold-standard artifacts are generated through two explicit steps:

🧱 Stage 1 — Conversion

RAW/WIFF → mzML using ProteoWizard msconvert in Docker

🧩 Stage 2 — Extraction

mzML → MS1/MS2 using xtract_xml.exe (external dependency, long-lived binary)

✅ 3) Deterministic runtime via bundle isolation

EpiProfile_PLANTS is organized into species bundles.

🚫 Invariant: load one bundle at a time in MATLAB path.

✅ Why this matters:

avoids function collisions

guarantees deterministic behavior

makes audits clean and attributable

🧠 Concepts & terminology (repo-wide)
Concept	Meaning	Output
hDP	histone-derived peptide backbone	hDP × sample matrix
hPF	peptideform (PTM combinatorial state)	peptideform-level artifacts
hPTM	site-level PTM derived from hPFs	hPTM × sample matrix

This separation prevents ambiguity and supports both peptide-level and site-level analytics.

🧪 Experimental scope & assumptions

EpiProfile_PLANTS assumes:

✅ bottom-up histone proteomics
✅ anhydride derivatization (BottomMap-like)
✅ DDA-style LC–MS/MS acquisition
✅ curated RT-aware layouts

⚠️ If your chemistry/acquisition differs, results may require catalog/layout adaptation.

⚙️ Preprocessing workflow (gold-standard input generation)
🧱 Stage 1 — RAW/WIFF → mzML (Docker + msconvert)

Conversion is performed using ProteoWizard msconvert inside Docker.

✅ Recommended settings:

64-bit mode

peakPicking vendor

zlib compression

📁 Scripts:

workflows/00_convert_raw_or_wiff_to_mzml.*

🧩 Stage 2 — mzML → MS1/MS2 (xtract_xml.exe)

✅ Generates:

*.MS1

*.MS2

📁 Scripts:

workflows/01_extract_ms1_ms2_from_mzml.*

⚠️ Licensing note:
xtract_xml.exe is an external dependency (not redistributed unless explicitly permitted).

🗂️ Minimal folder recipe (recommended template)
dataset/
├── raw/                          # optional (source only)
│   ├── sample1.raw / sample1.wiff
│   └── ...
├── mzML/
│   ├── sample1.mzML
│   └── ...
├── MS1_MS2/
│   ├── sample1.MS1
│   ├── sample1.MS2
│   └── ...
├── layouts/
│   ├── AT_histone_*.txt
│   └── ...
├── phenodata/
│   └── phenodata.tsv
├── output/
│   ├── epiprofile/
│   └── qc/
└── logs/
    ├── msconvert.log
    └── xtract_xml.log


✅ Rule: your quantification is defined by:

mzML/

MS1_MS2/

layouts/

phenodata/

Everything else is provenance / traceability.

⚡ Quickstart (reproducible, recommended)
0) Ensure gold-standard input exists (per run)

✅ run.mzML
✅ run.MS1
✅ run.MS2
🧯 run.raw (optional placeholder)

1) Run exactly one bundle in MATLAB

Example: Arabidopsis bundle

restoredefaultpath;
addpath(genpath("bundles/EpiProfile2.0_AT"));

which EpiProfile -all
EpiProfile;


✅ Sanity check:

which EpiProfile -all returns exactly one EpiProfile.m

otherwise your MATLAB path is contaminated

🧾 Repository structure (software-grade)
epiprofile-plants/
├── README.md
├── LICENSE
├── CITATION.cff
├── reference/                      # upstream intact baseline (EpiProfile 2.0)
├── bundles/                        # species bundles (isolated runtime)
├── metadata/                       # datasets + audit snapshots
├── docs/                           # manual + architecture + provenance
└── workflows/                      # reproducible scripts and recipes

📦 bundles/

Examples:

bundles/EpiProfile2.0_AT/ — Arabidopsis

bundles/EpiProfile2.0_MP/ — Marchantia

bundles/EpiProfile2.0_CR/ — Chlamydomonas

bundles/EpiProfile2.0_PP/ — Physcomitrella

Each bundle includes:

full MATLAB runtime

species catalogs (init_histone0_*)

species panels + RT-aware layouts

📤 Outputs (high-level)

EpiProfile_PLANTS produces standardized artifacts enabling:

🧬 hDP × sample matrices

🧩 hPF artifacts (peptideforms)

🧷 hPTM × sample matrices

🧾 QC & audit artifacts

📌 See: docs/MANUAL.md

🧾 Provenance & audit model (T1–T4)

All MATLAB functions are classified into four tiers:

T1 — Reused unchanged

T2 — Copied and modified

T3 — Newly implemented

T4 — Excluded intentionally (e.g., SILAC/C13-heavy)

📌 Sources of truth:

docs/tiers.md

docs/provenance.md

metadata/audit_master.tsv

✅ Policy: any code change must update provenance + docs.

✅ Reproducibility checklist (paper/tesis-ready)

Before claiming reproducibility, ensure:

 ✅ Bundle isolation: only 1 bundle in MATLAB path

 ✅ Each run has mzML + MS1 + MS2

 ✅ paras.txt and layouts archived with the run

 ✅ Conversion logs archived (msconvert, xtract_xml)

 ✅ phenodata.tsv pinned to exact run IDs / filenames

 ✅ QC artifacts included (not only matrices)

 ✅ Repo version/commit hash recorded in metadata

⚠️ Known failure modes (save yourself hours)
🧨 1) MATLAB path contamination (multiple bundles loaded)

Symptom: which EpiProfile -all returns >1 entry
Fix:

restoredefaultpath;
addpath(genpath("bundles/EpiProfile2.0_AT"));

🧨 2) xtract_xml.exe generates empty MS1/MS2

Causes: wrong mzML type, missing vendor peak-picking, corrupted conversion
Fix: regenerate mzML with recommended msconvert settings.

🧨 3) Layout mismatch (wrong organism / wrong peptide catalog)

Fix: confirm correct bundle + corresponding layouts.

🧨 4) RT window drift (layout RT windows too strict)

Fix: widen/curate RT windows using representative runs.

🧨 5) Paths with spaces / unicode characters

Fix: move dataset to a simple path:

C:\data\epiprofile\dataset\
/mnt/data/epiprofile/dataset/

🧨 6) Dummy .raw placeholder issues

If a legacy path expects .raw, ensure:

file exists (can be empty)

name matches expected run ID

correct path in paras.txt

🧩 Design principles & invariants (core contract)

These are hard rules. Breaking them means the run is not considered valid/reproducible.

🔒 Invariants

🧱 Bundle isolation

exactly one bundle loaded in MATLAB path

🧬 Quantification input is explicit

mzML + MS1 + MS2 define the run (RAW/WIFF are upstream sources only)

🧾 Traceability is mandatory

preprocessing logs + parameters + layouts are archived

🧩 Species awareness is explicit

catalogs and layouts are organism-specific (no silent cross-species reuse)

🧪 BottomMap scope

supports bottom-up anhydride derivatization workflows; excluded modes are explicit (T4)

🧠 QC philosophy (what “quality” means here)

QC in EpiProfile_PLANTS is not “one plot”. It is a multi-layer audit:

✅ QC layers

🧪 Acquisition-level QC (instrument reports, chromatography stability)

🧱 Conversion QC (msconvert settings + mzML integrity)

🧩 Extraction QC (MS1/MS2 not empty; expected scan counts)

🧬 Peptide-level QC (expected hDP presence, XIC sanity checks)

🧷 PTM-level QC (expected ladders, isotopic patterns, RT alignment)

📊 Matrix-level QC (missingness, zeros, distributions, batch effects)

✅ QC principle

If the pipeline cannot explain why a value is missing/zero, it is not auditable.

⚖️ Licensing & dependencies (matrix)

This repo is designed to keep licensing clean and explicit.

Component	Role	Distributed in repo?	License / Notes
EpiProfile_PLANTS MATLAB code	core quantification	✅ yes	GPL-3.0-only (derived from upstream)
Upstream EpiProfile 2.0 baseline	reference/audit	✅ yes	GPL-3.0 (upstream)
MATLAB	runtime	❌ no	proprietary (institutional license recommended)
Docker	reproducible conversion	❌ no	external dependency
ProteoWizard msconvert	RAW/WIFF → mzML	❌ no	used via Docker, separate distribution
xtract_xml.exe	mzML → MS1/MS2	❌ no (default)	treated as external binary; placed locally by user
Datasets (PRIDE PXDs)	validation	❌ no	referenced in metadata/datasets.tsv

✅ Principle: this repo ships code + documentation, not vendor binaries.
External tools are described and pinned via workflows/logs.

🧾 Audit snapshot (example snippet)

The authoritative audit is stored in:

metadata/audit_master.tsv

A typical record includes:

function	tier	origin	notes
EpiProfile.m	T2	copied+modified	bundle orchestration changes
ReadInput.m	T2	copied+modified	plant-aware parsing / artifacts
init_histone0_AT.m	T3	new	Arabidopsis catalog
get_rts_AT.m	T3	new	RT library for AT
RawToMS1.exe	T1	reused	upstream binary
Extract_C13_*.m	T4	excluded	SILAC/C13-heavy mode not used
✅ Tier policy

T1/T2: must preserve upstream attribution + diff notes

T3: must include standalone docs and rationale

T4: must include explicit reason for exclusion

🧭 Roadmap (high-level)
🎯 Short-term

strengthen QC/audit artifacts per run

finalize stable bundle set (AT + MP + CR + PP)

formalize minimal manifests: RUN_MANIFEST.tsv/json

🧬 Mid-term

cross-species “core peptide panel” bundle (init_histone0_core)

standardized outputs for downstream R workflows (EDA + differential)

📦 Long-term

DOI-tagged releases (Zenodo)

containerized reproducible workflows (where possible)

optional Shiny/interactive QC viewer

📝 Citation
Cite EpiProfile_PLANTS

Use CITATION.cff (GitHub: “Cite this repository”).
Stable releases will have a Zenodo DOI.

Cite upstream EpiProfile 2.0

This software is derived from EpiProfile 2.0.
Please cite the upstream publication when using EpiProfile_PLANTS.

⚖️ License

Distributed under GPL-3.0-only (see LICENSE).
Derived work of EpiProfile 2.0 → upstream licensing applies.

🤝 Contributing (strict rules)

Contributions welcome under these constraints:

🧱 Bundle isolation is sacred

🧾 Update provenance (metadata/audit_master.tsv)

📚 Document new/modified functions (docs/functions/ EN/ES)

🆘 Support

When opening an issue, include:

OS + MATLAB version

bundle used

preprocessing logs (msconvert, xtract_xml)

input artifacts (mzML + MS1 + MS2) or minimal subset

paras.txt and layout files used

✅ Final note

This README defines the contract of EpiProfile_PLANTS:
explicit inputs, isolated bundles, auditable provenance, and reproducible artifacts.
