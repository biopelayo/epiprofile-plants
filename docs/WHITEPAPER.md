# EpiProfile_PLANTS — Methodological Whitepaper

**Adapting EpiProfile 2.0 for plant histone PTM quantification: a systematic,
MSA-driven approach**

---

## 1. Motivation

Histone post-translational modification (PTM) profiling by mass spectrometry
is well-established in human and mouse chromatin biology, largely thanks to
dedicated quantification tools such as EpiProfile (Yuan et al., 2015; 2018).
However, plant histones carry amino acid substitutions relative to the human
reference that make the original peptide panels inapplicable without
modification. Specific residues that serve as PTM sites or propionylation
targets differ between species, and new peptide sequences appear that have no
counterpart in the human panel.

EpiProfile_PLANTS addresses this gap by extending EpiProfile 2.0 basic to
support plant species through a tiered, auditable adaptation process. The
methodology is designed so that (a) every change relative to the upstream
codebase is classified and traceable, (b) new species can be added
systematically by repeating the same workflow, and (c) the resulting bundles
are standalone MATLAB packages that a reviewer or collaborator can run
end-to-end.

---

## 2. Upstream baseline: EpiProfile 2.0 basic

The starting point is the EpiProfile 2.0 basic codebase by Yuan et al.
(publicly available on GitHub). This MATLAB pipeline:

1. Loads MS1/MS2 extracted data from LC-MS/MS runs.
2. Uses a library of histone-derived peptide (hDP) "layouts" — fixed
   combinations of sequence, charge, and modification mass — to define
   XIC extraction targets.
3. Builds a retention-time (RT) reference from the top-N most intense runs
   (via `DrawISOProfile0`), storing it in `0_ref_info.mat`.
4. Quantifies each peptide across all runs by extracting ion chromatograms
   within mass-tolerance and RT windows (`DrawISOProfile1`).
5. Aggregates results into cohort-level ratio tables (`OutputTogether`,
   `OutputSinglePTMs`) and optional QC figures (`OutputFigures`).

The human peptide catalog is hard-coded in `init_histone0.m`, which registers
each expected peptide with its sequence, m/z values, RT hints, and the
quantification module responsible for it.

---

## 3. The tier system (T1-T4)

To maintain full traceability between the plant adaptation and the upstream
codebase, every file in a species bundle is assigned a provenance tier:

| Tier | Meaning | Example |
|------|---------|---------|
| **T1** | Identical to upstream — copied without modification | `Getaamass.m`, `find_pair.m`, `H3_01_3_8.m` |
| **T2** | Copied from upstream and modified | `init_histone0.m`, `get_rts.m`, `DrawISOProfile1.m` |
| **T3** | New — created specifically for plants | `H3_11_3_8.m`, `H3_14_27_40.m`, `H3_Snapshot.m` |
| **T4** | Upstream-only — not ported (out of scope) | `Extract_SILAC.m`, `GetaamassH.m`, `DrawISOProfile3.m` |

This classification is recorded in `metadata/audit_master.tsv` and described
in `docs/tiers.md`. The tier assignment is the primary mechanism by which a
reviewer can assess how much of the pipeline is novel vs. reused.

---

## 4. Adaptation methodology (step by step)

### Step 0 — Sequence retrieval and Multiple Sequence Alignment (MSA)

For each target species, retrieve the full protein sequences for histone H3
(all variants: H3.1, H3.3, and any species-specific isoforms) and H4 from
UniProt or Phytozome.

Perform a Multiple Sequence Alignment using ClustalW, MUSCLE, or equivalent,
aligning each plant histone against the human H3.1 reference
(UniProt: P68431) and H4 reference (UniProt: P62805).

The MSA reveals:
- **Conserved regions**: peptides identical to human after propionylation +
  trypsin digestion. These are handled by existing T1 modules.
- **Substituted regions**: peptides where one or more residues differ. These
  require new T3 modules or T2 modifications.
- **Insertions/deletions**: rare for core histones but checked systematically.

### Step 1 — In-silico digestion and peptide catalog

Apply the propionylation + trypsin (Arg-C-like) digestion rules to each
aligned histone sequence:

1. Propionylate all unmodified and monomethylated lysines.
2. Cleave at the C-terminal side of arginine (R) residues.
3. Predict the resulting tryptic peptides and their masses.

For each peptide region (e.g., residues 3-8, 9-17, 27-40), compare the
plant-specific sequence(s) against the human reference:

- If identical: the upstream quantification module (T1) applies directly.
- If different: a new module is needed (T3), or the existing module must be
  modified (T2) to handle the variant mass ladder.

Record the full peptide catalog in `init_histone0.m`, including:
- `new_seq`: propionylated peptide sequence
- `new_mz`: monoisotopic m/z at the expected charge state
- `seq_godel`: deterministic hash for RT reference lookup, calculated as
  `sum((new_seq-'0'+49).*log(2:1+length(new_seq)))`

### Step 2 — Writing new quantification modules (TIER3)

Each T3 module follows the upstream H3_XX pattern. The anatomy of a module:

```matlab
function H3_XX_start_end(MS1_index, MS1_peaks, MS2_index, MS2_peaks, ...
                         ptol, cur_outpath, special)
```

Inside the module:

1. **Define the peptide panel**: list all peptideforms (hPFs) for this
   region — unmodified + every combination of PTMs (me1, me2, me3, ac,
   etc.) applicable to the sequence.

2. **Calculate m/z for each hPF**: using `Getaamass.m` and
   `calculate_pepmz.m` to compute the monoisotopic mass, then dividing
   by charge state.

3. **Retrieve RT reference**: call `get_rts.m` to look up the expected
   RT from `0_ref_info.mat` (if available) using `seq_godel`.

4. **Extract XICs from MS1**: for each hPF, extract the ion chromatogram
   within +/- `ptol` ppm around the expected m/z.

5. **Quantify**: call `get_histone0.m` (for unmodified) or the
   appropriate `get_histoneXX.m` helper to integrate areas and pick
   the best RT window.

6. **Relocate (if needed)**: for modified forms where signal may be
   weak, use `relocateD.m` to shift the search window relative to the
   unmodified peptide's RT.

7. **Validate with MS2**: use MS2 fragmentation data to confirm peptide
   identity when available (DDA mode).

8. **Write results**: save per-hPF areas, ratios, and RT to a `.mat`
   file in `cur_outpath/detail/`.

For plant-specific modules, the key difference is the mass ladder: each
substituted amino acid changes the peptide mass, the isotope envelope, and
potentially the fragmentation pattern. The module must encode the correct
mass and any new PTM combinations that become possible (e.g., a K->R
substitution removes an acetylation site).

### Step 3 — Modifying existing modules (TIER2)

Some upstream functions require changes beyond what a new T3 module provides:

- **`init_histone0.m`**: Expand the peptide catalog to include all
  plant-specific sequences. This is the central registry where every hDP
  is defined with its sequence, m/z, charge, and module assignment.

- **`DrawISOProfile1.m`**: The main quantification runner. Must call all
  T1 modules, all T2 modules, and all T3 modules in sequence. In the
  upstream codebase, this function was missing from the basic distribution
  and had to be reconstructed.

- **`get_rts.m`**: Add robustness guards for cases where the RT reference
  is empty or the peptide is not found (common when first running on a
  new species without prior RT data).

- **`get_histone0.m`**: Add gradient-length awareness and guards against
  empty MS1/MS2 arrays that can occur with plant samples of lower
  abundance than human cell-line extracts.

- **`draw_layout.m`**: Fix XIC plotting for edge cases where no signal is
  found (prevents MATLAB crashes on empty arrays).

### Step 4 — Assembling the species bundle

A complete bundle for species XX is a self-contained MATLAB directory:

```
bundles/XX/src/
  TIER1/    72 files (unchanged upstream)
  TIER2/    ~10 files (modified)
  TIER3/    variable (one per variant peptide region + snapshot aggregator)
```

The user runs the pipeline by calling `EpiProfile.m` from `TIER2/`, which:
1. Calls `init_histone0.m` to register the species-specific peptide catalog.
2. Calls `DrawISOProfile0` to build the RT reference.
3. Calls `DrawISOProfile1` to run quantification.
4. Calls `OutputTogether` + `OutputSinglePTMs` + `OutputFigures`.

### Step 5 — Validation

Each bundle is validated against a reference LC-MS/MS dataset:

- **AT (Arabidopsis)**: PRIDE PXD014739 (Lochmanová et al. 2019, *IJMS* 20:5093) — acid-extracted histones, DDA
  acquisition on LTQ Orbitrap Elite. Study of genetic/chemical HDAC downregulation in *A. thaliana*.
- **MP (Marchantia)** and **CR (Chlamydomonas)**: pending validation datasets.

Validation checks:
1. All expected peptides are detected (non-zero area in the ratio table).
2. RT consistency across replicates (CV < 5%).
3. Ratio distributions match biological expectations (e.g., H3K27me3
   enrichment in heterochromatin fractions).
4. XIC traces show clean peaks with correct isotope envelopes.

---

## 5. The AT bundle in detail

The Arabidopsis thaliana bundle is the reference implementation. Below is the
complete inventory of plant-specific changes.

### 5.1 TIER3 modules — sequence variants from MSA

| Module | Residues | Human reference | AT variants | Substitutions |
|--------|----------|----------------|-------------|---------------|
| H3_11_3_8 | 3-8 | TKQTAR | TKQSAR, SNQTAR | T->S (pos 5), TK->SN (pos 3-4) |
| H3_12_9_17 | 9-17 | KSTGGKAPR | KSTGGKGPR, KSHGGKAPR, ISTGGKAPR | A->G (pos 15), A->H (pos 12), K->I (pos 9) |
| H3_13_18_26 | 18-26 | KQLATKAAR | KELATKAAR, TLLATKAAR, KQLAPKAAR | Q->E (pos 19), K->T+Q->L (pos 18-19), T->P (pos 22) |
| H3_14_27_40 | 27-40 | KSAPATGGVKKPHR | KSAPTTGGVKKPHR, QSAPATGGVKKPHR | A->T (pos 31, H3.3), K->Q (pos 27, H3.1-like) |
| H3_16_53_63 | 53-63 | KYQKSTELLIR | KYQKSTELLNR | I->N (pos 62) |
| H3_17_73_83 | 73-83 | EIAQDFKTDLR | EIAQDYKTDLR | F->Y (pos 78) |
| H3_18_117_128 | 117-128 | VTIMPKDIQLAR | VTIMPKDVQLAR, VTIMPKEIQLAR | I->V (pos 123), D->E (pos 123) |

Additionally, `H3_06_53_63.m` provides the full PTM panel with relocate
logic for the KYQKSTELLIR sequence region (complementing the variant module
H3_16_53_63).

### 5.2 TIER2 modifications summary

| File | Change type | Rationale |
|------|-------------|-----------|
| init_histone0.m | Catalog expansion | 56 peptide entries including all AT variants |
| DrawISOProfile1.m | New (gap-fill) | Main runner calling T1+T2+T3 modules |
| get_rts.m | Robustness | Guards for empty RT reference |
| get_histone0.m | Robustness | Gradient-aware extraction, empty-array guards |
| get_histone11.m | Robustness | GetTopBottom11 variant |
| draw_layout.m | Robustness | XIC plotting fixes |
| H3_04v3_27_40.m | Variant | KSAPTTGGVKKPHR (H3.3-specific) |
| OutputFigures.m | Adaptation | QC figures for plant PTM panel |
| check_otherparas.m | Defaults | ptol=10, ndebug=0 |
| EpiProfile.m | Entry point | Dispatch logic (minimal changes) |

---

## 6. Extending to new species

The methodology for adding a new species (e.g., MP or CR) follows the same
5-step process:

1. **MSA**: Align the new species' histone sequences against human + AT.
2. **Digest**: Predict propionylation + trypsin peptides.
3. **Compare**: Identify which peptide regions match AT (reuse T3 modules),
   which match human (reuse T1 modules), and which are unique (write new
   T3 modules).
4. **Fork**: Copy the AT bundle as a starting point. Update
   `init_histone0.m` with the new peptide catalog. Add/replace T3 modules
   as needed.
5. **Validate**: Run against a reference dataset for the new species.

Known differences for planned species:

| Species | Key variant (H3 27-40) | Other known changes |
|---------|----------------------|---------------------|
| *M. polymorpha* (MP) | KSAPSTGGVKKPHR | H4 20-23: KVFR (L->F) |
| *C. reinhardtii* (CR) | KTPATGGVKKPHR | Pending full MSA |

---

## 7. Data model

The pipeline produces three levels of data:

| Level | Abbreviation | Description |
|-------|-------------|-------------|
| Histone-Derived Peptide | hDP | Raw XIC area for a specific peptide + charge + isotope |
| Histone Peptideform | hPF | A specific modification state of a peptide (e.g., H3K27me3 on KSAPATGGVKKPHR) |
| Histone PTM summary | hPTM | Site-level PTM abundance (ratio), aggregated across peptideforms |

The final cohort export (`histone_ratios.xls`, actually TSV) contains hPF-level
ratios, areas, and RTs in a two-header-row format. See `docs/MANUAL.md`
section 15 for the exact output contract.

---

## 8. Retention-time reference system

The RT reference mechanism is central to EpiProfile's accuracy:

1. **`DrawISOProfile0`** processes the top-3 most intense runs to establish
   RT anchors for each peptide, stored in `0_ref_info.mat`.
2. Each anchor is keyed by `seq_godel` — a deterministic hash of the peptide
   sequence: `sum((new_seq-'0'+49).*log(2:1+length(new_seq)))`.
3. During quantification, `get_rts.m` retrieves the expected RT for each
   peptide and defines the extraction window.
4. The `relocateD.m` function shifts the search window for modified peptides
   relative to their unmodified counterpart's RT.

The reference is "global by folder" (`raw_path`), meaning all runs in the same
folder share the same RT anchors. This requires that runs within a folder share
compatible LC conditions (same column, gradient, instrument).

---

## 9. Quality control

The pipeline includes several QC layers:

1. **XIC visual inspection**: `draw_layout.m` produces per-peptide XIC plots
   stored in `detail/` for manual review.
2. **Snapshot aggregation**: `H3_Snapshot.m` and `H4_Snapshot.m` compile
   per-run summaries that can be inspected before cohort aggregation.
3. **OutputFigures**: Generates cohort-level QC plots (bar charts, heatmaps)
   for rapid assessment of data quality.
4. **Survey template**: `docs/SURVEY.md` provides a structured checklist for
   documenting each analysis run's parameters and outcomes.

---

## 10. Repository structure

```
epiprofile-plants-main/
  bundles/
    AT/src/TIER1/       72 files — upstream unchanged (T1)
    AT/src/TIER2/       10 files — modified for plants (T2)
    AT/src/TIER3/        9 files — new for AT variants (T3)
    MP/                  planned
    CR/                  planned
  docs/
    MANUAL.md            end-to-end usage guide
    SURVEY.md            structured run-report template
    tiers.md             tier definitions
    provenance.md        file-level provenance tracking
    WHITEPAPER.md        this document
  metadata/
    audit_master.tsv     per-file tier audit (97 entries)
    species.tsv          species bundle registry
  workflows/
    PXD014739/           reference validation dataset (AT)
  CITATION.cff           citation metadata
  LICENSE                GPL-3.0
```

---

## 11. References

- Yuan, Z.F. et al. (2015). EpiProfile quantifies histone peptides with
  modifications by extracting retention time and intensity in
  high-resolution mass spectra. *Mol. Cell. Proteomics* 14, 1696-1707.
- Yuan, Z.F. et al. (2018). EpiProfile 2.0: A Computational Platform for
  Processing Epi-Proteomics Mass Spectrometry Data. *J. Proteome Res.* 17,
  2533-2541.
- EpiProfile 2.0 basic source code: upstream GitHub repository.

---

*Document version: 1.0 — Generated 2026-02-18*
