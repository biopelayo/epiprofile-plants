# bundles/

Each subdirectory is a **standalone species bundle** for EpiProfile_PLANTS.

A bundle contains all the MATLAB functions needed to quantify histone PTMs from
LC-MS/MS data for one plant species.

## Directory layout

```
bundles/
  AT/src/            Arabidopsis thaliana  (active — validated with PXD014739)
    TIER1/           Upstream EpiProfile 2.0 functions (unchanged)
    TIER2/           Functions modified for plant histones
    TIER3/           New functions for plant sequence variants (MSA-derived)
  MP/src/            Marchantia polymorpha  (assembled — pending validation)
    TIER1/           Upstream EpiProfile 2.0 functions (unchanged)
    TIER2/           Functions modified for MP (KSAPSTGGVKKPHR, KVFR)
    TIER3/           AT-derived + MP-specific modules
  CR/src/            Chlamydomonas reinhardtii  (assembled — pending validation)
    TIER1/           Upstream EpiProfile 2.0 functions (unchanged)
    TIER2/           Functions modified for CR (KTPATGGVKKPHR)
    TIER3/           AT-derived + CR-specific modules
```

## How bundles are built

1. **TIER1** — Copy every function from upstream EpiProfile 2.0 basic that is
   reused without changes (72 files for AT).
2. **TIER2** — Copy functions that need modification (e.g. init_histone0,
   DrawISOProfile1, get_rts, draw_layout) and adapt them to the plant peptide
   catalog and robustness requirements.
3. **TIER3** — Write new quantification modules for peptide sequences that
   differ from the human reference, identified through Multiple Sequence
   Alignment of histone protein FASTA files.

See `docs/tiers.md` for full tier definitions and `docs/provenance.md` for
file-level provenance tracking.
