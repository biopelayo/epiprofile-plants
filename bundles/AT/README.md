# AT — *Arabidopsis thaliana* bundle

First and reference species bundle for EpiProfile_PLANTS.

## Status

**Active.** Runs end to end on *Arabidopsis* data (validated 2026-08-18 on 5 SCIEX ZenoTOF
7600 DDA runs, and through the 2026 thesis runs on PXD046034, PXD046788 and PXD014739).

## Source data

| Item | Value |
|------|-------|
| Organism | *Arabidopsis thaliana* (Col-0) |
| Histones quantified | H3 (H3.1/H3.3), H4, H2A (canonical, H2A.Z, H2A.X, H2A.W), H2B, H1 (H1.1/H1.2) |
| Derivatization | Propionylation + trypsin (Arg-C-like specificity) |
| Instrument modes | DDA and DIA supported |

## Key sequence differences vs. human

- H3 27-40 carries two plant-specific isoforms (H3.3: KSAPTTGGVKKPHR; H3.1-like:
  KSAPATGGVKKPHR / QSAPATGGVKKPHR) handled by TIER2/TIER3 modules; other substitutions
  across H3 3-8, 9-17, 18-26, 53-63, 73-83, 117-128 by dedicated TIER3 modules.
- H4 68-78 is `DAVTYTEHAR` (R77, no internal K) and H4 24-102 carries I61: the upstream
  human modules were adapted (`H4_05_68_78`, `H4_07u_24_102`, `H4_Snapshot`).
- H2A/H2B/H1: AT-specific modules `HH2A_AT01..08`, `HH2B_AT01..07`, `HH1_AT01..04`
  (peptides verified against UniProt, see `../../assets/sequences/`).

## Layout

```
src/TIER1..TIER4/   MATLAB code by tier (see src/README.md)
metadata/           module manifest, nucleosome master table, RT tables, per-bundle audit
scripts/            provenance / catalogue tools (Python; not needed to run the pipeline)
```

## Running

From the repository root, in MATLAB:

```matlab
restoredefaultpath;
addpath(genpath('bundles/AT/src'));
which EpiProfile -all            % must be exactly one
EpiProfile('C:\data\MS1_MS2\paras.txt');   % see ../../paras.example.txt
```

Defaults live in `src/TIER2/check_otherparas.m` (`def_ptol = 10` ppm, `soutput = '31'`,
`nfigure = 0`, `ndebug = 0`). See `docs/MANUAL.md` for the data layout, the RT reference
system and the output contract.
