# AT/src — Source code

All MATLAB `.m` files for the *Arabidopsis thaliana* bundle live here, organised by
tier (see `docs/tiers.md`). Counts from `tools/tier_audit.py` against EpiProfile 2.0
basic (2026-08-18):

| Directory | Tier | Count | Description |
|-----------|------|------:|-------------|
| `TIER1/`  | T1   |  46 | Upstream EpiProfile 2.0 basic, unchanged (normalised for line endings) |
| `TIER2/`  | T2   |  30 | Upstream files modified for plants and hardened (guards, AT sequences, ordering, thresholds) |
| `TIER3/`  | T3   |  32 | New in the PLANTS family: H3 sequence-variant modules, AT H2A/H2B/H1 modules and their Snapshots, `HH2A_01u_1_7`, `run_calibrate` |
| `TIER4/`  | T4   | 119 | Present but **not invoked** by `DrawISOProfile1.m`: 5 unused helpers, 5 legacy upstream H2A/H2B modules, the DRAFT `H3_05b_41_49.m`, and the 108 generated H1/H2A/H2B catalogue modules (`HH1_*oA*`, `HH2A_*oA*`, `HH2B_*oA*`, `H2A_Snapshot`, `H2B_Snapshot`; `rt_ref = 0`, uncalibrated, described in `../metadata/`) |

T4 takes precedence over T1-T3, so the four folders add up to the bundle (227 files).
No basename appears twice, so `which <name> -all` returns one hit for every function once
`addpath(genpath('bundles/AT/src'))` has run.

## Entry point

`TIER2/EpiProfile.m` reads `paras.txt`, builds the run list (`check_otherparas.m`), parses
MS1/MS2 into `.mat` caches, builds the RT reference (`DrawISOProfile0`) and dispatches the
quantification loop:

```
EpiProfile.m  ->  DrawISOProfile0  (RT reference from the runs, histone_layouts/0_ref_info.mat)
              ->  DrawISOProfile1  (per-run modules: H3, H4, H2A, H2B, H1 + Snapshots)
              ->  OutputTogether (histone_ratios.xls) / OutputSinglePTMs / OutputFigures (nfigure=1 only)
```

`DrawISOProfile1.m` selects modules with `soutput` (default `'31'` in `check_otherparas.m`:
all H3/H4 panels + H1/H2A/H2B). The H2A/H2B/H1 block runs the AT-specific modules
(`HH2A_AT01..08`, `HH2B_AT01..07`, `HH1_AT01..04`) plus two conserved upstream ones
(`HH2A_06m1_72_77`, `HH2BMo_03u_1_100`).

## Provenance of this tree

This tree is the AT bundle used for the thesis runs (`revision_peptidos/src_AT` +
`run_PXD046034`, June 2026), including the 2026-05-24 code-review fixes and the AT
sequences for H4 68-78 (`DAVTYTEHAR`) and 24-102 (I61). It replaces the tree published
until commit `eb89002`, which crashed on real data (no empty-window guards in
`get_histone12` and siblings) and still quantified the human H4 sequences. Validated on
2026-08-18: 5 AT runs, exit 0, 14.5 min, 306-line `histone_ratios.xls`.

## Adding the path

```matlab
restoredefaultpath;
addpath(genpath('bundles/AT/src'));   % TIER1..TIER4
which EpiProfile -all                 % exactly one hit
```
