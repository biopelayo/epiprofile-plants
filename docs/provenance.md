# Provenance & scope

## Upstream reference
- The upstream EpiProfile 2.0 *basic* distribution (Yuan et al., 2018) is the provenance baseline.
- TIER1 files in each bundle are copies of upstream functions, used unchanged.
- Audit comparisons use the upstream commit/release as reference point.
- The upstream is **not redistributed** in this repository (GPL-3.0 licensing constraints). Instead, each TIER1 file is documented in `metadata/audit_master.tsv` with its origin.

## PLANTS family bundles
- `bundles/<SPECIES_CODE>/src/` contains **standalone executable bundles**, one per species.
- Each bundle is organized in four tiers under `src/` (one folder per tier, so the folder
  layout and the T1-T4 arithmetic of `docs/tiers.md` coincide):
  - `TIER1/` — upstream functions reused unchanged (T1)
  - `TIER2/` — upstream functions modified for plant workflows (T2)
  - `TIER3/` — new functions specific to plants (T3)
  - `TIER4/` — present but not invoked (T4; precedence over T1-T3). In AT this holds the
    108 generated H1/H2A/H2B catalogue modules (`rt_ref = 0`), a few unused helpers and
    legacy upstream modules, and drafts.
- Each bundle must run independently without depending on other species folders.

## Species codes
- AT: *Arabidopsis thaliana* (active, validated on data)
- MP: *Marchantia polymorpha* (assembled, not yet validated on data)
- CR: *Chlamydomonas reinhardtii* (assembled, not yet validated on data)
- PP: *Physcomitrella patens* (if/when included)

Current status is tracked in `metadata/species.tsv`.

## Audit principle
For every file/function we classify:
- Tier (T1-T4)
- Short description
- Provenance note (where it came from and why changes were introduced)

Source of truth: `metadata/audit_master.tsv`, regenerated with
`python tools/build_audit_master.py --upstream <EpiProfile2.0_1Basic>` (which calls
`tools/tier_audit.py` per bundle: normalised-content comparison against upstream for T1/T2/T3,
invocation search ignoring comments and strings for T4). Per-bundle detail:
`bundles/<XX>/metadata/tier_audit_<XX>.tsv`.

Audit comparisons are primarily:
- A (upstream): the original EpiProfile 2.0 basic function
- C (bundle): `bundles/<XX>/src/TIER<N>/<file>`

## T4 — Present but not invoked
Functions that **ship in the bundle** but no other file calls, so `DrawISOProfile1.m` never reaches
them. They are kept as reference or for future development — for example upstream quantification
modules the plant workflow stopped calling (`HH2A_07v_1_88.m`, `HH2B_02v_1_29.m`) and helpers
superseded by a plant-specific variant (`check_layout.m`, `get_main_ch.m`, `output_histone2.m`).

## Upstream functions that were not ported
These exist in upstream but were **left out** of the bundle, so they are **not a tier** — they are
simply out of scope:
- SILAC modules (`Extract_SILAC.m`, `Extract_SILAC_1.m`)
- C13/N15/13CD3 runners (`DrawISOProfile3-5.m`)
- Heavy amino acid masses (`GetaamassH.m`)

These six are tagged **`NOT_PORTED`** in `audit_master.tsv` (both copies: `metadata/` and
`bundles/AT/metadata/`). They used to be tagged `T4` under the previous definition; re-tagged
2026-07-21 so that `T4` keeps its new meaning of "present but not invoked". Backups:
`.bak_tier_20260721`.
