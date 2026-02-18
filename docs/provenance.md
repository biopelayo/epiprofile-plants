# Provenance & scope

## Upstream reference
- The upstream EpiProfile 2.0 *basic* distribution (Yuan et al., 2018) is the provenance baseline.
- TIER1 files in each bundle are copies of upstream functions, used unchanged.
- Audit comparisons use the upstream commit/release as reference point.
- The upstream is **not redistributed** in this repository (GPL-3.0 licensing constraints). Instead, each TIER1 file is documented in `metadata/audit_master.tsv` with its origin.

## PLANTS family bundles
- `bundles/<SPECIES_CODE>/src/` contains **standalone executable bundles**, one per species.
- Each bundle is organized in three tiers under `src/`:
  - `TIER1/` — upstream functions reused unchanged (T1)
  - `TIER2/` — upstream functions modified for plant workflows (T2)
  - `TIER3/` — new functions specific to plant sequence variants (T3)
- Each bundle must run independently without depending on other species folders.

## Species codes
- AT: *Arabidopsis thaliana* (active)
- MP: *Marchantia polymorpha* (planned)
- CR: *Chlamydomonas reinhardtii* (planned)
- PP: *Physcomitrella patens* (if/when included)

Current status is tracked in `metadata/species.tsv`.

## Audit principle
For every file/function we classify:
- Tier (T1-T4)
- Short description
- Provenance note (where it came from and why changes were introduced)

Source of truth: `metadata/audit_master.tsv`

Audit comparisons are primarily:
- A (upstream): the original EpiProfile 2.0 basic function
- C (bundle): `bundles/<XX>/src/TIER<N>/<file>`

## T4 — Excluded modules
Functions that exist in upstream but are **not ported** to the PLANTS family:
- SILAC modules (`Extract_SILAC.m`, `Extract_SILAC_1.m`)
- C13/N15/13CD3 runners (`DrawISOProfile3-5.m`)
- Heavy amino acid masses (`GetaamassH.m`)

These are documented as T4 entries in `audit_master.tsv`.
