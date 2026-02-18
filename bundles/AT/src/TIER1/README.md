# TIER1 — Upstream functions (unchanged)

These 72 `.m` files are identical copies from Yuan et al.'s EpiProfile 2.0
basic (commit used as upstream baseline — see `docs/provenance.md`).

## Why they are here

Each species bundle must be self-contained so it can run without external
dependencies on the upstream repository. TIER1 files are copied verbatim;
no lines have been added, removed, or modified.

## Contents by functional block

| Block | Files | Role |
|-------|-------|------|
| H3 quantification | H3_01 through H3_09u | Per-peptide hDP quantification for H3 |
| H4 quantification | H4_01 through H4_04 | Per-peptide hDP quantification for H4 |
| Snapshots | H3_Snapshot, H4_Snapshot | Aggregate per-run results |
| Output | OutputTogether, OutputSinglePTMs | Cohort-level tables |
| Core utilities | Getaamass, calculate_pepmz, find_pair, find_pair_new, get_area, relocateD, GetMods | Shared helpers |
| Reference | check_ref, DrawISOProfile0 | RT reference building |
| I/O | ReadInput, DrawISOProfile2 | Data loading, SILAC runner |

## Provenance audit

Every TIER1 file is registered in `metadata/audit_master.tsv` with tier = T1.
To verify integrity against upstream, compare SHA-256 checksums (see
`docs/provenance.md` Verification section).
