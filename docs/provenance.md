# Provenance & scope

## Upstream reference
- `upstream/epiprofile2.0_basic/` contains an **exact**, untouched copy of the EpiProfile 2.0 *basic* distribution.
- This folder is used strictly as a provenance reference for audits.
- No edits should be made inside `upstream/` after import.

## PLANTS family bundles
- `family/EpiProfile2.0_XX/` contains **standalone executable bundles**, one per species.
- Each bundle must run independently without depending on other species folders.

## Species codes (current plan)
- AT: *Arabidopsis thaliana*
- MP: *Marchantia polymorpha*
- CR: *Chlamydomonas reinhardtii*
- PP: *Physcomitrella patens* (if/when included)

## Audit principle
For every file/function we classify:
- Tier (T1–T4)
- Short description (EN/ES)
- Provenance note (where it came from and why changes were introduced)

Audit comparisons are primarily:
- A: `upstream/epiprofile2.0_basic/<file>`
- C: `family/EpiProfile2.0_XX/<file>`
