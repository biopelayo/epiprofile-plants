# Tiers (T1–T4) — Function provenance & modification level

This repository audits and documents the origin and modification level of every file/function.

## T1 — Identical (reused)
- File is identical to upstream (EpiProfile 2.0 basic).
- No modifications (byte-to-byte, or functionally identical).
- Example: vendor binaries reused as-is, or unchanged MATLAB functions.

## T2 — Copied & modified
- File originates from upstream but has been edited.
- Changes may include: parameters, logic, bug fixes, new features, refactors, new outputsinia.

## T3 — New (created in PLANTS family)
- File does not exist in upstream.
- New functionality specific to plants/species, new panels, new init files, helpers, etc.

## T4 — Upstream-only / not ported
- File exists in upstream but is not used/ported in the PLANTS family bundles.
- Typically: SILAC/C13 modules or organism blocks not needed for our scope.

## Notes
- Each audited file must have: Tier, short description, and provenance notes.
- The audit is performed by comparing: `upstream/epiprofile2.0_basic/` vs `family/EpiProfile2.0_XX/`.
