# Tiers (T1–T4) — Function provenance & invocation status

This repository audits and documents the origin, modification level and invocation status of every
file/function in the PLANTS bundles.

Two independent questions are being answered:

1. **Where does this file come from?** → T1 / T2 / T3 (provenance).
2. **Does the pipeline actually call it?** → T4 (invocation).

**Precedence rule: T4 wins.** A file that ships in the bundle but is never invoked is classified as
**T4**, regardless of its provenance. This keeps the four tiers mutually exclusive, so that
`T1 + T2 + T3 + T4` equals the total number of `.m` files in the bundle.

---

## T1 — Identical (reused)

- File is identical to upstream (EpiProfile 2.0 basic).
- No modifications (byte-to-byte, or functionally identical).
- Example: unchanged MATLAB helper functions.

## T2 — Copied & modified

- File originates from upstream but has been edited.
- Changes may include: parameters, logic, bug fixes, new features, refactors, new outputs.
- Example: `get_histone13.m` and `get_histone22.m`, which add empty-vector guards
  (`if isempty(p); return; end;`) that upstream does not have.

## T3 — New (created in PLANTS family)

- File does not exist in upstream.
- New functionality specific to plants/species: new quantification modules, new init files, helpers.

## T4 — Present but not invoked

- File **ships in the bundle** but is never called by any other file in it — it does not appear
  anywhere outside its own definition, so `DrawISOProfile1.m` never reaches it.
- Kept as reference or for future development.
- Typically: upstream modules the plant workflow stopped calling (for example
  `HH2A_07v_1_88.m`, `HH2B_02v_1_29.m`), plus helpers superseded by a plant-specific variant.

> **⚠ Changed 2026-07-21.** T4 used to mean *"upstream-only / not ported"* — files that exist upstream
> but were **left out** of the bundle. That definition contradicted the thesis (Chapter 2, §2.2.4),
> and it also made the tiers unusable as a breakdown of the bundle: files that are not in the bundle
> cannot be a share of it, which is how percentages ended up being computed over a total larger than
> the bundle itself. T4 now means *present but not invoked*, matching the thesis.
>
> Files that exist upstream and were **not ported** are simply **out of the bundle**. They are not a
> tier; describe them as "upstream, not ported" if they need to be mentioned.

---

## Notes

- Each audited file must have: Tier, short description, and provenance notes.
- The provenance audit is performed by comparing `upstream/epiprofile2.0_basic/` against the bundle,
  hashing file contents **normalised for whitespace and line endings** (otherwise line-ending
  differences alone show up as modifications).
- **State the reference version explicitly** in any reported figure. The classification is highly
  sensitive to it: against EpiProfile **2.0 basic** the AT bundle splits roughly 60/26/33, while
  against 2.1, 2.2 or 2.2_AT it collapses to about 2/54/63, because those releases rewrote nearly
  every file. 2.0 basic is the correct reference — it is the one the bundle derives from.
- The invocation audit searches each file's base name across every other file in the bundle,
  ignoring comments and string literals.
- The verified per-tier counts for the bundle used in the thesis live in Chapter 2 (§2.2.4 and
  Table A2.1); do not duplicate them here, to avoid the two drifting apart.
