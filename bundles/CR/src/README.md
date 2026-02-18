# CR/src — Source code

All MATLAB `.m` files for the Chlamydomonas reinhardtii bundle, organised by
provenance tier.

## Subdirectories

| Directory | Tier | Count | Description |
|-----------|------|-------|-------------|
| `TIER1/`  | T1   | 72    | Upstream EpiProfile 2.0 basic — unchanged |
| `TIER2/`  | T2   | 10    | Modified for plant histones (forked from AT) |
| `TIER3/`  | T3   | 9     | Plant sequence-variant modules (AT-derived + CR-specific) |

## CR-specific changes vs. AT

- `init_histone0.m`: KTPATGGVKKPHR registered for H3 27-40 (instead of AT's
  KSAPTTGGVKKPHR and QSAPATGGVKKPHR).
- `H3_14_27_40.m` (TIER3): uses KTPATGGVKKPHR (13-residue variant).
- H4 20-23: uses canonical KVLR (same as human/AT; no L->F substitution).

## Entry point

```matlab
cd bundles/CR/src
EpiProfile
```
