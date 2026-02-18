# MP/src — Source code

All MATLAB `.m` files for the Marchantia polymorpha bundle, organised by
provenance tier.

## Subdirectories

| Directory | Tier | Count | Description |
|-----------|------|-------|-------------|
| `TIER1/`  | T1   | 72    | Upstream EpiProfile 2.0 basic — unchanged |
| `TIER2/`  | T2   | 10    | Modified for plant histones (forked from AT) |
| `TIER3/`  | T3   | 10    | Plant sequence-variant modules (AT-derived + MP-specific) |

## MP-specific changes vs. AT

- `init_histone0.m`: KSAPSTGGVKKPHR registered (instead of AT H3.3 KSAPTTGGVKKPHR
  and AT QSAPATGGVKKPHR); KVFR registered for H4 20-23.
- `H3_14_27_40.m` (TIER3): uses KSAPSTGGVKKPHR variant.
- `H4_02v_20_23.m` (TIER3): new module for KVFR (L->F at H4K20 region).
- `DrawISOProfile1.m`: calls H4_02v_20_23 in addition to upstream H4_02_20_23.

## Entry point

```matlab
cd bundles/MP/src
EpiProfile
```
