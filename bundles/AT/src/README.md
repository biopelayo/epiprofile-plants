# AT/src — Source code

All MATLAB `.m` files for the Arabidopsis thaliana bundle live here,
organised by provenance tier.

## Subdirectories

| Directory | Tier | Count | Description |
|-----------|------|-------|-------------|
| `TIER1/`  | T1   | 72    | Upstream EpiProfile 2.0 basic — unchanged |
| `TIER2/`  | T2   | 10    | Modified for plant histones (catalog, runner, QC) |
| `TIER3/`  | T3   | 9     | New modules for MSA-derived sequence variants |

## Entry point

`TIER2/EpiProfile.m` is the main entry point. It reads acquisition
parameters, converts RAW files (if needed), and dispatches the pipeline:

```
EpiProfile.m  ->  DrawISOProfile0  (build RT reference from top-3 runs)
              ->  DrawISOProfile1  (main quantification loop)
              ->  OutputTogether / OutputSinglePTMs / OutputFigures
```

## Adding the path

MATLAB needs all three tier directories on its path:

```matlab
addpath('TIER1','TIER2','TIER3');
```

`EpiProfile.m` does this automatically via `init_histone0.m`.
