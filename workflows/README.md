# workflows/

Reproducible analysis workflows demonstrating EpiProfile_PLANTS on real
datasets.

## Available workflows

| Directory | Dataset | Species | Acquisition | Status |
|-----------|---------|---------|-------------|--------|
| `PXD014739/` | PRIDE PXD014739 | *A. thaliana* | DDA | Reference workflow |

## Structure convention

Each workflow directory contains:

- `README` — Dataset description, parameters used, and how to reproduce.
- Output files (`.xls` ratio tables, `.pdf` XIC galleries).
- Reference publication PDF (when openly available).

## Adding a new workflow

1. Create a directory named after the dataset accession (e.g., `PXD099999/`).
2. Document the species, instrument, acquisition mode, and parameters.
3. Include the output ratio table and representative XIC gallery.
4. Add a README following the PXD014739 template.
