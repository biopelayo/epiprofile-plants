# AT — *Arabidopsis thaliana* bundle

First and reference species bundle for EpiProfile_PLANTS.

## Status

**Active** — validated against PRIDE dataset PXD014739 (Arabidopsis histone
extracts, DDA acquisition).

## Source data

| Item | Value |
|------|-------|
| Organism | *Arabidopsis thaliana* (Col-0) |
| Histone variants | H3.1, H3.3, H4 |
| Derivatization | Propionylation + trypsin (Arg-C-like specificity) |
| Instrument modes | DDA and DIA supported |

## Key sequence differences vs. human

H3 27-40 carries two plant-specific isoforms (H3.3: KSAPTTGGVKKPHR;
AtH3.1-like: QSAPATGGVKKPHR) that required new TIER2/TIER3 modules.
Multiple other substitutions across H3 3-8, 9-17, 18-26, 53-63, 73-83,
117-128 are handled by dedicated TIER3 modules.

## Running

```matlab
cd bundles/AT/src
EpiProfile      % follow GUI prompts or set params in check_otherparas.m
```

See `docs/MANUAL.md` for detailed instructions.
