# assets/sequences/

Reference histone protein sequences used to derive the peptide panels for each
species bundle.

## Files

| File | Description |
|------|-------------|
| `H3_AT_MP_CR.fasta` | Histone H3 full-length sequences (human ref + AT + MP + CR) |
| `H4_AT_MP_CR.fasta` | Histone H4 full-length sequences (human ref + AT + MP + CR) |
| `canonical_AT_MP_CR.md` | Human-readable variant map: per-region substitution table |

## Usage

These files are the input to Step 0 (MSA) of the adaptation methodology
described in `docs/WHITEPAPER.md`. The workflow is:

1. Align sequences in the FASTA files using ClustalW or MUSCLE.
2. Identify substituted peptide regions from the alignment.
3. Use the variant map (`canonical_AT_MP_CR.md`) to verify which TIER3
   modules are needed.
4. Update `init_histone0.m` in the species bundle with the new peptide catalog.

## Notes

- AT sequences are confirmed from MSA against UniProt references.
- MP and CR sequences are **preliminary** — marked in FASTA headers and
  pending full characterisation.
- The FASTA files contain full-length protein sequences. The tryptic peptide
  sequences used by the MATLAB modules are derived from these after
  propionylation + trypsin digestion (see `docs/WHITEPAPER.md` Step 1).
