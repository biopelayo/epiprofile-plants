# assets/sequences/

Reference histone protein sequences used to derive the peptide panels for each
species bundle.

## Files

| File | Description |
|------|-------------|
| `H3_AT_MP_CR.fasta` | Histone H3 full-length sequences (human ref + AT + MP + CR) |
| `H4_AT_MP_CR.fasta` | Histone H4 full-length sequences (human ref + AT + MP + CR) |
| `canonical_AT_MP_CR.md` | Human-readable variant map for H3/H4: per-region substitution table |
| `canonical_H2A_H2B_H1_AT.md` | Peptide/variant map for H2A, H2B, H1 (AT), derived from the bundle modules + UniProt anchoring |
| `H2A_AT_peptides.fasta` | H2A peptide catalog (variants .1/.W/.X/.Z), one record per module |
| `H2B_AT_peptides.fasta` | H2B peptide catalog, one record per module |
| `H1_AT_peptides.fasta` | H1 peptide catalog (variants H1.1/H1.2), one record per module |
| `H2A_H2B_H1_AT_fulllength.fasta` | Full-length UniProt sequences anchoring each variant (verified by substring) |

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
- The `*_AT_MP_CR.fasta` files contain full-length protein sequences. The tryptic
  peptide sequences used by the MATLAB modules are derived from these after
  propionylation + trypsin digestion (see `docs/WHITEPAPER.md` Step 1).
- The `*_AT_peptides.fasta` files (H2A/H2B/H1) go the other way: they are peptide
  catalogs reverse-documented from the bundle modules, one record per `pep_seq`,
  because these families had no full-length reference. Full-length UniProt
  anchoring for them is a pending follow-up (see `canonical_H2A_H2B_H1_AT.md`).
