# Canonical sequence notes — AT, MP, CR

This document summarises the key histone sequence variants across the three
target species, with positions mapped to the propionylation + trypsin peptide
regions used by EpiProfile.

**Status**: AT variants confirmed from MSA; MP and CR variants are preliminary
(pending full characterisation).

---

## H3 peptide regions — variant map

### Region 3-8 (H3_01 / H3_11)

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | TKQTAR | Reference |
| AT (canonical) | TKQTAR | Same as human |
| AT (variant 1) | TKQSAR | T->S at pos 5 |
| AT (variant 2) | SNQTAR | TK->SN at pos 3-4 |
| MP | TKQTAR | Same as human (preliminary) |
| CR | TKQTAR | Same as human (preliminary) |

### Region 9-17 (H3_02 / H3_12)

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | KSTGGKAPR | Reference |
| AT (canonical) | KSTGGKAPR | Same as human |
| AT (variant 1) | KSTGGKGPR | A->G at pos 15 |
| AT (variant 2) | KSHGGKAPR | T->H at pos 12 (likely H3.1 isoform) |
| AT (variant 3) | ISTGGKAPR | K->I at pos 9 |
| MP | pending | |
| CR | pending | |

### Region 18-26 (H3_03 / H3_13)

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | KQLATKAAR | Reference |
| AT (canonical) | KQLATKAAR | Same as human |
| AT (variant 1) | KELATKAAR | Q->E at pos 19 |
| AT (variant 2) | TLLATKAAR | K->T at pos 18, Q->L at pos 19 |
| AT (variant 3) | KQLAPKAAR | T->P at pos 22 |
| MP | pending | |
| CR | pending | |

### Region 27-40 (H3_04 / H3_14) — most divergent region

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | KSAPATGGVKKPHR | Reference |
| AT H3.1 | KSAPATGGVKKPHR | Same as human |
| AT H3.3 | KSAPTTGGVKKPHR | A->T at pos 31 |
| AT (variant) | QSAPATGGVKKPHR | K->Q at pos 27 |
| **MP** | **KSAPSTGGVKKPHR** | A->S at pos 31 |
| **CR** | **KTPATGGVKKPHR** | SA->T at pos 28-29 (shorter motif) |

This is the region with the most inter-species variation and requires a
dedicated TIER3 module for each species.

### Region 41-49 (H3_05) — H3.1 vs H3.3 discriminator

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | YRPGTVALR | Reference (Y41) |
| AT H3.1 | FRPGTVALR | **Y41F** — canonical plant H3.1 residue (UniProt P59226, HTR1/2/3/9/13) |
| AT H3.3 | YRPGTVALR | Y41 — same as human (UniProt P59169, HTR4/5/8) |
| MP | pending | |
| CR | pending | |

Position 41 is the canonical H3.1/H3.3 discriminator in plants (F in H3.1, Y in
H3.3). The existing module `H3_05_41_49.m` uses `YRPGTVALR`, i.e. it tracks the
**H3.3** form only. The H3.1 sibling `FRPGTVALR` has no module yet — see note
below. This peptide carries no PTM site, so it is used mainly as a reference /
normaliser; the missing H3.1 form means the H3.1 pool is not captured at 41-49.

### Region 53-63 (H3_06 / H3_16)

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | RYQKSTELLIR | Reference (R53; peptide 54-63 after tryptic cut) |
| AT (canonical) | KYQKSTELLIR | **R53K** — K53 is propionylated, blocking the cut, so the peptide extends to 53-63. Both AT H3.1 (P59226) and H3.3 (P59169) carry K53. Handled by `H3_06` (K variant) + `H3_06a` (R variant) |
| AT (variant) | KYQKSTELLNR | I->N at pos 62 |
| MP | pending | |
| CR | pending | |

### Region 73-83 (H3_07 / H3_17)

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | EIAQDFKTDLR | Reference |
| AT (canonical) | EIAQDFKTDLR | Same as human |
| AT (variant) | EIAQDYKTDLR | F->Y at pos 78 |
| MP | pending | |
| CR | pending | |

### Region 117-128 (H3_08 / H3_18)

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | VTIMPKDIQLAR | Reference |
| AT (canonical) | VTIMPKDIQLAR | Same as human |
| AT (variant 1) | VTIMPKDVQLAR | I->V at pos 123 |
| AT (variant 2) | VTIMPKEIQLAR | D->E at pos 123 |
| MP | pending | |
| CR | pending | |

> **H3 globular domain (AT) — resolved.** Only the N-terminal tail regions had
> been hand-curated; the globular-domain residues were left as the human
> sequence. Both records are now corrected against UniProt:
> - `H3.1_AT` = P59226 (HTR1/2/3/9/13): **F41, K53, A90, A96**.
> - `H3.3_AT` = P59169 (HTR4/5/8): **T31, Y41, K53, H87, L90, A96**.
>
> The four residues that separate plant H3.1 from H3.3 are positions **31 (A/T),
> 41 (F/Y), 87 (S/H), 90 (A/L)**; both carry K53 and A96. (An earlier note about a
> `VTIIMPK` stretch in H3.3 was a transcription error — P59169 has the normal
> `VTIMPK`.)
>
> **One open item — H3.1 41-49 module.** `H3_05_41_49.m` encodes only the H3.3
> form `YRPGTVALR`. A draft H3.1 sibling `H3_05b_41_49.m` (`FRPGTVALR`) has been
> added but is **not registered** in the runner: `FRPGTVALR` has no PTM site
> (unmod only, used as a normaliser), so it needs an RT calibration and a PSM
> check that the F41 peptide appears in the data before being enabled.

---

## H4 peptide regions — variant map

### Region 4-17 (H4_01)

All three species share the human sequence: **GKGGKGLGKGGAKR**

### Region 20-23 (H4_02)

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | KVLR | Reference |
| AT | KVLR | Same as human |
| **MP** | **KVFR** | **L->F at pos 21** |
| CR | KVLR | Same as human (preliminary) |

### Regions 24-35, 40-45, 79-92

No known variants across the three species. Human sequences apply.

### Region 56-67 (H4 core) — V60I in AT

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | GVLKVFLENVIR | Reference (V60) |
| **AT** | **GVLKIFLENVIR** | **V->I at pos 60** (UniProt P59259). MS²-level correction; +45 PSMs empirically |
| MP | pending | |
| CR | pending | |

### Region 68-78 (H4_05) — K77R in AT

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | DAVTYTEHAKR | Reference (K77, tryptic 68-78) |
| **AT** | **DAVTYTEHAR** | **K77->R** (UniProt P59259). R77 introduces an Arg-C/tryptic cut site, so the peptide shortens to 68-77 (drops the terminal R78). Empirically 0 -> 207 PSMs. Implemented in `H4_05_68_78.m` (`His.pep_seq = 'DAVTYTEHAR'`) |
| MP | pending | |
| CR | pending | |

> **Note (H4 AT):** UniProt P62805 (human) vs P59259 (AT) differ at exactly two
> positions in the mature protein: V60I and K77R. The K77R substitution is
> functionally relevant because it changes the trypsin/propionylation digestion
> pattern (shorter 68-77 peptide), not just an MS² mass.

---

## Summary: modules needed per species

| Region | AT | MP | CR |
|--------|-----|-----|-----|
| H3 3-8 | T3 (3 sequences) | T1 (pending MSA) | T1 (pending MSA) |
| H3 9-17 | T3 (4 sequences) | pending | pending |
| H3 18-26 | T3 (4 sequences) | pending | pending |
| H3 27-40 | T2+T3 (3 sequences) | **new T3 needed** | **new T3 needed** |
| H3 41-49 | T1 (H3.3 `YRPGTVALR` + H3.1 `FRPGTVALR` **draft, unregistered**) | pending | pending |
| H3 53-63 | T3 (2 sequences) | pending | pending |
| H3 73-83 | T3 (2 sequences) | pending | pending |
| H3 117-128 | T3 (3 sequences) | pending | pending |
| H4 20-23 | T1 | **new T3 needed** | T1 (pending) |
| H4 56-67 | **T3 (V60I)** | pending | pending |
| H4 68-78 | **T3 (K77R, 68-77)** | pending | pending |
| All others | T1 | T1 | T1 |

---

## UniProt references

| Histone | Species | UniProt ID | Notes |
|---------|---------|------------|-------|
| H3.1 | Human | P68431 | Alignment reference |
| H4 | Human | P62805 | Alignment reference |
| H3.1 | AT | multiple | See Phytozome / TAIR |
| H3.3 | AT | multiple | HTR4, HTR5 loci |
| H3.1 | MP | pending | |
| H3 | CR | pending | |

---

*This file is a human-readable reference note. Do not treat it as a
machine-readable table. For programmatic use, see `init_histone0.m` in the
species bundle.*
