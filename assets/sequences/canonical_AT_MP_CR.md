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

### Region 53-63 (H3_06 / H3_16)

| Species | Sequence | Notes |
|---------|----------|-------|
| Human | KYQKSTELLIR | Reference |
| AT (canonical) | KYQKSTELLIR | Same as human |
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

### Regions 24-35, 40-45, 68-78, 79-92

No known variants across the three species. Human sequences apply.

---

## Summary: modules needed per species

| Region | AT | MP | CR |
|--------|-----|-----|-----|
| H3 3-8 | T3 (3 sequences) | T1 (pending MSA) | T1 (pending MSA) |
| H3 9-17 | T3 (4 sequences) | pending | pending |
| H3 18-26 | T3 (4 sequences) | pending | pending |
| H3 27-40 | T2+T3 (3 sequences) | **new T3 needed** | **new T3 needed** |
| H3 53-63 | T3 (2 sequences) | pending | pending |
| H3 73-83 | T3 (2 sequences) | pending | pending |
| H3 117-128 | T3 (3 sequences) | pending | pending |
| H4 20-23 | T1 | **new T3 needed** | T1 (pending) |
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
