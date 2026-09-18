# Canonical peptide map — H2A, H2B, H1 (Arabidopsis thaliana)

Reference peptide catalog for the linker (H1) and core (H2A, H2B) histone
families, **derived directly from the `pep_seq` fields of the MATLAB modules**
in `bundles/AT/src/`. Unlike `H3_AT_MP_CR.fasta` / `H4_AT_MP_CR.fasta` (which are
full-length MSA inputs), these families had **no reference sequence file at all**
even though the AT bundle ships extensive modules for them. This document closes
that gap by reverse-documenting what the code actually targets, one peptide per
module, so each peptide is traceable back to its `.m` source.

**Status**: AT only (the modules live in the AT bundle). MP and CR pending.

**Provenance note**: these are peptide-level catalogs, not reconstructed
full-length proteins — the modules cover discontinuous regions, so gaps between
listed regions are simply not targeted, not deletions. A full-length,
UniProt-anchored FASTA per variant is a recommended follow-up (see end).

Companion FASTA files (one record per module peptide):
`H2A_AT_peptides.fasta`, `H2B_AT_peptides.fasta`, `H1_AT_peptides.fasta`.

---

## H2A — four variant classes

The module suffix encodes the H2A variant: `oA1` = H2A type 1 (canonical),
`oAW` = H2A.W, `oAX` = H2A.X, `oAZ` = H2A.Z. The 45-84 core peptide is the
cleanest discriminator between them.

### H2A.1 (canonical, `oA1`)

| Region | Peptide | Module |
|--------|---------|--------|
| 1-7 | unmod | `HH2A_01u_1_7` (TIER2) |
| 7-14 | TLGSGVAK | `HH2A_01oA1_7_14` |
| 23-31 | AGLQFPVGR | `HH2A_02oA1_23_31` |
| 45-73 | VGAGAPVYLAAVLEYLAAEVLELAGNAAR | `HH2A_03oA1_45_73` |
| 84-90 | HIQLAVR | `HH2A_04oA1_84_90` |
| 91-97 | NDEELSK | `HH2A_05oA1_91_97` |
| 98-120 | LLGDVTIANGGVMPNIHSLLLPK | `HH2A_06oA1_98_120` |
| 122-132 | AGASKPSADED | `HH2A_07oA1_122_132` |

### H2A.W (`oAW`)

| Region | Peptide | Module |
|--------|---------|--------|
| 1-12 | MESSQATTKPTR | `HH2A_01oAW_1_12` |
| 13-17 | GAGGR | `HH2A_02oAW_13_17` |
| 32-40 | AGLQFPVGR | `HH2A_03oAW_32_40` |
| 54-82 | YGSGAPVYLAAVLEYLAAEVLELAGNAAR | `HH2A_04oAW_54_82` |
| 93-99 | HLCLAIR | `HH2A_05oAW_93_99` |
| 100-106 | NDEELGR | `HH2A_06oAW_100_106` |
| 107-129 | LLHGVTIASGGVLPNINPVLLPK | `HH2A_07oAW_107_129` |
| 131-140 | STASSSQAEK | `HH2A_08oAW_131_140` |
| 141-145 | ASATK | `HH2A_09oAW_141_145` |

### H2A.X (`oAX`)

| Region | Peptide | Module |
|--------|---------|--------|
| 1-11 | MSTGAGSGTTK | `HH2A_01oAX_1_11` |
| 29-37 | AGLQFPVGR | `HH2A_02oAX_29_37` |
| 51-79 | VGAGAPVYLSAVLEYLAAEVLELAGNAAR | `HH2A_03oAX_51_79` |
| 90-96 | HIQLAVR | `HH2A_04oAX_90_96` |
| 97-103 | NDEELSK | `HH2A_05oAX_97_103` |
| 104-127 | LLGSVTIANGGVLPNIHQTLLPSK | `HH2A_06oAX_104_127` |
| 133-142 | GDIGSASQEF | `HH2A_07oAX_133_142` |

### H2A.Z (`oAZ`)

| Region | Peptide | Module |
|--------|---------|--------|
| 8-19 | GLIMGKPSGSDK | `HH2A_01oAZ_8_19` |
| 25-29 | KPITR | `HH2A_02oAZ_25_29` |
| 33-41 | AGLQFPVGR | `HH2A_03oAZ_33_41` |
| 50-55 | STAHGR | `HH2A_04oAZ_50_55` |
| 56-85 | VGATAAAVYTAAILEYLTAEVLELAGNASK | `HH2A_05oAZ_56_84` |
| 102-111 | GDEELDTLIK | `HH2A_06oAZ_102_111` |
| 112-125 | GTIAGGGVIPHIHK | `HH2A_07oAZ_112_125` |
| 126-130 | SLINK | `HH2A_08oAZ_126_130` |

---

## H2B (single canonical form, `oA`)

Modules `HH2B_01..13` and `HH2B_14..26` are two passes over the same peptide set
(unique peptides listed once).

| Region | Peptide | Module |
|--------|---------|--------|
| 104-145 | unmod | `HH2B_01u_104_145` (TIER2) |
| 8-12 | KPAEK | `HH2B_01oA_8_12` / `HH2B_14oA_8_12` |
| 14-24 | TAAERPVEENK | `HH2B_02oA_14_24` / `HH2B_15oA_14_24` |
| 29-33 | APAEK | `HH2B_03oA_29_33` / `HH2B_16oA_29_33` |
| 45-49 | EAGDK | `HH2B_04oA_45_49` / `HH2B_17oA_45_49` |
| 57-62 | NVETYK | `HH2B_05oA_57_62` / `HH2B_18oA_57_62` |
| 63-67 | IYIFK | `HH2B_06oA_63_67` / `HH2B_19oA_63_67` |
| 71-81 | QVHPDIGISSK | `HH2B_07oA_71_81` / `HH2B_20oA_71_81` |
| 82-96 | AMGIMNSFINDIFEK | `HH2B_08oA_82_96` / `HH2B_21oA_82_96` |
| 97-103 | LAQESSK | `HH2B_09oA_97_103` / `HH2B_22oA_97_103` |
| 110-116 | KPTITSR | `HH2B_10oA_110_116` / `HH2B_23oA_110_116` |
| 117-123 | EIQTAVR | `HH2B_11oA_117_123` / `HH2B_24oA_117_123` |
| 124-132 | LVLPGELAK | `HH2B_12oA_124_132` / `HH2B_25oA_124_132` |
| 133-140 | HAVSEGTK | `HH2B_13oA_133_140` / `HH2B_26oA_133_140` |

---

## H1 — two variant classes

`oA1` and `oA2` are two AT H1 variants. They diverge substantially in both the
globular and C-terminal regions (e.g. 97-103 `ELPPTFR` vs `SLPPTFR`; 105-111
`LLLLNLK` vs `LLLVNLK`; 113-118 `LVASGK` vs `LVASEK`) and have different lengths
(oA1 reaches ~259, oA2 ~250), so they need separate modules throughout.

### H1 variant 1 (`oA1`)

| Region | Peptide | Module |
|--------|---------|--------|
| 1-29 | MSEVEIENAATIEGNTAADAPVTDAAVEK | `HH1_01oA1_1_29` |
| 30-34 | KPAAK | `HH1_02oA1_30_34` |
| 49-55 | TVAAAPK | `HH1_03oA1_49_55` |
| 58-71 | TVSHSHPTYEEMIK | `HH1_04oA1_58_70` |
| 71-77 | DAIVTLK | `HH1_05oA1_71_77` |
| 80-89 | TGSSQYAIQK | `HH1_06oA1_80_89` |
| 90-94 | FIEEK | `HH1_07oA1_90_94` |
| 97-103 | ELPPTFR | `HH1_08oA1_97_103` |
| 105-111 | LLLLNLK | `HH1_09oA1_105_111` |
| 113-118 | LVASGK | `HH1_10oA1_113_118` |
| 128-134 | LPSASAK | `HH1_11oA1_128_134` |
| 135-139 | ASSPK | `HH1_12oA1_135_139` |
| 140-144 | AAAEK | `HH1_13oA1_140_144` |
| 145-149 | SAPAK | `HH1_14oA1_145_149` |
| 151-159 | KPATVAVTK | `HH1_15oA1_151_159` |
| 164-169 | VAAASK | `HH1_16oA1_164_169` |
| 173-179 | TIAVKPK | `HH1_17oA1_173_179` |
| 180-184 | TAAAK | `HH1_18oA1_180_184` |
| 192-197 | AKPVPR | `HH1_19oA1_192_197` |
| 198-204 | ATAAATK | `HH1_20oA1_198_204` |
| 207-213 | AVDAKPK | `HH1_21oA1_207_213` |
| 216-220 | ARPAK | `HH1_22oA1_216_220` |
| 227-232 | VTSPAK | `HH1_23oA1_227_232` |
| 234-239 | AVAATK | `HH1_24oA1_234_239` |
| 241-247 | VATVATK | `HH1_25oA1_241_247` |
| 255-259 | VVKPK | `HH1_26oA1_255_259` |

### H1 variant 2 (`oA2`)

| Region | Peptide | Module |
|--------|---------|--------|
| 1-22 | MSIEEENVPTTVDSGAADTTVK | `HH1_01oA2_1_22` |
| 38-42 | TTTAK | `HH1_03oA2_38_42` |
| 50-55 | AAAPTK | `HH1_04oA2_50_55` |
| 58-70 | TTSSHPTYEEMIK | `HH1_05oA2_58_70` |
| 71-77 | DAIVTLK | `HH1_06oA2_71_77` |
| 80-89 | TGSSQYAIQK | `HH1_07oA2_80_89` |
| 90-94 | FIEEK | `HH1_08oA2_90_94` |
| 97-103 | SLPPTFR | `HH1_09oA2_97_103` |
| 105-111 | LLLVNLK | `HH1_10oA2_105_111` |
| 113-118 | LVASEK | `HH1_11oA2_113_118` |
| 128-132 | IPSAR | `HH1_12oA2_128_132` |
| 133-144 | SAATPKPAAPVK | `HH1_13oA2_133_144` |
| 147-154 | ATVVAKPK | `HH1_14oA2_147_154` |
| 157-165 | VAAAVAPAK | `HH1_15oA2_157_165` |
| 175-179 | KPAAK | `HH1_02oA2_175_179` |
| 186-191 | VTAKPK | `HH1_16oA2_186_191` |
| 194-200 | VTAAKPK | `HH1_17oA2_194_200` |
| 203-209 | SVAAVSK | `HH1_18oA2_203_209` |
| 212-218 | AVAAKPK | `HH1_19oA2_212_218` |
| 221-225 | ERPAK | `HH1_20oA2_221_225` |
| 233-237 | TSPGK | `HH1_21oA2_233_237` |
| 239-244 | VAAPAK | `HH1_22oA2_239_244` |
| 246-250 | VAVTK | `HH1_23oA2_246_250` |

---

## UniProt full-length anchoring

Each variant was matched to a reviewed UniProt AT entry by checking that the
module peptides are exact substrings. Full-length sequences are in
`H2A_H2B_H1_AT_fulllength.fasta`. Verification summary:

| Variant | Module code | UniProt | Gene | Peptides exact | Notes |
|---------|-------------|---------|------|----------------|-------|
| H2A.1 | oA1 | **Q9LHQ5** | HTA13 | 7/7 | clean |
| H2A.Z | oAZ | **Q9C944** | HTA9 | 8/8 | clean (56-85 corrected 2026-07-12: module `VGATAAVYT…` -> `VGATAAAVYT…`, +A) |
| H2A.W | oAW | **Q94F49** | HTA7 | 9/9 | clean |
| H2A.X | oAX | **O04848** | HTA5 (At1g08880) | 7/7 | clean; resolved by UniProt peptide search on `GDIGSASQEF`. N-term `MSTGAGSGTTK` matches HTA5, not HTA3 (Q9S9K7, `MSSGAGSGTTK`) |
| H2B | oA | **Q9LQQ4** | HTB1 | 13/13 | clean |
| H1.1 | oA1 | **P26568** | His1.1 | 26/26 | clean (58-71 corrected 2026-07-12: module `TVSSHPTYEEMIK` -> `TVSHSHPTYEEMIK`, +H) |
| H1.2 | oA2 | **P26569** | His1.2 | 23/23 | clean |

**Cross-check with `bundles/AT/metadata/build_AT_nucleosome_master.py`.** That
file already carried a `UNIPROT` dict, but with errors that this anchoring
corrected: H2A.W and H2A.X were swapped (Q94F49 is H2A.W, not H2A.X), and the
H3.1/H3.3 slots held H2A accessions as placeholders (now P59226 / P59169). The
dict has been updated accordingly. H2A_can (Q9LHQ5), H2A.Z (Q9C944), H2B.1
(Q9LQQ4), H1.1 (P26568) and H1.2 (P26569) were already correct and are confirmed.

## Remaining follow-up

1. **MP / CR.** No H2A/H2B/H1 modules exist for Marchantia or Chlamydomonas yet;
   extend once those bundles are built.

*(Resolved: H2A.X anchored to O04848/HTA5 via UniProt peptide search; the two
1-residue module discrepancies — H2A.Z 56-85 and H1.1 58-71 — were curation slips
and have been corrected in the modules to match Q9C944 / P26568. Every variant now
verifies at 100% of its module peptides.)*
