function His = init_histone0(special)
%%
% EpiProfile-PLANTS init_histone0 for AT
% Generated: 2026-03-28, Updated: 2026-04-16
% 20 peptides, 75 peptidoforms
% Added: KGKYAER (H2A AT-specific N->K), LARYNKK (H2B AT-specific H->R, R->K)
% Sistema Pelamovic

no = 0;

% H3.1 18-26 (K18,K23)
no = no + 1;
His.out_filename{no,1} = 'H3_03_19_27';
His.pep_seq{no,1} = 'KQLATKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.1 27-40 (K27,K36)
no = no + 1;
His.out_filename{no,1} = 'H3_04_28_41';
His.pep_seq{no,1} = 'KSAPATGGVKKPHR';
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.1 3-8 (K4)
no = no + 1;
His.out_filename{no,1} = 'H3_01_4_9';
His.pep_seq{no,1} = 'TKQTAR';
His.mod_type{no,1} = '0,pr;2,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.1 9-17 (K9,K14)
no = no + 1;
His.out_filename{no,1} = 'H3_02_10_18';
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.3 27-40 (K27,K36)
no = no + 1;
His.out_filename{no,1} = 'H3_04_28_41';
His.pep_seq{no,1} = 'KSAPTTGGVKKPHR';
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 20-23 (K20)
no = no + 1;
His.out_filename{no,1} = 'H4_02_21_24';
His.pep_seq{no,1} = 'KVLR';
His.mod_type{no,1} = '0,pr;1,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 4-17 (K5,K8,K12,K16)
no = no + 1;
His.out_filename{no,1} = 'H4_01_5_18';
His.pep_seq{no,1} = 'GKGGKGLGKGGAKR';
His.mod_type{no,1} = '0,pr;2,pr;5,pr;9,pr;13,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 20-31 (K positions: [22])
no = no + 1;
His.out_filename{no,1} = 'HH2A_02_AT_20_31';
His.pep_seq{no,1} = 'SSKAGLQFPVGR';
His.mod_type{no,1} = '0,pr;3,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 35-44 (K positions: [37, 38, 40])
no = no + 1;
His.out_filename{no,1} = 'HH2A_03_AT_35_44';
His.pep_seq{no,1} = 'FLKKGKYAER';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 5-19 (K positions: [6, 14, 15])
no = no + 1;
His.out_filename{no,1} = 'HH2A_01_AT_5_19';
His.pep_seq{no,1} = 'GKTLGSGSAKKATTR';
His.mod_type{no,1} = '0,pr;2,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 3;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 74-79 (K positions: [76, 77])
no = no + 1;
His.out_filename{no,1} = 'HH2A_04_AT_74_79';
His.pep_seq{no,1} = 'DNKKTR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2B 109-118 (K positions)
no = no + 1;
His.out_filename{no,1} = 'HH2B_01_AT_109_118';
His.pep_seq{no,1} = 'YNKKPTITSR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2B 126-150 (K positions)
no = no + 1;
His.out_filename{no,1} = 'HH2B_02_AT_126_150';
His.pep_seq{no,1} = 'LVLPGELAKHAVSEGTKAVTKFTSS';
His.mod_type{no,1} = '0,pr;9,pr;17,pr;21,pr;';
His.pep_ch(no,1) = 4;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.1 73-83 (K79)
no = no + 1;
His.out_filename{no,1} = 'H3_07_74_84';
His.pep_seq{no,1} = 'EIAQDFKTDLR';
His.mod_type{no,1} = '0,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 24-35 (none (control))
no = no + 1;
His.out_filename{no,1} = 'H4_03_25_36';
His.pep_seq{no,1} = 'DNIQGITKPAIR';
His.mod_type{no,1} = '0,pr;8,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 79-92 (K91)
no = no + 1;
His.out_filename{no,1} = 'H4_05_80_93';
His.pep_seq{no,1} = 'KTVTAMDVVYALKR';
His.mod_type{no,1} = '0,pr;1,pr;13,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.1 117-128 (K122)
no = no + 1;
His.out_filename{no,1} = 'H3_08_118_129';
His.pep_seq{no,1} = 'VTIMPKDIQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 40-45 (none)
no = no + 1;
His.out_filename{no,1} = 'H4_04_42_46';
His.pep_seq{no,1} = 'GGVKR';
His.mod_type{no,1} = '0,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 36-42 AT-specific (K36,K38,Y39) — KGKYAER (human has KGNYAER, N->K)
% UniProt A0A178WDN2. K38 is a NEW modifiable site in AT.
no = no + 1;
His.out_filename{no,1} = 'HH2A_AT09_can_36_42';
His.pep_seq{no,1} = 'KGKYAER';
His.mod_type{no,1} = '0,pr;1,pr;3,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2B 80-86 AT-specific (Y) — LARYNKK (human has LAHYNKR, H->R + R->K)
% All 9 AT H2B variants. 100% conserved in AT.
no = no + 1;
His.out_filename{no,1} = 'HH2B_AT08_shared_80_86';
His.pep_seq{no,1} = 'LARYNKK';
His.mod_type{no,1} = '0,pr;6,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

%------------------------------------
function pep_mz = calculate_pepmz0(His,hno,special)
%%

Mods = GetMods();
if 4==special.nsource && (1==special.nsubtype || 3==special.nsubtype)
    aamass = GetaamassH();
else
    aamass = Getaamass();
end;
element = [12 1.0078246 14.0030732 15.9949141 31.972070];% element mass
mH2O = element(2)*2 + element(4);
pmass = 1.007276;

c_seq = His.pep_seq{hno};
idx = c_seq-'A'+1;
residuemass = aamass(idx,1)';
c_mod = His.mod_type{hno};
deltam = get_mod_mass(c_seq,c_mod,Mods);
% peptide+modification
residuemass_new = residuemass + deltam;
Mr = sum(residuemass_new)+mH2O;
c_ch = His.pep_ch(hno);
pep_mz = (Mr+c_ch*pmass)/c_ch;
