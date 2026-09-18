function His = init_histone0(special)
%%
% EpiProfile-PLANTS init_histone0 for AT
% Generated: 2026-03-28, Repaired: 2026-04-22, Cleaned: 2026-05-24
% 18 peptides, ~66 peptidoforms
% Retired 2026-05-24 (Tribunal de Codigo Capa B + Art. 16 FN-1):
%   - KGKYAER (HH2A_AT09): sufijo 36-42 inconsistente con A0A178WDN2
%     (posicion real 38-44 1-based); ausente de
%     ArgC_digestion_all_histones.tsv -> peptido fantasma.
%   - LARYNKK (HH2B_AT08): termina en KK (no Arg-C); posicion 80-86
%     no existe en ninguna de las 9 isoformas H2B AT; ausente del TSV
%     teorico -> peptido fantasma.
%   Decision del usuario: retirar sin verificacion empirica previa.
% Sistema Pelamovic
%
% Cambios 2026-04-22 (repair RT=0 calibration bugs):
%   - Orden Yuan restaurado: H3_01, H3_02, H3_03, H3_04 en posiciones 1..4
%     (DrawISOProfile0 usa ipos=[1 2 3 4] para rankear rawfiles).
%   - out_filenames 0-based matching los modulos .m del bundle AT
%     (ej: H3_04_27_40 en lugar de H3_04_28_41).
%   - Duplicate 'H3_04_28_41' resuelto: H3.3 pasa a 'H3_04v3_27_40'.
%   - TKQTAR y KVLR: pep_ch=1 (era 2). En plantas z=1 es la carga dominante
%     para peptidos cortos; con z=2 el m/z no se detectaba y rt_ref salia 0.
%   - H4_04_40_45 secuencia corregida: 'GGVKR' (5aa) -> 'RGGVKR' (6aa) con
%     mod '0,pr;5,pr;'. El peptido real tras digestion Arg-C incluye la R(40)
%     terminal del peptido anterior. Esta seq coincide con His.pep_seq del
%     modulo H4_04_40_45.m (que tiene 'RGGVKR').

no = 0;

% H3.1 3-8 (K4). Carga z=2 es mas robusta entre datasets (z=1 fallo calibracion
% en PXD046034; los modulos .m aceptan [1 2]).
no = no + 1;
His.out_filename{no,1} = 'H3_01_3_8';
His.pep_seq{no,1} = 'TKQTAR';
His.mod_type{no,1} = '0,pr;2,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.1 9-17 (K9,K14)
no = no + 1;
His.out_filename{no,1} = 'H3_02_9_17';
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.1 18-26 (K18,K23)
no = no + 1;
His.out_filename{no,1} = 'H3_03_18_26';
His.pep_seq{no,1} = 'KQLATKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.1 27-40 (K27,K36,K37)
no = no + 1;
His.out_filename{no,1} = 'H3_04_27_40';
His.pep_seq{no,1} = 'KSAPATGGVKKPHR';
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 4-17 (K5,K8,K12,K16)
no = no + 1;
His.out_filename{no,1} = 'H4_01_4_17';
His.pep_seq{no,1} = 'GKGGKGLGKGGAKR';
His.mod_type{no,1} = '0,pr;2,pr;5,pr;9,pr;13,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 20-23 (K20)
no = no + 1;
His.out_filename{no,1} = 'H4_02_20_23';
His.pep_seq{no,1} = 'KVLR';
His.mod_type{no,1} = '0,pr;1,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3.3 27-40 (K27,K36,K37) — variante H3.3 (A31T)
no = no + 1;
His.out_filename{no,1} = 'H3_04v3_27_40';
His.pep_seq{no,1} = 'KSAPTTGGVKKPHR';
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3 73-83 (K79)
no = no + 1;
His.out_filename{no,1} = 'H3_07_73_83';
His.pep_seq{no,1} = 'EIAQDFKTDLR';
His.mod_type{no,1} = '0,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H3 117-128 (K122)
no = no + 1;
His.out_filename{no,1} = 'H3_08_117_128';
His.pep_seq{no,1} = 'VTIMPKDIQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 24-35 (none — control)
no = no + 1;
His.out_filename{no,1} = 'H4_03_24_35';
His.pep_seq{no,1} = 'DNIQGITKPAIR';
His.mod_type{no,1} = '0,pr;8,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 40-45 (K44). Seq real tras digestion Arg-C: RGGVKR (no GGVKR).
no = no + 1;
His.out_filename{no,1} = 'H4_04_40_45';
His.pep_seq{no,1} = 'RGGVKR';
His.mod_type{no,1} = '0,pr;5,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H4 79-92 (K91)
no = no + 1;
His.out_filename{no,1} = 'H4_06_79_92';
His.pep_seq{no,1} = 'KTVTAMDVVYALKR';
His.mod_type{no,1} = '0,pr;1,pr;13,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 20-31 (K22)
no = no + 1;
His.out_filename{no,1} = 'HH2A_AT01_can_20_31';
His.pep_seq{no,1} = 'SSKAGLQFPVGR';
His.mod_type{no,1} = '0,pr;3,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 5-19 (K6,K14,K15)
no = no + 1;
His.out_filename{no,1} = 'HH2A_AT02_can1_5_19';
His.pep_seq{no,1} = 'GKTLGSGSAKKATTR';
His.mod_type{no,1} = '0,pr;2,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 3;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 72-77 (K74,K75) — DNKKTR
no = no + 1;
His.out_filename{no,1} = 'HH2A_06m1_72_77';
His.pep_seq{no,1} = 'DNKKTR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2A 35-44 canonical variant (K37,K38,K40) — FLKKGKYAER. No hay modulo .m
% para esta seq; el init solo calcula m/z+rt para calibracion (sin uso directo).
no = no + 1;
His.out_filename{no,1} = 'HH2A_AT_can_35_44';
His.pep_seq{no,1} = 'FLKKGKYAER';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% [Retirado 2026-05-24] HH2A_AT09_can_36_42 (KGKYAER) - peptido fantasma.
% Ausente de ArgC_digestion_all_histones.tsv. Sufijo de posicion
% inconsistente con A0A178WDN2 (era 36-42, real 38-44 1-based).

% H2B 108-117 (K110,K111) — YNKKPTITSR
no = no + 1;
His.out_filename{no,1} = 'HH2B_AT02_shared_108_117';
His.pep_seq{no,1} = 'YNKKPTITSR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% H2B 126-150 (K134,K142,K146) — peptido largo 25 aa. Sin modulo .m exacto
% (modulo HH2B_AT03 tiene variante truncada de 17 aa, residuos 118-134).
no = no + 1;
His.out_filename{no,1} = 'HH2B_AT_shared_126_150';
His.pep_seq{no,1} = 'LVLPGELAKHAVSEGTKAVTKFTSS';
His.mod_type{no,1} = '0,pr;9,pr;17,pr;21,pr;';
His.pep_ch(no,1) = 4;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

% [Retirado 2026-05-24] HH2B_AT08_shared_80_86 (LARYNKK) - peptido fantasma.
% Termina en KK (no Arg-C). Ausente de ArgC_digestion_all_histones.tsv.
% Posicion 80-86 no existe en ninguna isoforma H2B AT.

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
