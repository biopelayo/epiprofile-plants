function His = init_histone0(special)
%%

no = 0;

% checklist
no = no + 1;
His.out_filename{no,1} = 'H3_01_3_8';
His.pep_seq{no,1} = 'TKQTAR';
His.mod_type{no,1} = '0,pr;2,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_02_9_17';
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_03_18_26';
His.pep_seq{no,1} = 'KQLATKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_04_27_40';
His.pep_seq{no,1} = 'KSAPATGGVKKPHR';
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_01_4_17';
His.pep_seq{no,1} = 'GKGGKGLGKGGAKR';
His.mod_type{no,1} = '0,pr;2,pr;5,pr;9,pr;13,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_02_20_23';
His.pep_seq{no,1} = 'KVLR';
His.mod_type{no,1} = '0,pr;1,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

%------------------H3------------------
no = no + 1;
His.out_filename{no,1} = 'H3_02a_9_17';
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;2,ph;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_02b_9_17';
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;2,ac;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_05_41_49';
His.pep_seq{no,1} = 'YRPGTVALR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_06_54_63';
His.pep_seq{no,1} = 'YQKSTELLIR';
His.mod_type{no,1} = '0,pr;3,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_06a_53_63';
His.pep_seq{no,1} = 'RYQKSTELLIR';
His.mod_type{no,1} = '0,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_07_73_83';
His.pep_seq{no,1} = 'EIAQDFKTDLR';
His.mod_type{no,1} = '0,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_08_117_128';
His.pep_seq{no,1} = 'VTIMPKDIQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_09u_64_135';
His.pep_seq{no,1} = 'KLPFQR';
His.mod_type{no,1} = '0,pr;1,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_09u_64_135';
His.pep_seq{no,1} = 'RIRGER';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_09u_64_135';
His.pep_seq{no,1} = 'IRGERA';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

%------------------H4------------------
no = no + 1;
His.out_filename{no,1} = 'H4_02a_18_23';
His.pep_seq{no,1} = 'HRKVLR';
His.mod_type{no,1} = '0,pr;3,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_02b_20_35';
His.pep_seq{no,1} = 'KVLRDNIQGITKPAIR';
His.mod_type{no,1} = '0,pr;1,pr;12,pr;';
His.pep_ch(no,1) = 3;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_02c_20_36';
His.pep_seq{no,1} = 'KVLRDNIQGITKPAIRR';
His.mod_type{no,1} = '0,pr;1,pr;12,pr;';
His.pep_ch(no,1) = 3;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_03_24_35';
His.pep_seq{no,1} = 'DNIQGITKPAIR';
His.mod_type{no,1} = '0,pr;8,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_04_40_45';
His.pep_seq{no,1} = 'RGGVKR';
His.mod_type{no,1} = '0,pr;5,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_05_68_78';
His.pep_seq{no,1} = 'DAVTYTEHAKR';
His.mod_type{no,1} = '0,pr;10,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_06_79_92';
His.pep_seq{no,1} = 'KTVTAMDVVYALKR';
His.mod_type{no,1} = '0,pr;1,pr;13,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_07u_24_102';
His.pep_seq{no,1} = 'DNIQGITKPAIRR';
His.mod_type{no,1} = '0,pr;8,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_07u_24_102';
His.pep_seq{no,1} = 'ISGLIYEETR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_07u_24_102';
His.pep_seq{no,1} = 'GVLKVFLENVIR';
His.mod_type{no,1} = '0,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_07u_24_102';
His.pep_seq{no,1} = 'TLYGFGG';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

%------------------H1------------------
if 1==special.norganism
    no = no + 1;
    His.out_filename{no,1} = 'HH1_01o4_25_32';
    His.pep_seq{no,1} = 'KSAGAAKR';
    His.mod_type{no,1} = '0,pr;1,pr;7,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_02m2_33_53';
    His.pep_seq{no,1} = 'KASGPPVSELITKAVAASKER';
    His.mod_type{no,1} = '0,pr;1,pr;13,pr;19,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_03o5_36_56';
    His.pep_seq{no,1} = 'KATGPPVSELITKAVAASKER';
    His.mod_type{no,1} = '0,pr;1,pr;13,pr;19,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_04v_1_35';
    His.pep_seq{no,1} = 'SETAPAAPAAAPPAEKAPVKKKAAKKAGGTPR';
    His.mod_type{no,1} = '0,pr;16,pr;20,pr;21,pr;22,pr;25,pr;26,pr;';
    His.pep_ch(no,1) = 4;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_04v_1_35';
    His.pep_seq{no,1} = 'SETAPLAPTIPAPAEKTPVKKKAKKAGATAGKR';
    His.mod_type{no,1} = '0,pr;16,pr;20,pr;21,pr;22,pr;24,pr;25,pr;32,pr;';
    His.pep_ch(no,1) = 4;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_04v_1_35';
    His.pep_seq{no,1} = 'SETAPAAPAAPAPAEKTPVKKKAR';
    His.mod_type{no,1} = '0,pr;16,pr;20,pr;21,pr;22,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_04v_1_35';
    His.pep_seq{no,1} = 'SETAPAETATPAPVEKSPAKKKATKKAAGAGAAKR';
    His.mod_type{no,1} = '0,pr;16,pr;20,pr;21,pr;22,pr;25,pr;26,pr;34,pr;';
    His.pep_ch(no,1) = 4;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_05v_54_81';
    His.pep_seq{no,1} = 'GGVSLAALKKALAAAGYDVEKNNSR';
    His.mod_type{no,1} = '0,pr;9,pr;10,pr;21,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_05v_54_81';
    His.pep_seq{no,1} = 'SGVSLAALKKALAAAGYDVEKNNSR';
    His.mod_type{no,1} = '0,pr;9,pr;10,pr;21,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_05v_54_81';
    His.pep_seq{no,1} = 'NGLSLAALKKALAAGGYDVEKNNSR';
    His.mod_type{no,1} = '0,pr;9,pr;10,pr;21,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_05v_54_81';
    His.pep_seq{no,1} = 'VGMSLVALKKALAAAGYDVEKNNSR';
    His.mod_type{no,1} = '0,pr;9,pr;10,pr;21,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_06u_1_207';
    His.pep_seq{no,1} = 'TENSTSAPAAKPKR';
    His.mod_type{no,1} = '0,pr;11,pr;13,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_06u_1_207';
    His.pep_seq{no,1} = 'KQGGAAKDTR';
    His.mod_type{no,1} = '0,pr;1,pr;7,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_06u_1_207';
    His.pep_seq{no,1} = 'GRKPAGLISASR';
    His.mod_type{no,1} = '0,pr;3,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_06u_1_207';
    His.pep_seq{no,1} = 'KKLEGGGER';
    His.mod_type{no,1} = '0,pr;1,pr;2,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
else
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_01o4_25_32';
    His.pep_seq{no,1} = 'KAAGGAKR';
    His.mod_type{no,1} = '0,pr;1,pr;7,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1_02m2_33_53';
    His.pep_seq{no,1} = 'KASGPPVSELITKAVAASKER';
    His.mod_type{no,1} = '0,pr;1,pr;13,pr;19,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_03o5_33_53';
    His.pep_seq{no,1} = 'KATGPPVSELITKAVSASKER';
    His.mod_type{no,1} = '0,pr;1,pr;13,pr;19,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_04v_1_35';
    His.pep_seq{no,1} = 'SETAPVAQAASTATEKPAAAKKTKKPAKAAAPR';
    His.mod_type{no,1} = '0,pr;16,pr;21,pr;22,pr;24,pr;25,pr;28,pr;';
    His.pep_ch(no,1) = 4;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_04v_1_35';
    His.pep_seq{no,1} = 'SEAAPAAPAAAPPAEKAPAKKKAAKKPAGVR';
    His.mod_type{no,1} = '0,pr;16,pr;20,pr;21,pr;22,pr;25,pr;26,pr;';
    His.pep_ch(no,1) = 4;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_04v_1_35';
    His.pep_seq{no,1} = 'SETAPAAPAAPAPVEKTPVKKKAKKTGAAAGKR';
    His.mod_type{no,1} = '0,pr;16,pr;20,pr;21,pr;22,pr;24,pr;25,pr;32,pr;';
    His.pep_ch(no,1) = 4;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_04v_1_35';
    His.pep_seq{no,1} = 'SETAPAAPAAPAPAEKTPVKKKAR';
    His.mod_type{no,1} = '0,pr;16,pr;20,pr;21,pr;22,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_04v_1_35';
    His.pep_seq{no,1} = 'SETAPAETAAPAPVEKSPAKKKTTKKAGAAKR';
    His.mod_type{no,1} = '0,pr;16,pr;20,pr;21,pr;22,pr;25,pr;26,pr;31,pr;';
    His.pep_ch(no,1) = 4;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_05v_54_81';
    His.pep_seq{no,1} = 'SGVSLAALKKSLAAAGYDVEKNNSR';
    His.mod_type{no,1} = '0,pr;9,pr;10,pr;21,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_05v_54_81';
    His.pep_seq{no,1} = 'SGVSLAALKKALAAAGYDVEKNNSR';
    His.mod_type{no,1} = '0,pr;9,pr;10,pr;21,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_05v_54_81';
    His.pep_seq{no,1} = 'GGVSLPALKKALAAGGYDVEKNNSR';
    His.mod_type{no,1} = '0,pr;9,pr;10,pr;21,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH1Mo_05v_54_81';
    His.pep_seq{no,1} = 'AGMSLAALKKALAAAGYDVEKNNSR';
    His.mod_type{no,1} = '0,pr;9,pr;10,pr;21,pr;';
    His.pep_ch(no,1) = 3;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
end;

%------------------H2A------------------
no = no + 1;
His.out_filename{no,1} = 'HH2A_01m1_36_42';
His.pep_seq{no,1} = 'KGNYAER';
His.mod_type{no,1} = '0,pr;1,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01m3_36_42';
His.pep_seq{no,1} = 'KGNYSER';
His.mod_type{no,1} = '0,pr;1,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01oX_36_42';
His.pep_seq{no,1} = 'KGHYAER';
His.mod_type{no,1} = '0,pr;1,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_02m1_4_11';
His.pep_seq{no,1} = 'GKQGGKAR';
His.mod_type{no,1} = '0,pr;2,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_02oJ_4_11';
His.pep_seq{no,1} = 'GKQGGKVR';
His.mod_type{no,1} = '0,pr;2,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_02oX_4_11';
His.pep_seq{no,1} = 'GKTGGKAR';
His.mod_type{no,1} = '0,pr;2,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_03m1_1_11';
His.pep_seq{no,1} = 'SGRGKQGGKAR';
His.mod_type{no,1} = '0,pr;5,pr;9,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_04oV_1_19';
His.pep_seq{no,1} = 'AGGKAGKDSGKAKAKAVSR';
His.mod_type{no,1} = '0,pr;4,pr;7,pr;11,pr;13,pr;15,pr;';
His.pep_ch(no,1) = 3;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_04oZ_1_19';
His.pep_seq{no,1} = 'AGGKAGKDSGKAKTKAVSR';
His.mod_type{no,1} = '0,pr;4,pr;7,pr;11,pr;13,pr;15,pr;';
His.pep_ch(no,1) = 3;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_05m1_12_17';
His.pep_seq{no,1} = 'AKAKTR';
His.mod_type{no,1} = '0,pr;2,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_05m3_12_17';
His.pep_seq{no,1} = 'AKAKSR';
His.mod_type{no,1} = '0,pr;2,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_06m1_72_77';
His.pep_seq{no,1} = 'DNKKTR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_07v_1_88';
His.pep_seq{no,1} = 'HLQLAIR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_07v_1_88';
His.pep_seq{no,1} = 'AGGKAGKDSGKAKTKAVSR';
His.mod_type{no,1} = '0,pr;4,pr;7,pr;11,pr;13,pr;15,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

if 1==special.norganism
    no = no + 1;
    His.out_filename{no,1} = 'HH2A_07v_1_88';
    His.pep_seq{no,1} = 'SAKAGVIFPVGR';
    His.mod_type{no,1} = '0,pr;3,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
end;

no = no + 1;
His.out_filename{no,1} = 'HH2A_07v_1_88';
His.pep_seq{no,1} = 'GKTGGKAR';
His.mod_type{no,1} = '0,pr;2,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_08u_4_99';
His.pep_seq{no,1} = 'HLQLAVR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

if 1==special.norganism
    no = no + 1;
    His.out_filename{no,1} = 'HH2A_08u_4_99';
    His.pep_seq{no,1} = 'GGKKKSTKTSR';
    His.mod_type{no,1} = '0,pr;3,pr;4,pr;5,pr;8,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH2A_08u_4_99';
    His.pep_seq{no,1} = 'SGKKKMSKLSR';
    His.mod_type{no,1} = '0,pr;3,pr;4,pr;5,pr;8,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
end;

no = no + 1;
His.out_filename{no,1} = 'HH2A_08u_4_99';
His.pep_seq{no,1} = 'IHRHLKTR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_08u_4_99';
His.pep_seq{no,1} = 'IHRHLKSR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_08u_4_99';
His.pep_seq{no,1} = 'NDEELNKLLGR';
His.mod_type{no,1} = '0,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_08u_4_99';
His.pep_seq{no,1} = 'AGLQFPVGR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_08u_4_99';
His.pep_seq{no,1} = 'VHRLLR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

if 1==special.norganism
    no = no + 1;
    His.out_filename{no,1} = 'HH2A_08u_4_99';
    His.pep_seq{no,1} = 'IHPELLAKKR';
    His.mod_type{no,1} = '0,pr;8,pr;9,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH2A_08u_4_99';
    His.pep_seq{no,1} = 'YIKKGHPKYR';
    His.mod_type{no,1} = '0,pr;3,pr;4,pr;8,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
end;

%------------------H2B------------------
no = no + 1;
His.out_filename{no,1} = 'HH2B_01m1B_80_86';
His.pep_seq{no,1} = 'LAHYNKR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

if 1==special.norganism
    no = no + 1;
    His.out_filename{no,1} = 'HH2B_01o1A_81_87';
    His.pep_seq{no,1} = 'LAHYSKR';
    His.mod_type{no,1} = '0,pr;6,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
end;

no = no + 1;
His.out_filename{no,1} = 'HH2B_02v_1_29';
His.pep_seq{no,1} = 'PEPAKSAPAPKKGSKKAVTKAQKKDGKKR';
His.mod_type{no,1} = '0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;';
His.pep_ch(no,1) = 4;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

if 1==special.norganism
    no = no + 1;
    His.out_filename{no,1} = 'HH2B_03u_1_100';
    His.pep_seq{no,1} = 'LPHYNKR';
    His.mod_type{no,1} = '0,pr;6,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    
    no = no + 1;
    His.out_filename{no,1} = 'HH2B_03u_1_100';
    His.pep_seq{no,1} = 'LRTEVPRLPR';
    His.mod_type{no,1} = '0,pr;';
    His.pep_ch(no,1) = 2;
    His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
    new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
    His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
end;

no = no + 1;
His.out_filename{no,1} = 'HH2B_03u_1_100';
His.pep_seq{no,1} = 'EIQTAVR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2B_03u_1_100';
His.pep_seq{no,1} = 'EVQTAVR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2B_03u_1_100';
His.pep_seq{no,1} = 'IAGEASR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2B_03u_1_100';
His.pep_seq{no,1} = 'IASEASR';
His.mod_type{no,1} = '0,pr;';
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