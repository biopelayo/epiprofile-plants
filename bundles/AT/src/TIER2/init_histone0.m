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
His.out_filename{no,1} = 'H3_06_53_63';
His.pep_seq{no,1} = 'KYQKSTELLIR';
His.mod_type{no,1} = '0,pr;1,pr;4,pr;';
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

no = no + 1;
His.out_filename{no,1} = 'H3_11_3_8';
His.pep_seq{no,1} = 'TKQTAR';
His.mod_type{no,1} = '0,pr;2,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_11_3_8';
His.pep_seq{no,1} = 'TKQSAR';
His.mod_type{no,1} = '0,pr;2,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_11_3_8';
His.pep_seq{no,1} = 'SNQTAR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_12_9_17';
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_12_9_17';
His.pep_seq{no,1} = 'KSTGGKGPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_12_9_17';
His.pep_seq{no,1} = 'KSHGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_12_9_17';
His.pep_seq{no,1} = 'ISTGGKAPR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_13_18_26';
His.pep_seq{no,1} = 'KQLATKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_13_18_26';
His.pep_seq{no,1} = 'KELATKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_13_18_26';
His.pep_seq{no,1} = 'TLLATKAAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_13_18_26';
His.pep_seq{no,1} = 'KQLAPKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_14_27_40';
His.pep_seq{no,1} = 'KSAPATGGVKKPHR';
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_14_27_40';
His.pep_seq{no,1} = 'KSAPTTGGVKKPHR';
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_14_27_40';
His.pep_seq{no,1} = 'QSAPATGGVKKPHR';
His.mod_type{no,1} = '0,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_16_53_63';
His.pep_seq{no,1} = 'KYQKSTELLIR';
His.mod_type{no,1} = '0,pr;1,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_16_53_63';
His.pep_seq{no,1} = 'KYQKSTELLNR';
His.mod_type{no,1} = '0,pr;1,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_17_73_83';
His.pep_seq{no,1} = 'EIAQDFKTDLR';
His.mod_type{no,1} = '0,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_17_73_83';
His.pep_seq{no,1} = 'EIAQDYKTDLR';
His.mod_type{no,1} = '0,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_18_117_128';
His.pep_seq{no,1} = 'VTIMPKDIQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_18_117_128';
His.pep_seq{no,1} = 'VTIMPKDVQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_18_117_128';
His.pep_seq{no,1} = 'VTIMPKEIQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
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

%------------------H2A------------------
no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';
His.pep_seq{no,1} = 'DNKKSR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';
His.pep_seq{no,1} = 'DNKKTR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';
His.pep_seq{no,1} = 'DNKKNR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';
His.pep_seq{no,1} = 'HLLLAIR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';
His.pep_seq{no,1} = 'HLQLAIR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';
His.pep_seq{no,1} = 'HLCLAIR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';
His.pep_seq{no,1} = 'HIQLAVR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';
His.pep_seq{no,1} = 'HVLLAVR';
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