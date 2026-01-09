function Extract_C13_1(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath)
%%

%----------------------------
% H3K9ac|K14ac
His.out_filename = 'H3_0206_9_17';
His.pep_seq = 'KSTGGKAPR';
His.mod_short = {'K9ac|K14ac_12C2';
    'K9ac|K14ac_13C2'};
His.mod_type = {'0,pr;1,pr;6,ac;';
    '0,pr;1,pr;6,hac;'};
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H3K9me1K14ac
His.out_filename = 'H3_0207_9_17';
His.pep_seq = 'KSTGGKAPR';
His.mod_short = {'K9me1K14ac_12C2';
    'K9me1K14ac_13C2'};
His.mod_type = {'0,pr;1,me1;6,ac;';
    '0,pr;1,me1;6,hac;'};
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H3K9me2K14ac
His.out_filename = 'H3_0208_9_17';
His.pep_seq = 'KSTGGKAPR';
His.mod_short = {'K9me2K14ac_12C2';
    'K9me2K14ac_13C2'};
His.mod_type = {'0,pr;1,me2;6,ac;';
    '0,pr;1,me2;6,hac;'};
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H3K9me3K14ac
His.out_filename = 'H3_0209_9_17';
His.pep_seq = 'KSTGGKAPR';
His.mod_short = {'K9me3K14ac_12C2';
    'K9me3K14ac_13C2'};
His.mod_type = {'0,pr;1,me3;6,ac;';
    '0,pr;1,me3;6,hac;'};
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H3K9acK14ac
His.out_filename = 'H3_0210_9_17';
His.pep_seq = 'KSTGGKAPR';
His.mod_short = {'K9acK14ac_(12C2)2';
    'K9acK14ac_(13C2)1';
    'K9acK14ac_(13C2)2'};
His.mod_type = {'0,pr;1,ac;6,ac;';
    '0,pr;1,ac;6,hac;';
    '0,pr;1,hac;6,hac;';};
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

%----------------------------
% H3K18ac|K23ac
His.out_filename = 'H3_0306_18_26';
His.pep_seq = 'KQLATKAAR';
His.mod_short = {'K18ac|K23ac_12C2';
    'K18ac|K23ac_13C2'};
His.mod_type = {'0,pr;1,pr;6,ac;';
    '0,pr;1,pr;6,hac;'};
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H3K18acK23ac
His.out_filename = 'H3_0307_18_26';
His.pep_seq = 'KQLATKAAR';
His.mod_short = {'K18acK23ac_(12C2)2';
    'K18acK23ac_(13C2)1';
    'K18acK23ac_(13C2)2'};
His.mod_type = {'0,pr;1,ac;6,ac;';
    '0,pr;1,ac;6,hac;';
    '0,pr;1,hac;6,hac;'};
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

%----------------------------
% H3K27ac
His.out_filename = 'H3_0415_27_40';
His.pep_seq = 'KSAPATGGVKKPHR';
His.mod_short = {'K27ac_12C2';
    'K27ac_13C2'};
His.mod_type = {'0,pr;1,ac;10,pr;11,pr;';
    '0,pr;1,hac;10,pr;11,pr;'};
His.pep_ch = repmat([2 3 4],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 1;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

%----------------------------
% H33K27ac
His.out_filename = 'H3_04v315_27_40';
His.pep_seq = 'KSAPSTGGVKKPHR';
His.mod_short = {'K27ac_12C2';
    'K27ac_13C2'};
His.mod_type = {'0,pr;1,ac;10,pr;11,pr;';
    '0,pr;1,hac;10,pr;11,pr;'};
His.pep_ch = repmat([2 3 4],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 1;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

%----------------------------
% H4K5K8K12K16-1ac
His.out_filename = 'H4_0105_4_17';
His.pep_seq = 'GKGGKGLGKGGAKR';
His.mod_short = {'K5K8K12K16-1ac_12C2';
    'K5K8K12K16-1ac_13C2'};
His.mod_type = {'0,pr;2,pr;5,pr;9,pr;13,ac;';
    '0,pr;2,pr;5,pr;9,pr;13,hac;'};
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H4K5K8K12K16-2ac
His.out_filename = 'H4_0111_4_17';
His.pep_seq = 'GKGGKGLGKGGAKR';
His.mod_short = {'K5K8K12K16-2ac_(12C2)2';
    'K5K8K12K16-2ac_(13C2)1';
    'K5K8K12K16-2ac_(13C2)2'};
His.mod_type = {'0,pr;2,pr;5,pr;9,ac;13,ac;';
    '0,pr;2,pr;5,pr;9,ac;13,hac;';
    '0,pr;2,pr;5,pr;9,hac;13,hac;'};
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H4K5K8K12K16-3ac
His.out_filename = 'H4_0115_4_17';
His.pep_seq = 'GKGGKGLGKGGAKR';
His.mod_short = {'K5K8K12K16-3ac_(12C2)3';
    'K5K8K12K16-3ac_(13C2)1';
    'K5K8K12K16-3ac_(13C2)2';
    'K5K8K12K16-3ac_(13C2)3'};
His.mod_type = {'0,pr;2,pr;5,ac;9,ac;13,ac;';
    '0,pr;2,pr;5,ac;9,ac;13,hac;';
    '0,pr;2,pr;5,ac;9,hac;13,hac;';
    '0,pr;2,pr;5,hac;9,hac;13,hac;'};
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H4K5K8K12K16-4ac
His.out_filename = 'H4_0116_4_17';
His.pep_seq = 'GKGGKGLGKGGAKR';
His.mod_short = {'K5K8K12K16-4ac_(12C2)4';
    'K5K8K12K16-4ac_(13C2)1';
    'K5K8K12K16-4ac_(13C2)2';
    'K5K8K12K16-4ac_(13C2)3';
    'K5K8K12K16-4ac_(13C2)4'};
His.mod_type = {'0,pr;2,ac;5,ac;9,ac;13,ac;';
    '0,pr;2,ac;5,ac;9,ac;13,hac;';
    '0,pr;2,ac;5,ac;9,hac;13,hac;';
    '0,pr;2,ac;5,hac;9,hac;13,hac;';
    '0,pr;2,hac;5,hac;9,hac;13,hac;'};
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

%----------------------------
% H2AK5ac|K9ac
His.out_filename = 'HH2A_02m103_4_11';
His.pep_seq = 'GKQGGKAR';
His.mod_short = {'K5ac|K9ac_12C2';
    'K5ac|K9ac_13C2'};
His.mod_type = {'0,pr;2,pr;6,ac;';
    '0,pr;2,pr;6,hac;'};
His.pep_ch = repmat([1 2],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H2AK5acK9ac
His.out_filename = 'HH2A_02m104_4_11';
His.pep_seq = 'GKQGGKAR';
His.mod_short = {'K5acK9ac_(12C2)2';
    'K5acK9ac_(13C2)1';
    'K5acK9ac_(13C2)2'};
His.mod_type = {'0,pr;2,ac;6,ac;';
    '0,pr;2,ac;6,hac;';
    '0,pr;2,hac;6,hac;'};
His.pep_ch = repmat([1 2],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

%----------------------------
% H2AVK4K7K11K15-1ac
His.out_filename = 'HH2A_04oV05_1_19';
His.pep_seq = 'AGGKAGKDSGKAKAKAVSR';
His.mod_short = {'K4K7K11K15-1ac_12C2';
    'K4K7K11K15-1ac_13C2'};
His.mod_type = {'0,pr;4,pr;7,pr;11,pr;13,pr;15,ac;';
    '0,pr;4,pr;7,pr;11,pr;13,pr;15,hac;'};
His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H2AVK4K7K11K15-2ac
His.out_filename = 'HH2A_04oV11_1_19';
His.pep_seq = 'AGGKAGKDSGKAKAKAVSR';
His.mod_short = {'K4K7K11K15-2ac_(12C2)2';
    'K4K7K11K15-2ac_(13C2)1';
    'K4K7K11K15-2ac_(13C2)2'};
His.mod_type = {'0,pr;4,pr;7,pr;11,ac;13,pr;15,ac;';
    '0,pr;4,pr;7,pr;11,ac;13,pr;15,hac;';
    '0,pr;4,pr;7,pr;11,hac;13,pr;15,hac;'};
His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H2AVK4K7K11K15-3ac
His.out_filename = 'HH2A_04oV12_1_19';
His.pep_seq = 'AGGKAGKDSGKAKAKAVSR';
His.mod_short = {'K4K7K11K15-3ac_(12C2)3';
    'K4K7K11K15-3ac_(13C2)1';
    'K4K7K11K15-3ac_(13C2)2';
    'K4K7K11K15-3ac_(13C2)3'};
His.mod_type = {'0,pr;4,pr;7,ac;11,ac;13,pr;15,ac;';
    '0,pr;4,pr;7,ac;11,ac;13,pr;15,hac;';
    '0,pr;4,pr;7,ac;11,hac;13,pr;15,hac;';
    '0,pr;4,pr;7,hac;11,hac;13,pr;15,hac;'};
His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H2AVK4K7K11K15-4ac
His.out_filename = 'HH2A_04oV16_1_19';
His.pep_seq = 'AGGKAGKDSGKAKAKAVSR';
His.mod_short = {'K4K7K11K15-4ac_(12C2)4';
    'K4K7K11K15-4ac_(13C2)1';
    'K4K7K11K15-4ac_(13C2)2';
    'K4K7K11K15-4ac_(13C2)3';
    'K4K7K11K15-4ac_(13C2)4'};
His.mod_type = {'0,pr;4,ac;7,ac;11,ac;13,pr;15,ac;';
    '0,pr;4,ac;7,ac;11,ac;13,pr;15,hac;';
    '0,pr;4,ac;7,ac;11,hac;13,pr;15,hac;';
    '0,pr;4,ac;7,hac;11,hac;13,pr;15,hac;';
    '0,pr;4,hac;7,hac;11,hac;13,pr;15,hac;'};
His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

%----------------------------
% H2AZK4K7K11K15-1ac
His.out_filename = 'HH2A_04oZ05_1_19';
His.pep_seq = 'AGGKAGKDSGKAKTKAVSR';
His.mod_short = {'K4K7K11K15-1ac_12C2';
    'K4K7K11K15-1ac_13C2'};
His.mod_type = {'0,pr;4,pr;7,pr;11,pr;13,pr;15,ac;';
    '0,pr;4,pr;7,pr;11,pr;13,pr;15,hac;'};
His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H2AZK4K7K11K15-2ac
His.out_filename = 'HH2A_04oZ11_1_19';
His.pep_seq = 'AGGKAGKDSGKAKTKAVSR';
His.mod_short = {'K4K7K11K15-2ac_(12C2)2';
    'K4K7K11K15-2ac_(13C2)1';
    'K4K7K11K15-2ac_(13C2)2'};
His.mod_type = {'0,pr;4,pr;7,pr;11,ac;13,pr;15,ac;';
    '0,pr;4,pr;7,pr;11,ac;13,pr;15,hac;';
    '0,pr;4,pr;7,pr;11,hac;13,pr;15,hac;'};
His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H2AZK4K7K11K15-3ac
His.out_filename = 'HH2A_04oZ12_1_19';
His.pep_seq = 'AGGKAGKDSGKAKTKAVSR';
His.mod_short = {'K4K7K11K15-3ac_(12C2)3';
    'K4K7K11K15-3ac_(13C2)1';
    'K4K7K11K15-3ac_(13C2)2';
    'K4K7K11K15-3ac_(13C2)3'};
His.mod_type = {'0,pr;4,pr;7,ac;11,ac;13,pr;15,ac;';
    '0,pr;4,pr;7,ac;11,ac;13,pr;15,hac;';
    '0,pr;4,pr;7,ac;11,hac;13,pr;15,hac;';
    '0,pr;4,pr;7,hac;11,hac;13,pr;15,hac;'};
His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

% H2AZK4K7K11K15-4ac
His.out_filename = 'HH2A_04oZ16_1_19';
His.pep_seq = 'AGGKAGKDSGKAKTKAVSR';
His.mod_short = {'K4K7K11K15-4ac_(12C2)4';
    'K4K7K11K15-4ac_(13C2)1';
    'K4K7K11K15-4ac_(13C2)2';
    'K4K7K11K15-4ac_(13C2)3';
    'K4K7K11K15-4ac_(13C2)4'};
His.mod_type = {'0,pr;4,ac;7,ac;11,ac;13,pr;15,ac;';
    '0,pr;4,ac;7,ac;11,ac;13,pr;15,hac;';
    '0,pr;4,ac;7,ac;11,hac;13,pr;15,hac;';
    '0,pr;4,ac;7,hac;11,hac;13,pr;15,hac;';
    '0,pr;4,hac;7,hac;11,hac;13,pr;15,hac;'};
His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

%----------------------------
% H2AK13ac|K15ac
His.out_filename = 'HH2A_05m103_12_17';
His.pep_seq = 'AKAKTR';
His.mod_short = {'K13ac|K15ac_12C2';
    'K13ac|K15ac_13C2'};
His.mod_type = {'0,pr;2,pr;4,ac;';
    '0,pr;2,pr;4,hac;'};
His.pep_ch = repmat([1 2],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.display = ones(length(His.mod_type),1);
main_ch_idx = 2;
get_C13_info(MS1_index,MS1_peaks,MS2_index,MS2_peaks,special,ptol,cur_outpath,His,main_ch_idx);

fprintf(1,'\n');