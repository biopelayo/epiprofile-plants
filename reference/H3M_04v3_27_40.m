function H3M_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_04v3_27_40';
% fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% init
His = init_histone();

% relocate
His = relocate(cur_outpath,out_filename,His);

% calculate
unitdiff = 1.0032;
[pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,ptol,unitdiff,His,special);

% output
output_histone2(cur_outpath,out_filename,His,pep_intens,pep_rts,special.onlyme);

% draw
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS2_index,MS2_peaks,special);

% Get PSM
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

function His = init_histone()
%%

His.pep_seq = 'KSAPSTGGVKKPHR';
His.mod_short = {'K27me11';
    'K27me21';
    'K27me22';
    'K36me21';
    'K36me22';
    'K27me31';
    'K27me32';
    'K27me33';
    'K27me2K36me11';
    'K27me21K36me11';
    'K27me22K36me11';
    'K27me11K36me2';
    'K27me11K36me21';
    'K27me11K36me22';
    'K27me11K36me1';
    'K27me11K36me11';
    'K27me3K36me11';
    'K27me31K36me11';
    'K27me32K36me11';
    'K27me33K36me11';
    'K27me11K36me3';
    'K27me11K36me31';
    'K27me11K36me32';
    'K27me11K36me33';
    'K27me21K36me2';
    'K27me22K36me2';
    'K27me22K36me21';
    'K27me22K36me22';
    'K27me3K36me21';
    'K27me3K36me22';
    'K27me33K36me2';
    'K27me33K36me21';
    'K27me33K36me22'};
His.mod_type = {'0,pr;1,me11;10,pr;11,pr;';
    '0,pr;1,me21;10,pr;11,pr;';
    '0,pr;1,me22;10,pr;11,pr;';
    '0,pr;1,pr;10,me21;11,pr;';
    '0,pr;1,pr;10,me22;11,pr;';
    '0,pr;1,me31;10,pr;11,pr;';
    '0,pr;1,me32;10,pr;11,pr;';
    '0,pr;1,me33;10,pr;11,pr;';
    '0,pr;1,me2;10,me11;11,pr;';
    '0,pr;1,me21;10,me11;11,pr;';
    '0,pr;1,me22;10,me11;11,pr;';
    '0,pr;1,me11;10,me2;11,pr;';
    '0,pr;1,me11;10,me21;11,pr;';
    '0,pr;1,me11;10,me22;11,pr;';
    '0,pr;1,me11;10,me1;11,pr;';
    '0,pr;1,me11;10,me11;11,pr;';
    '0,pr;1,me3;10,me11;11,pr;';
    '0,pr;1,me31;10,me11;11,pr;';
    '0,pr;1,me32;10,me11;11,pr;';
    '0,pr;1,me33;10,me11;11,pr;';
    '0,pr;1,me11;10,me3;11,pr;';
    '0,pr;1,me11;10,me31;11,pr;';
    '0,pr;1,me11;10,me32;11,pr;';
    '0,pr;1,me11;10,me33;11,pr;';
    '0,pr;1,me21;10,me2;11,pr;';
    '0,pr;1,me22;10,me2;11,pr;';
    '0,pr;1,me22;10,me21;11,pr;';
    '0,pr;1,me22;10,me22;11,pr;';
    '0,pr;1,me3;10,me21;11,pr;';
    '0,pr;1,me3;10,me22;11,pr;';
    '0,pr;1,me33;10,me2;11,pr;';
    '0,pr;1,me33;10,me21;11,pr;';
    '0,pr;1,me33;10,me22;11,pr;'};

His.pep_ch = repmat([2 3 4],length(His.mod_type),1);
%{
His.pep_mz = [816.4574	408.7323
    830.4730	415.7402
    788.4625	394.7349
    802.4781	401.7427
    802.4417	401.7245];
%}
His.pep_mz = calculate_pepmz(His);
His.rt_ref = zeros(length(His.mod_type),1);
His.display = ones(length(His.mod_type),1);

% main ch
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = 2:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end;

function His0 = relocate(cur_outpath,out_filename,His0)
%%

[path1,name1] = fileparts(cur_outpath);
[path2,name2] = fileparts(path1);
out_file1 = fullfile(fullfile(fileparts(path2),name2,name1),[out_filename,'.mat']);
if 0==exist(out_file1,'file')
    fprintf(1,'%s: not exist.\n',out_file1);
    return;
end;
load(out_file1);
His0.rt_ref(1,1) = repmat(auc(3,1),[1,1]);%K27me1
His0.rt_ref(2:3,1) = repmat(auc(4,1),[2,1]);%K27me2
His0.rt_ref(4:5,1) = repmat(auc(5,1),[2,1]);%K36me2
His0.rt_ref(6:8,1) = repmat(auc(6,1),[3,1]);%K27me3
His0.rt_ref(9:11,1) = repmat(auc(8,1),[3,1]);%K27me2K36me1
His0.rt_ref(12:14,1) = repmat(auc(9,1),[3,1]);%K27me1K36me2
His0.rt_ref(15:16,1) = repmat(auc(10,1),[2,1]);%K27me1K36me1
His0.rt_ref(17:20,1) = repmat(auc(11,1),[4,1]);%K27me3K36me1
His0.rt_ref(21:24,1) = repmat(auc(12,1),[4,1]);%K27me1K36me3
His0.rt_ref(25:28,1) = repmat(auc(13,1),[4,1]);%K27me2K36me2
His0.rt_ref(29:33,1) = repmat(auc(14,1),[5,1]);%K27me3K36me2

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,ptol,unitdiff,His,special)%#ok
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

for hno=1:npep
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;