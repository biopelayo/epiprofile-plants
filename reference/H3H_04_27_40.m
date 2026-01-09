function H3H_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_04_27_40';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% init
His = init_histone(cur_outpath,out_filename,special);

% relocate
His = relocate(cur_outpath,out_filename,His);

% calculate
unitdiff = 1.0032;
Mods = GetMods();
[pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

% output
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% draw
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS2_index,MS2_peaks,special);

% Get PSM
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

function His = init_histone(cur_outpath,out_filename,special)
%%

His.pep_seq = 'KSAPATGGVKKPHR';
His.mod_short = {'unmod';
    'K36me1';
    'K27me1';
    'K27me2';
    'K36me2';
    'K27me3';
    'K36me3';
    'K27me2K36me1';
    'K27me1K36me2';
    'K27me1K36me1';
    'K27me3K36me1';
    'K27me1K36me3';
    'K27me2K36me2';
    'K27me3K36me2';
    'K27ac'};
if 5==special.nsubtype
    His.mod_type = {'0,pr;1,pr;10,pr;11,pr;';
        '0,pr;1,pr;10,hme1;11,pr;';
        '0,pr;1,hme1;10,pr;11,pr;';
        '0,pr;1,hme2;10,pr;11,pr;';
        '0,pr;1,pr;10,hme2;11,pr;';
        '0,pr;1,hme3;10,pr;11,pr;';
        '0,pr;1,pr;10,hme3;11,pr;';
        '0,pr;1,hme2;10,hme1;11,pr;';
        '0,pr;1,hme1;10,hme2;11,pr;';
        '0,pr;1,hme1;10,hme1;11,pr;';
        '0,pr;1,hme3;10,hme1;11,pr;';
        '0,pr;1,hme1;10,hme3;11,pr;';
        '0,pr;1,hme2;10,hme2;11,pr;';
        '0,pr;1,hme3;10,hme2;11,pr;';
        '0,pr;1,hac;10,pr;11,pr;'};
else
    His.mod_type = {'0,pr;1,pr;10,pr;11,pr;';
        '0,pr;1,pr;10,me1;11,pr;';
        '0,pr;1,me1;10,pr;11,pr;';
        '0,pr;1,me2;10,pr;11,pr;';
        '0,pr;1,pr;10,me2;11,pr;';
        '0,pr;1,me3;10,pr;11,pr;';
        '0,pr;1,pr;10,me3;11,pr;';
        '0,pr;1,me2;10,me1;11,pr;';
        '0,pr;1,me1;10,me2;11,pr;';
        '0,pr;1,me1;10,me1;11,pr;';
        '0,pr;1,me3;10,me1;11,pr;';
        '0,pr;1,me1;10,me3;11,pr;';
        '0,pr;1,me2;10,me2;11,pr;';
        '0,pr;1,me3;10,me2;11,pr;';
        '0,pr;1,ac;10,pr;11,pr;'};
end;

His.pep_ch = repmat([2 3 4],length(His.mod_type),1);
%{
His.pep_mz = [829.4728	553.3177	415.2401
    836.4807	557.9895	418.7440
    836.4807	557.9895	418.7440
    815.4754	543.9860	408.2413
    815.4754	543.9860	408.2413
    822.4832	548.6579	411.7452
    822.4832	548.6579	411.7452
    822.4832	548.6579	411.7452
    822.4832	548.6579	411.7452
    843.4885	562.6614	422.2479
    829.4910	553.3298	415.2492
    829.4910	553.3298	415.2492
    801.4779	534.6544	401.2426
    808.4858	539.3263	404.7465
    822.4650	548.6458	411.7361];
%}
His.pep_mz = calculate_pepmzH(His);
His.rt_ref = [26.08
    27.81
    27.94
    21.40
    22.68
    21.40
    22.66
    22.64
    24.44
    29.15
    22.64
    24.44
    18.29
    18.20
    24.99];
His.display = ones(length(His.mod_type),1);

His.outpath = cur_outpath;
His.outfile = out_filename;

% main ch
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);%#ok
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = [4 5 6 7 8 9 11 12 13 14];%1:npep;
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
for ino=1:length(His0.rt_ref)
    His0.rt_ref(ino) = auc(ino,1);
end;

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special)
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% unmod
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K36me1/K27me1
if His.rt_ref(3)-His.rt_ref(2)>0.4
    for hno=2:3
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end;
    end;
else
    hno = 2;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone22(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end;
end;

% K27me2
% K36me2
% K27me3
for hno=[4 6]
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;
hno = 5;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K36me3/K27me2K36me1
if His.rt_ref(7)-His.rt_ref(8)>0.4
    for hno=7:8
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end;
    end;
else
    hno = 7;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end;
end;

% K27me1K36me2
% K27me1K36me1
% K27me3K36me1
% K27me1K36me3
% K27me2K36me2
for hno=9:13
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

% K27me3K36me2
hno = 14;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K27ac
hno = 15;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;