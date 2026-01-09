function H3H_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_02_9_17';
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

His.pep_seq = 'KSTGGKAPR';
His.mod_short = {'unmod';
    'K9me1';
    'K9me2';
    'K9me3';
    'K9ac';
    'K14ac';
    'K9me1K14ac';
    'K9me2K14ac';
    'K9me3K14ac';
    'K9acK14ac'};
if 5==special.nsubtype
    His.mod_type = {'0,pr;1,pr;6,pr;';
        '0,pr;1,hme1;6,pr;';
        '0,pr;1,hme2;6,pr;';
        '0,pr;1,hme3;6,pr;';
        '0,pr;1,hac;6,pr;';
        '0,pr;1,pr;6,hac;';
        '0,pr;1,hme1;6,hac;';
        '0,pr;1,hme2;6,hac;';
        '0,pr;1,hme3;6,hac;';
        '0,pr;1,hac;6,hac;'};
else
    His.mod_type = {'0,pr;1,pr;6,pr;';
        '0,pr;1,me1;6,pr;';
        '0,pr;1,me2;6,pr;';
        '0,pr;1,me3;6,pr;';
        '0,pr;1,ac;6,pr;';
        '0,pr;1,pr;6,ac;';
        '0,pr;1,me1;6,ac;';
        '0,pr;1,me2;6,ac;';
        '0,pr;1,me3;6,ac;';
        '0,pr;1,ac;6,ac;'};
end;

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
%{
His.pep_mz = [1069.6	535.3037	357.2049
    1083.6157	542.3115	361.8767
    1041.6051	521.3062	347.8732
    1055.6208	528.3140	352.5451
    1055.5844	528.2958	352.533
    1055.5844	528.2958	352.533
    1069.6000	535.3037	357.2049
    1027.5895	514.2984	343.2013
    1041.6051	521.3062	347.8732
    1041.5687	521.288	347.8611];
%}
His.pep_mz = calculate_pepmzH(His);
His.rt_ref = [23.20
    25.91
    16.59
    16.46
    21.90
    21.91
    24.68
    14.95
    14.86
    20.23];
His.display = ones(length(His.mod_type),1);

His.outpath = cur_outpath;
His.outfile = out_filename;

% main ch
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = 1:npep;
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
% K9me1
% K9me2
% K9me3
for hno=1:4
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

% K9ac/K14ac
if His.rt_ref(6)-His.rt_ref(5)>0.2
    for hno=5:6
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone13(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end;
    end;
else
    hno = 5;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end;
end;

% K9me1K14ac
% K9me2K14ac
% K9me3K14ac
% K9acK14ac
for hno=7:10
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;