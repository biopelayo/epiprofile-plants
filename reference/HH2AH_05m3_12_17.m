function HH2AH_05m3_12_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'HH2A_05m3_12_17';
fprintf(1,'%s..',out_filename(2:end));
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

His.pep_seq = 'AKAKSR';
His.mod_short = {'unmod';
    'K13ac';
    'K15ac';
    'K15me1';
    'K13me1'};
if 5==special.nsubtype
    His.mod_type = {'0,pr;2,pr;4,pr;';
        '0,pr;2,hac;4,pr;';
        '0,pr;2,pr;4,hac;';
        '0,pr;2,pr;4,hme1;';
        '0,pr;2,hme1;4,pr;'};
else
    His.mod_type = {'0,pr;2,pr;4,pr;';
        '0,pr;2,ac;4,pr;';
        '0,pr;2,pr;4,ac;';
        '0,pr;2,pr;4,me1;';
        '0,pr;2,me1;4,pr;'};
end;

His.pep_ch = repmat([1 2],length(His.mod_type),1);
%{
His.pep_mz = [828.4938	414.7505
    814.4781	407.7427
    814.4781	407.7427
    842.5094	421.7584
    842.5094	421.7584];
%}
His.pep_mz = calculate_pepmzH(His);
His.rt_ref = [26.63
    24.56
    24.62
    28.11
    28.22];
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
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K13ac/K15ac
if His.rt_ref(3)-His.rt_ref(2)>0.4
    for hno=2:3
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end;
    end;
else
    hno = 2;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end;
end;

% K15me1/K13me1
if His.rt_ref(5)-His.rt_ref(4)>0.4
    for hno=4:5
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end;
    end;
else
    hno = 4;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end;
end;