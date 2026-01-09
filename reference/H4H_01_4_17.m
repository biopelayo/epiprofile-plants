function H4H_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H4_01_4_17';
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
% change order
[His,pep_rts,pep_intens,mono_isointens] = change_order(His,pep_rts,pep_intens,mono_isointens);
%

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

His.pep_seq = 'GKGGKGLGKGGAKR';
His.mod_short = {'unmod';
    'K5ac';
    'K8ac';
    'K12ac';
    'K16ac';
    'K5acK8ac';
    'K5acK12ac';
    'K5acK16ac';
    'K8acK12ac';
    'K8acK16ac';
    'K12acK16ac';
    'K8acK12acK16ac';
    'K5acK12acK16ac';
    'K5acK8acK16ac';
    'K5acK8acK12ac';
    'K5acK8acK12acK16ac'};
if 5==special.nsubtype
    His.mod_type = {'0,pr;2,pr;5,pr;9,pr;13,pr;';
        '0,pr;2,hac;5,pr;9,pr;13,pr;';
        '0,pr;2,pr;5,hac;9,pr;13,pr;';
        '0,pr;2,pr;5,pr;9,hac;13,pr;';
        '0,pr;2,pr;5,pr;9,pr;13,hac;';
        '0,pr;2,hac;5,hac;9,pr;13,pr;';
        '0,pr;2,hac;5,pr;9,hac;13,pr;';
        '0,pr;2,hac;5,pr;9,pr;13,hac;';
        '0,pr;2,pr;5,hac;9,hac;13,pr;';
        '0,pr;2,pr;5,hac;9,pr;13,hac;';
        '0,pr;2,pr;5,pr;9,hac;13,hac;';
        '0,pr;2,pr;5,hac;9,hac;13,hac;';
        '0,pr;2,hac;5,pr;9,hac;13,hac;';
        '0,pr;2,hac;5,hac;9,pr;13,hac;';
        '0,pr;2,hac;5,hac;9,hac;13,pr;';
        '0,pr;2,hac;5,hac;9,hac;13,hac;'};
else
    His.mod_type = {'0,pr;2,pr;5,pr;9,pr;13,pr;';
        '0,pr;2,ac;5,pr;9,pr;13,pr;';
        '0,pr;2,pr;5,ac;9,pr;13,pr;';
        '0,pr;2,pr;5,pr;9,ac;13,pr;';
        '0,pr;2,pr;5,pr;9,pr;13,ac;';
        '0,pr;2,ac;5,ac;9,pr;13,pr;';
        '0,pr;2,ac;5,pr;9,ac;13,pr;';
        '0,pr;2,ac;5,pr;9,pr;13,ac;';
        '0,pr;2,pr;5,ac;9,ac;13,pr;';
        '0,pr;2,pr;5,ac;9,pr;13,ac;';
        '0,pr;2,pr;5,pr;9,ac;13,ac;';
        '0,pr;2,pr;5,ac;9,ac;13,ac;';
        '0,pr;2,ac;5,pr;9,ac;13,ac;';
        '0,pr;2,ac;5,ac;9,pr;13,ac;';
        '0,pr;2,ac;5,ac;9,ac;13,pr;';
        '0,pr;2,ac;5,ac;9,ac;13,ac;'};
end;

His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);
%{
His.pep_mz = [1550.9013	775.9543	517.6386	388.4808
    1536.8857	768.9465	512.9667	384.9769
    1536.8857	768.9465	512.9667	384.9769
    1536.8857	768.9465	512.9667	384.9769
    1536.8857	768.9465	512.9667	384.9769
    1522.8700	761.9386	508.2949	381.473
    1522.8700	761.9386	508.2949	381.473
    1522.8700	761.9386	508.2949	381.473
    1522.8700	761.9386	508.2949	381.473
    1522.8700	761.9386	508.2949	381.473
    1522.8700	761.9386	508.2949	381.473
    1508.8544	754.9308	503.6230	377.969
    1508.8544	754.9308	503.6230	377.969
    1508.8544	754.9308	503.6230	377.969
    1508.8544	754.9308	503.6230	377.969
    1494.8387	747.9230	498.9511	374.4651];
%}
His.pep_mz = calculate_pepmzH(His);
His.rt_ref = [30.08
    28.80
    28.80
    28.79
    29.11
    27.78
    27.72
    27.94
    27.78
    27.94
    27.88
    26.72
    26.65
    26.76
    26.55
    25.35];
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

% K5ac/K8ac/K12ac/K16ac
%{
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
%}
%
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone4(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+3,1:ncharge) = cur_rts(1:4,:);
    pep_intens(hno:hno+3,1:ncharge) = cur_intens(1:4,:);
    mono_isointens(1:num_MS1,hno:hno+3) = cur_mono_isointens(:,1:4);
end;
%}

% K5acK8ac/K5acK12ac/K5acK16ac/K8acK12ac/K8acK16ac/K12acK16ac
%{
hno = 6;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
%}
%
hno = 6;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone6(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+5,1:ncharge) = cur_rts(1:6,:);
    pep_intens(hno:hno+5,1:ncharge) = cur_intens(1:6,:);
    mono_isointens(1:num_MS1,hno:hno+5) = cur_mono_isointens(:,1:6);
end;
%}

% K8acK12acK16ac/K5acK12acK16ac/K5acK8acK16ac/K5acK8acK12ac
%{
hno = 12;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
%}
%
hno = 12;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone4(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+3,1:ncharge) = cur_rts(1:4,:);
    pep_intens(hno:hno+3,1:ncharge) = cur_intens(1:4,:);
    mono_isointens(1:num_MS1,hno:hno+3) = cur_mono_isointens(:,1:4);
end;
%}

% K5acK8acK12acK16ac
hno = 16;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

function [His,pep_rts,pep_intens,mono_isointens] = change_order(His,pep_rts,pep_intens,mono_isointens)
%% 12,13,14,15->15,14,13,12

% copy
b1 = His.mod_short{15};
b2 = His.mod_short{14};
b3 = His.mod_short{13};
b4 = His.mod_short{12};

c1 = His.mod_type{15};
c2 = His.mod_type{14};
c3 = His.mod_type{13};
c4 = His.mod_type{12};

d1 = His.rt_ref(15);
d2 = His.rt_ref(14);
d3 = His.rt_ref(13);
d4 = His.rt_ref(12);

e1 = pep_rts(15,:);
e2 = pep_rts(14,:);
e3 = pep_rts(13,:);
e4 = pep_rts(12,:);

f1 = pep_intens(15,:);
f2 = pep_intens(14,:);
f3 = pep_intens(13,:);
f4 = pep_intens(12,:);

g1 = mono_isointens(:,15);
g2 = mono_isointens(:,14);
g3 = mono_isointens(:,13);
g4 = mono_isointens(:,12);

% paste
His.mod_short{12} = b1;
His.mod_short{13} = b2;
His.mod_short{14} = b3;
His.mod_short{15} = b4;

His.mod_type{12} = c1;
His.mod_type{13} = c2;
His.mod_type{14} = c3;
His.mod_type{15} = c4;

His.rt_ref(12) = d1;
His.rt_ref(13) = d2;
His.rt_ref(14) = d3;
His.rt_ref(15) = d4;

pep_rts(12,:) = e1;
pep_rts(13,:) = e2;
pep_rts(14,:) = e3;
pep_rts(15,:) = e4;

pep_intens(12,:) = f1;
pep_intens(13,:) = f2;
pep_intens(14,:) = f3;
pep_intens(15,:) = f4;

mono_isointens(:,12) = g1;
mono_isointens(:,13) = g2;
mono_isointens(:,14) = g3;
mono_isointens(:,15) = g4;