function H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H4_01_4_17';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% init
His = init_histone(cur_outpath,out_filename);

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

function His = init_histone(cur_outpath,out_filename)
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
His.pep_mz = calculate_pepmz(His);
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

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special)
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% unmod
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        [His.rt_ref(1),special.ndebug] = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
    else
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
        else
            delta = 5;
            t1 = His.rt_ref(1)-delta;
            t2 = His.rt_ref(1)+delta;
        end;
        hno = 1;% unmod
        [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,1,t1,t2,nhmass);%#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end;
    end;
end;

hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% calibrate the rt_ref
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end;
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

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1;
nsplit = 1;

% K5ac
hno = 2;
t1 = His.rt_ref(1)-11;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K8ac/K12ac/K16ac
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:5) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:5) = His.rt_ref(hno-1);
end;

% K5acK8ac
hno = 6;
t1 = His.rt_ref(1)-17;
if 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K5acK12ac/K5acK16ac/K8acK12ac/K8acK16ac/K12acK16ac
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:11) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:11) = His.rt_ref(hno-1);
end;

% K8acK12acK16ac
hno = 12;
t1 = His.rt_ref(1)-23;
if 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts12,top1_rt12] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts12)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt12;
end;

% K5acK12acK16ac/K5acK8acK16ac/K5acK8acK12ac
hno = 13;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:15) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:15) = His.rt_ref(hno-1);
end;

% K5acK8acK12acK16ac
hno = 16;
t1 = His.rt_ref(1)-28;
if 0<His.rt_ref(12)
    t2 = His.rt_ref(12)-delta;
elseif 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts16,top1_rt16] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts16)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt16;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K5ac
hno = 2;
t1 = His.rt_ref(1)-11;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K8ac/K12ac/K16ac
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:5) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:5) = His.rt_ref(hno-1);
end;

% K5acK8ac
hno = 6;
t1 = His.rt_ref(1)-17;
if 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K5acK12ac/K5acK16ac/K8acK12ac/K8acK16ac/K12acK16ac
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:11) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:11) = His.rt_ref(hno-1);
end;

% K8acK12acK16ac
hno = 12;
t1 = His.rt_ref(1)-23;
if 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts12,top1_rt12] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts12)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt12;
end;

% K5acK12acK16ac/K5acK8acK16ac/K5acK8acK12ac
hno = 13;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:15) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:15) = His.rt_ref(hno-1);
end;

% K5acK8acK12acK16ac
hno = 16;
t1 = His.rt_ref(1)-28;
if 0<His.rt_ref(12)
    t2 = His.rt_ref(12)-delta;
elseif 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts16,top1_rt16] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts16)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt16;
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