function HH2A_AT08_W6_1_14(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H2A.W.6 peptide MESTGKVKKAFGGR (residues 1-14)
% K at positions 6, 8, 9 in peptide
% K6 modifications: me1, ac
% K8, K9 modifications: ac only
% Isobaric groups: {K6ac,K8ac,K9ac}, {K6acK8ac,K6acK9ac vs K8acK9ac separate}

% check
out_filename = 'HH2A_AT08_W6_1_14';
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

His.pep_seq = 'MESTGKVKKAFGGR';
His.mod_short = {'unmod';
    'K6me1';
    'K6ac';
    'K8ac';
    'K9ac';
    'K6acK8ac';
    'K6acK9ac';
    'K8acK9ac';
    'K6acK8acK9ac'};
His.mod_type = {'0,pr;6,pr;8,pr;9,pr;';
    '0,pr;6,me1;8,pr;9,pr;';
    '0,pr;6,ac;8,pr;9,pr;';
    '0,pr;6,pr;8,ac;9,pr;';
    '0,pr;6,pr;8,pr;9,ac;';
    '0,pr;6,ac;8,ac;9,pr;';
    '0,pr;6,ac;8,pr;9,ac;';
    '0,pr;6,pr;8,ac;9,ac;';
    '0,pr;6,ac;8,ac;9,ac;'};

His.pep_ch = repmat([2 3 4],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [33.90% unmod RT from 0_ref_info.mat
    36% K6me1 (estimate: me1 elutes after unmod)
    32% K6ac (estimate: ac elutes before unmod)
    32% K8ac (isobaric with K6ac)
    32% K9ac (isobaric with K6ac)
    31% K6acK8ac (estimate: double ac)
    31% K6acK9ac (isobaric with K6acK8ac)
    31% K8acK9ac (isobaric with K6acK8ac)
    29];% K6acK8acK9ac (triple ac)
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

% K6me1 (unique mass)
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K6ac/K8ac (isobaric pair, resolved by MS2)
hno = 3;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;

% K9ac (isobaric with K6ac/K8ac)
hno = 5;
if 0==His.rt_ref(3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = His.rt_ref(3);
end;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K6acK8ac/K6acK9ac (isobaric pair, resolved by MS2)
hno = 6;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;

% K8acK9ac (isobaric with K6acK8ac/K6acK9ac)
hno = 8;
if 0==His.rt_ref(6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = His.rt_ref(6);
end;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K6acK8acK9ac (unique mass)
hno = 9;
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

% K6me1 (me1 elutes after unmod)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K6ac (ac elutes before unmod)
hno = 3;
t1 = His.rt_ref(1)-11;
t2 = His.rt_ref(1)-delta;
[rts3,top1_rt3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end;

% K8ac (isobaric with K6ac)
hno = 4;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K9ac (isobaric with K6ac)
hno = 5;
His.rt_ref(hno) = His.rt_ref(3);

% K6acK8ac
hno = 6;
if 0<His.rt_ref(3)
    t1 = His.rt_ref(1)-17;
    t2 = His.rt_ref(3)-delta;
else
    t1 = His.rt_ref(1)-17;
    t2 = His.rt_ref(1)-delta;
end;
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K6acK9ac (isobaric with K6acK8ac)
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K8acK9ac (isobaric with K6acK8ac)
hno = 8;
His.rt_ref(hno) = His.rt_ref(6);

% K6acK8acK9ac
hno = 9;
if 0<His.rt_ref(6)
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(3)
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(3)-delta;
else
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(1)-delta;
end;
[rts9,top1_rt9] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts9)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt9;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K6me1 (me1 elutes after unmod)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K6ac (ac elutes before unmod)
hno = 3;
t1 = His.rt_ref(1)-11;
t2 = His.rt_ref(1)-delta;
[rts3,top1_rt3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end;

% K8ac (isobaric with K6ac)
hno = 4;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K9ac (isobaric with K6ac)
hno = 5;
His.rt_ref(hno) = His.rt_ref(3);

% K6acK8ac
hno = 6;
if 0<His.rt_ref(3)
    t1 = His.rt_ref(1)-17;
    t2 = His.rt_ref(3)-delta;
else
    t1 = His.rt_ref(1)-17;
    t2 = His.rt_ref(1)-delta;
end;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K6acK9ac (isobaric with K6acK8ac)
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K8acK9ac (isobaric with K6acK8ac)
hno = 8;
His.rt_ref(hno) = His.rt_ref(6);

% K6acK8acK9ac
hno = 9;
if 0<His.rt_ref(6)
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(3)
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(3)-delta;
else
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(1)-delta;
end;
[rts9,top1_rt9] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts9)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt9;
end;
