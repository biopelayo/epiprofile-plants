function HH2A_AT09_can_36_42(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% Adapted from HH2A_01m1_36_42.m (human KGNYAER)
% AT canonical H2A: KGKYAER (N->K at pos 3, UniProt A0A178WDN2)
% Extra K at position 38 creates new modifiable site
% Note: preceded by K in FASTA (not R); actual Arg-C peptide from upstream R
% would be FLKKGKYAER. This module monitors the 7-mer for consistency with
% the human EpiProfile convention (Yuan lab).
% Generated: 2026-04-16, Sistema Pelamovic

% check
out_filename = 'HH2A_AT09_can_36_42';
fprintf(1,'%s..',out_filename(2:end));
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% init
His = init_histone();

% calculate
unitdiff = 1.0032;
[pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special);

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

function His = init_histone()
%%
% AT canonical H2A: KGKYAER
% vs human: KGNYAER (N->K at residue 38)
% K38 is a NEW modifiable lysine in AT

His.pep_seq = 'KGKYAER';
His.mod_short = {'unmod';
    'K36ac';
    'K38ac';
    'K36acK38ac';
    'Y39ac';
    'Y39pr'};
His.mod_type = {'0,pr;1,pr;3,pr;';
    '0,pr;1,ac;3,pr;';
    '0,pr;1,pr;3,ac;';
    '0,pr;1,ac;3,ac;';
    '0,pr;1,pr;3,pr;4,ac;';
    '0,pr;1,pr;3,pr;4,pr;'};

His.pep_ch = repmat([1 2],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [28.01
    26.33
    26.33
    24.65
    32.1
    36.35];
His.display = ones(length(His.mod_type),1);
His.display(6) = 0;

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

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
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

% K36ac
% K38ac
% K36acK38ac
% Y39ac
% Y39pr
for hno=2:6
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1;
nsplit = 1;

% K36ac
hno = 2;
t1 = His.rt_ref(1)-12;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K38ac
hno = 3;
t1 = His.rt_ref(1)-12;
t2 = His.rt_ref(1)-delta;
[rts3,top1_rt3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end;

% K36acK38ac
hno = 4;
t1 = His.rt_ref(1)-20;
t2 = His.rt_ref(1)-delta;
[rts4,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt4;
end;

% Y39ac
hno = 5;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% Y39pr
hno = 6;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+40;
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

if His.rt_ref(6)>0 && His.rt_ref(5)>His.rt_ref(6)
    % Y39ac
    hno = 5;
    t1 = His.rt_ref(1)+delta;
    t2 = His.rt_ref(6)-delta;
    [rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

    if 1==isempty(rts5)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt5;
    end;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K36ac
hno = 2;
t1 = His.rt_ref(1)-12;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K38ac
hno = 3;
t1 = His.rt_ref(1)-12;
t2 = His.rt_ref(1)-delta;
[rts3,top1_rt3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end;

% K36acK38ac
hno = 4;
t1 = His.rt_ref(1)-20;
t2 = His.rt_ref(1)-delta;
[rts4,top1_rt4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt4;
end;

% Y39ac
hno = 5;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% Y39pr
hno = 6;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+40;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

if His.rt_ref(6)>0 && His.rt_ref(5)>His.rt_ref(6)
    % Y39ac
    hno = 5;
    t1 = His.rt_ref(1)+delta;
    t2 = His.rt_ref(6)-delta;
    [rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

    if 1==isempty(rts5)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt5;
    end;
end;
