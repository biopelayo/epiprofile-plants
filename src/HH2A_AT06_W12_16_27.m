function HH2A_AT06_W12_16_27(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H2A.W.12 peptide SGGGPKKKPVSR (residues 16-27)
% K at positions 6, 7, 8 in peptide
% K21, K22, K23 modifications: ac only
% Isobaric groups: {K21ac,K22ac,K23ac}, {K21acK22ac,K21acK23ac,K22acK23ac}
% Follows H4_01 acetylation-only pattern

% check
out_filename = 'HH2A_AT06_W12_16_27';
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

His.pep_seq = 'SGGGPKKKPVSR';
His.mod_short = {'unmod';
    'K21ac';
    'K22ac';
    'K23ac';
    'K21acK22ac';
    'K21acK23ac';
    'K22acK23ac';
    'K21acK22acK23ac'};
His.mod_type = {'0,pr;6,pr;7,pr;8,pr;';
    '0,pr;6,ac;7,pr;8,pr;';
    '0,pr;6,pr;7,ac;8,pr;';
    '0,pr;6,pr;7,pr;8,ac;';
    '0,pr;6,ac;7,ac;8,pr;';
    '0,pr;6,ac;7,pr;8,ac;';
    '0,pr;6,pr;7,ac;8,ac;';
    '0,pr;6,ac;7,ac;8,ac;'};

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [20.10% unmod RT from 0_ref_info.mat
    19% K21ac (estimate: ac elutes before unmod)
    19% K22ac (isobaric with K21ac)
    19% K23ac (isobaric with K21ac)
    18% K21acK22ac (estimate: more ac = earlier)
    18% K21acK23ac (isobaric with K21acK22ac)
    18% K22acK23ac (isobaric with K21acK22ac)
    17];% K21acK22acK23ac (triple ac)
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

% K21ac/K22ac (isobaric pair, resolved by MS2)
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;

% K23ac (isobaric with K21ac/K22ac)
hno = 4;
if 0==His.rt_ref(2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = His.rt_ref(2);
end;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K21acK22ac/K21acK23ac (isobaric pair, resolved by MS2)
hno = 5;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;

% K22acK23ac (isobaric with K21acK22ac/K21acK23ac)
hno = 7;
if 0==His.rt_ref(5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = His.rt_ref(5);
end;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K21acK22acK23ac (unique mass)
hno = 8;
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

% K21ac (ac elutes before unmod)
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

% K22ac (isobaric with K21ac)
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno) = His.rt_ref(hno-1);
end;

% K23ac (isobaric with K21ac)
hno = 4;
His.rt_ref(hno) = His.rt_ref(2);

% K21acK22ac
hno = 5;
if 0<His.rt_ref(2)
    t1 = His.rt_ref(1)-17;
    t2 = His.rt_ref(2)-delta;
else
    t1 = His.rt_ref(1)-17;
    t2 = His.rt_ref(1)-delta;
end;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% K21acK23ac (isobaric with K21acK22ac)
hno = 6;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno) = His.rt_ref(hno-1);
end;

% K22acK23ac (isobaric with K21acK22ac)
hno = 7;
His.rt_ref(hno) = His.rt_ref(5);

% K21acK22acK23ac
hno = 8;
if 0<His.rt_ref(5)
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(5)-delta;
elseif 0<His.rt_ref(2)
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(2)-delta;
else
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(1)-delta;
end;
[rts8,top1_rt8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt8;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K21ac (ac elutes before unmod)
hno = 2;
t1 = His.rt_ref(1)-11;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K22ac (isobaric with K21ac)
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno) = His.rt_ref(hno-1);
end;

% K23ac (isobaric with K21ac)
hno = 4;
His.rt_ref(hno) = His.rt_ref(2);

% K21acK22ac
hno = 5;
if 0<His.rt_ref(2)
    t1 = His.rt_ref(1)-17;
    t2 = His.rt_ref(2)-delta;
else
    t1 = His.rt_ref(1)-17;
    t2 = His.rt_ref(1)-delta;
end;
[rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% K21acK23ac (isobaric with K21acK22ac)
hno = 6;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno) = His.rt_ref(hno-1);
end;

% K22acK23ac (isobaric with K21acK22ac)
hno = 7;
His.rt_ref(hno) = His.rt_ref(5);

% K21acK22acK23ac
hno = 8;
if 0<His.rt_ref(5)
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(5)-delta;
elseif 0<His.rt_ref(2)
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(2)-delta;
else
    t1 = His.rt_ref(1)-23;
    t2 = His.rt_ref(1)-delta;
end;
[rts8,top1_rt8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt8;
end;
