function HH2A_AT05_Xa_41_50(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H2A.Xa peptide FLKSGKYAER (residues 41-50)
% K at positions 3, 6 in peptide
% K43 modifications: me1, me2, me3, ac
% K46 modifications: ac only
% Isobaric pairs: {K43ac,K46ac}, {K43me1K46ac unique}, {K43acK46ac unique}

% check
out_filename = 'HH2A_AT05_Xa_41_50';
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

His.pep_seq = 'FLKSGKYAER';
His.mod_short = {'unmod';
    'K43me1';
    'K43me2';
    'K43me3';
    'K43ac';
    'K46ac';
    'K43me1K46ac';
    'K43acK46ac'};
His.mod_type = {'0,pr;3,pr;6,pr;';
    '0,pr;3,me1;6,pr;';
    '0,pr;3,me2;6,pr;';
    '0,pr;3,me3;6,pr;';
    '0,pr;3,ac;6,pr;';
    '0,pr;3,pr;6,ac;';
    '0,pr;3,me1;6,ac;';
    '0,pr;3,ac;6,ac;'};

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [37.42% unmod RT from 0_ref_info.mat
    40% K43me1 (estimate: me1 elutes after unmod)
    30% K43me2 (estimate: me2 elutes before unmod)
    30% K43me3 (estimate: me3 co-elutes with me2)
    36% K43ac (estimate: between me2 and unmod)
    36% K46ac (isobaric with K43ac)
    39% K43me1K46ac (estimate: me1+ac combo)
    34];% K43acK46ac (estimate: double ac)
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

% K43me1
% K43me2
% K43me3
for hno=2:4
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

% K43ac/K46ac (isobaric pair, resolved by MS2)
hno = 5;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;

% K43me1K46ac (unique mass)
hno = 7;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K43acK46ac (unique mass)
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

% K43me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K43me2
hno = 3;
t1 = 4;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K43me3
hno = 4;
t1 = 4;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[His.rt_ref(3),His.rt_ref(4)] = find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3,rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K43ac
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% K46ac (isobaric with K43ac)
hno = 6;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K43me1K46ac
hno = 7;
t1 = His.rt_ref(1)+delta;
if 0==His.rt_ref(2)
    t2 = His.rt_ref(1)+9;
else
    t2 = His.rt_ref(2)-delta;
end;
[rts7,top1_rt7] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% K43acK46ac
hno = 8;
if 0==His.rt_ref(3)
    t1 = His.rt_ref(1)-18;
else
    t1 = His.rt_ref(3)+delta;
end;
if 0==His.rt_ref(5)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(5)-delta;
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

% K43me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K43me2
hno = 3;
t1 = 4;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K43me3
hno = 4;
t1 = 4;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[His.rt_ref(3),His.rt_ref(4)] = find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3,rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K46ac (resolve first, then K43ac uses same delta)
hno = 6;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K43ac (isobaric with K46ac)
hno = 5;
if 0==His.rt_ref(hno+1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno+1);
    d = His.rt_ref(hno+1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K43me1K46ac
hno = 7;
t1 = His.rt_ref(1)+delta;
if 0==His.rt_ref(2)
    t2 = His.rt_ref(1)+9;
else
    t2 = His.rt_ref(2)-delta;
end;
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% K43acK46ac
hno = 8;
if 0==His.rt_ref(3)
    t1 = His.rt_ref(1)-18;
else
    t1 = His.rt_ref(3)+delta;
end;
if 0==His.rt_ref(5)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(5)-delta;
end;
[rts8,top1_rt8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt8;
end;
