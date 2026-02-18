function H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_04_27_40';
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
    'K27me2K36me3';
    'K27me3K36me3';
    'K27ac'};
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
    '0,pr;1,me2;10,me3;11,pr;';
    '0,pr;1,me3;10,me3;11,pr;';
    '0,pr;1,ac;10,pr;11,pr;'};

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
    808.4858	539.3263	404.7465
    815.4936	543.9981	408.2504
    822.4650	548.6458	411.7361];
%}
His.pep_mz = calculate_pepmz(His);
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
    18.34
    18.20
    24.99];
His.display = ones(length(His.mod_type),1);
His.display([15 16]) = 0;

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
    tune = [4 5 6 7 8 9 11 12 13 14 15 16];%1:npep;
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
        if 2==special.ndebug
            hno = 1;% unmod
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);
            hno = 2;% K27me1
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2);%#ok
            hno = 4;% K27me2
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts4,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);%#ok
            if 0==isempty(top1_rt2) && 0==isempty(top1_rt4) && top1_rt2>top1_rt4
                if 1==isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt4+(top1_rt2-top1_rt4)*0.7;
                else
                    p = find(rts1>top1_rt4+(top1_rt2-top1_rt4)*0.5 & rts1<top1_rt2);
                    if 0==isempty(p)
                        [tmp,pp] = max(inten_sum1(p));%#ok
                        His.rt_ref(1) = rts1(p(pp));
                    end;
                end;
            else
                if 0==isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt1;
                end;
            end;
        end;
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

% K36me1/K27me1
if His.rt_ref(3)-His.rt_ref(2)>0.4
    for hno=2:3
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone10(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
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
hno = 7;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
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

% K27me3K36me2/K27me2K36me3
%
hno = 14;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
%}
%{
hno = 14;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;
%}

% K27me3K36me3
%{
hno = 16;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
%}

% K27ac
hno = 17;
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

% K36me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+15;
[rts2,top1_rt2,inten_sum2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2);

[tmp_sum,ix] = sort(inten_sum2,'descend');
tmp_rts = rts2(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
    His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end;
    
    % K27me1
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

% K27me2,K27me3,K27me3K36me1 (4,6,11)
hno = 4;
t1 = His.rt_ref(1)-35;
t2 = His.rt_ref(1)-delta;
[rts4,top1_rt4,inten_sum4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

hno = 6;
t1 = His.rt_ref(1)-35;
t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6,inten_sum6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);%#ok

hno = 11;
t1 = His.rt_ref(1)-35;
t2 = His.rt_ref(1)-delta;
[rts11,top1_rt11,inten_sum11] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);%#ok

xt = find_triple(rts4,top1_rt4,rts6,rts11,inten_sum4,inten_sum6,inten_sum11);

% K27me2K36me2,K27me3K36me2,K27me3K36me3 (13,14,16)
hno = 13;
t1 = His.rt_ref(1)-50;
if 0==xt(1,1)
    t2 = His.rt_ref(1)-4;
else
    t2 = xt(1,1)-delta-0.4;
end;
[rts13,top1_rt13,inten_sum13] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

hno = 14;
t1 = His.rt_ref(1)-50;
if 0==xt(1,1)
    t2 = His.rt_ref(1)-4;
else
    t2 = xt(1,1)-delta-0.4;
end;
[rts14,top1_rt14,inten_sum14] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

hno = 16;
t1 = His.rt_ref(1)-50;
if 0==xt(1,1)
    t2 = His.rt_ref(1)-4;
else
    t2 = xt(1,1)-delta-0.4;
end;
[rts16,top1_rt16,inten_sum16] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);%#ok

if 0==isempty(inten_sum13) || 0==isempty(inten_sum14)
    rt_array = [0 0];
    in_array = [0 0];
    if 0==isempty(inten_sum13)
        rt_array(1) = top1_rt13;
        in_array(1) = max(inten_sum13);
    end;
    if 0==isempty(inten_sum14)
        rt_array(2) = top1_rt14;
        in_array(2) = max(inten_sum14);
    end;
    [tmp,ix] = max(in_array);%#ok
    xt4 = rt_array(ix);
else
    xt4 = 0;
end;

% K27me2
hno = 4;
His.rt_ref(hno) = xt(1,1);

% K36me2
hno = 5;
His.rt_ref(hno) = xt(1,2);

% K27me3
hno = 6;
His.rt_ref(hno) = xt(2,1);

% K36me3
hno = 7;
His.rt_ref(hno) = xt(1,3);

% K27me2K36me1
hno = 8;
His.rt_ref(hno) = xt(2,2);

% K27me1K36me2
hno = 9;
His.rt_ref(hno) = xt(2,3);

% K27me1K36me1
hno = 10;
if 0==His.rt_ref(3)
    t1 = His.rt_ref(1)+delta;
else
    t1 = His.rt_ref(3)+delta;
end;
if 0==His.rt_ref(3)
    t2 = His.rt_ref(1)+18;
else
    t2 = His.rt_ref(3)+(His.rt_ref(3)-His.rt_ref(1))+2;
end;
[rts10,top1_rt10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts10)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt10;
end;

% K27me3K36me1
hno = 11;
His.rt_ref(hno) = xt(3,2);

% K27me1K36me3
hno = 12;
His.rt_ref(hno) = xt(3,3);

% K27me2K36me2
hno = 13;
ix = find(abs(rts13-xt4)<=delta+0.4);
if 0==isempty(ix)
    [tmp,id] = max(inten_sum13(ix));%#ok
    His.rt_ref(hno) = rts13(ix(id));
else
    His.rt_ref(hno) = 0;
end;

% K27me3K36me2
hno = 14;
old_t = His.rt_ref(hno);
ix = find(abs(rts14-xt4)<=delta+0.4);
if 0==isempty(ix)
    [tmp,id] = max(inten_sum14(ix));%#ok
    His.rt_ref(hno) = rts14(ix(id));
else
    His.rt_ref(hno) = 0;
end;

% K27me2K36me3
hno = 15;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K27me3K36me3
hno = 16;
ix = find(abs(rts16-xt4)<=delta+0.4);
if 0==isempty(ix)
    [tmp,id] = max(inten_sum16(ix));%#ok
    His.rt_ref(hno) = rts16(ix(id));
else
    His.rt_ref(hno) = 0;
end;

% K27ac
hno = 17;
if 0==xt(2,3)
    t1 = His.rt_ref(1)-10;
else
    t1 = xt(2,3)+delta;
end;
t2 = His.rt_ref(1)-delta;
[rts17,top1_rt17] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts17)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt17;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K36me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+15;
[rts2,top1_rt2,inten_sum2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,0,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K27me1
hno = 3;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+15;
[rts3,top1_rt3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,0,t1,t2,nhmass);

if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end;

if 0==isempty(rts2) && 0==isempty(rts3)
    if top1_rt2==top1_rt3
        hno = 2;
        [tmp_sum,ix] = sort(inten_sum2,'descend');
        tmp_rts = rts2(ix);
        if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
            His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
            His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
        else
            old_t = His.rt_ref(hno);
            if 1==isempty(rts2)
                His.rt_ref(hno) = 0;
            else
                His.rt_ref(hno) = top1_rt2;
            end;
            
            % K27me1
            hno = 3;
            if 0==His.rt_ref(hno-1)
                His.rt_ref(hno) = 0;
            elseif old_t~= His.rt_ref(hno-1);
                d = His.rt_ref(hno-1) - old_t;
                His.rt_ref(hno) = His.rt_ref(hno) + d;
            end;
        end;
    elseif top1_rt2>top1_rt3 && abs(top1_rt2-top1_rt3)<2
        hno = 2;
        His.rt_ref(hno) = top1_rt3;
        His.rt_ref(hno+1) = top1_rt2;
    end;
end;

% K27me2,K27me3,K27me3K36me1 (4,6,11)
hno = 4;
t1 = His.rt_ref(1)-35;
t2 = His.rt_ref(1)-delta;
[rts4,top1_rt4,inten_sum4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

hno = 6;
t1 = His.rt_ref(1)-35;
t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6,inten_sum6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);%#ok

hno = 11;
t1 = His.rt_ref(1)-35;
t2 = His.rt_ref(1)-delta;
[rts11,top1_rt11,inten_sum11] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);%#ok

xt = find_triple(rts4,top1_rt4,rts6,rts11,inten_sum4,inten_sum6,inten_sum11);

% K27me2K36me2,K27me3K36me2,K27me3K36me3 (13,14,16)
hno = 13;
t1 = His.rt_ref(1)-50;
if 0==xt(1,1)
    t2 = His.rt_ref(1)-4;
else
    t2 = xt(1,1)-delta-0.4;
end;
[rts13,top1_rt13,inten_sum13,top1_inten_sum13] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

hno = 14;
t1 = His.rt_ref(1)-50;
if 0==xt(1,1)
    t2 = His.rt_ref(1)-4;
else
    t2 = xt(1,1)-delta-0.4;
end;
[rts14,top1_rt14,inten_sum14,top1_inten_sum14] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

hno = 16;
t1 = His.rt_ref(1)-50;
if 0==xt(1,1)
    t2 = His.rt_ref(1)-4;
else
    t2 = xt(1,1)-delta-0.4;
end;
[rts16,top1_rt16,inten_sum16] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);%#ok

if 0==isempty(inten_sum13) || 0==isempty(inten_sum14)
    rt_array = [0 0];
    in_array = [0 0];
    if 0==isempty(inten_sum13)
        rt_array(1) = top1_rt13;
        in_array(1) = top1_inten_sum13;
    end;
    if 0==isempty(inten_sum14)
        rt_array(2) = top1_rt14;
        in_array(2) = top1_inten_sum14;
    end;
    [tmp,ix] = max(in_array);%#ok
    xt4 = rt_array(ix);
else
    xt4 = 0;
end;

% K27me2
hno = 4;
His.rt_ref(hno) = xt(1,1);

% K36me2
hno = 5;
His.rt_ref(hno) = xt(1,2);

% K27me3
hno = 6;
His.rt_ref(hno) = xt(2,1);

% K36me3
hno = 7;
His.rt_ref(hno) = xt(1,3);

% K27me2K36me1
hno = 8;
His.rt_ref(hno) = xt(2,2);

% K27me1K36me2
hno = 9;
His.rt_ref(hno) = xt(2,3);

% K27me1K36me1
hno = 10;
if 0==His.rt_ref(3)
    t1 = His.rt_ref(1)+delta;
else
    t1 = His.rt_ref(3)+delta;
end;
if 0==His.rt_ref(3)
    t2 = His.rt_ref(1)+18;
else
    t2 = His.rt_ref(3)+(His.rt_ref(3)-His.rt_ref(1))+2;
end;
[rts10,top1_rt10] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts10)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt10;
end;

% K27me3K36me1
hno = 11;
His.rt_ref(hno) = xt(3,2);

% K27me1K36me3
hno = 12;
His.rt_ref(hno) = xt(3,3);

% K27me2K36me2
hno = 13;
ix = find(abs(rts13-xt4)<=delta+0.4);
if 0==isempty(ix)
    [tmp,id] = max(inten_sum13(ix));%#ok
    His.rt_ref(hno) = rts13(ix(id));
else
    His.rt_ref(hno) = 0;
end;

% K27me3K36me2
hno = 14;
old_t = His.rt_ref(hno);
ix = find(abs(rts14-xt4)<=delta+0.4);
if 0==isempty(ix)
    [tmp,id] = max(inten_sum14(ix));%#ok
    His.rt_ref(hno) = rts14(ix(id));
else
    His.rt_ref(hno) = 0;
end;

% K27me2K36me3
hno = 15;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K27me3K36me3
hno = 16;
ix = find(abs(rts16-xt4)<=delta+0.4);
if 0==isempty(ix)
    [tmp,id] = max(inten_sum16(ix));%#ok
    His.rt_ref(hno) = rts16(ix(id));
else
    His.rt_ref(hno) = 0;
end;

% K27ac
hno = 17;
if 0==xt(2,3)
    t1 = His.rt_ref(1)-10;
else
    t1 = xt(2,3)+delta;
end;
t2 = His.rt_ref(1)-delta;
[rts17,top1_rt17] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts17)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt17;
end;