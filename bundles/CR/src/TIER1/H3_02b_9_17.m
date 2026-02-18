function H3_02b_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_02b_9_17';
%fprintf(1,'%s..',out_filename);
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

His.pep_seq = 'KSTGGKAPR';
His.mod_short = {'S10ac';
    'S10pr';
    'K9me1S10ac';
    'K9me2S10ac';
    'K9me3S10ac';
    'K9acS10ac';
    'S10acK14ac';
    'S10prK14ac';
    'K9me2S10acK14ac';
    'K9me3S10acK14ac';
    'K9acS10acK14ac'};
His.mod_type = {'0,pr;1,pr;2,ac;6,pr;';
    '0,pr;1,pr;2,pr;6,pr;';
    '0,pr;1,me1;2,ac;6,pr;';
    '0,pr;1,me2;2,ac;6,pr;';
    '0,pr;1,me3;2,ac;6,pr;';
    '0,pr;1,ac;2,ac;6,pr;';
    '0,pr;1,pr;2,ac;6,ac;';
    '0,pr;1,pr;2,pr;6,ac;';
    '0,pr;1,me2;2,ac;6,ac;';
    '0,pr;1,me3;2,ac;6,ac;';
    '0,pr;1,ac;2,ac;6,ac;'};

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
%{
His.pep_mz = [1111.6106	556.3089	371.2084
    1125.6263	563.3168	375.8803
    1125.6263	563.3168	375.8803
    1083.6157	542.3115	361.8767
    1097.6313	549.3193	366.5486
    1097.5950	549.3011	366.5365
    1097.5950	549.3011	366.5365
    1111.6106	556.3089	371.2084
    1069.6000	535.3037	357.2049
    1083.6157	542.3115	361.8767
    1083.5793	542.2933	361.8646];
%}
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [23.20
    25.91
    25.95
    16.59
    16.46
    21.90
    21.91
    24.68
    14.95
    14.86
    20.23];
His.display = ones(length(His.mod_type),1);
His.display([2 8]) = 0;

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
        if 2==special.ndebug
            hno = 1;% unmod
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);
            hno = 2;% K9me1
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);%#ok
            if 0==isempty(top1_rt2)
                if 1==isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt2-3;
                else
                    p = find(rts1>top1_rt2-16 & rts1<top1_rt2);
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

% S10pr/K9me1
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

% K9me2
% K9me3
for hno=4:5
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

% K9ac/K14ac
hno = 6;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;
% K9me1K14ac

% K9me2K14ac
% K9me3K14ac
% K9acK14ac
for hno=[9 11]
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;
hno = 10;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% correct K9me1K14ac
delta = 0.3;
newhno = 1;% unmod
rehno = 2;% K9me1
trihno = 8;% K9me1K14ac
% switch His of K9me1 and K9me1K14ac
tmp1 = His.mod_short{rehno};
tmp2 = His.mod_type{rehno};
tmp3 = His.pep_mz(rehno,1:ncharge);
tmp4 = His.rt_ref(rehno);
His.mod_short{rehno} = His.mod_short{trihno};
His.mod_type{rehno} = His.mod_type{trihno};
His.pep_mz(rehno,1:ncharge) = His.pep_mz(trihno,1:ncharge);
His.rt_ref(rehno) = His.rt_ref(trihno);
His.mod_short{trihno} = tmp1;
His.mod_type{trihno} = tmp2;
His.pep_mz(trihno,1:ncharge) = tmp3;
His.rt_ref(trihno) = tmp4;
His.rt_ref(rehno) = His.rt_ref(newhno)+delta;
% unmod/K9me1K14ac
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(newhno,1:ncharge) = cur_rts(1,:);
    pep_intens(newhno,1:ncharge) = cur_intens(1,:);
    mono_isointens(1:num_MS1,newhno) = cur_mono_isointens(:,1);
    pep_rts(trihno,1:ncharge) = cur_rts(2,:);
    pep_intens(trihno,1:ncharge) = cur_intens(2,:);
    mono_isointens(1:num_MS1,trihno) = cur_mono_isointens(:,2);
end;

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1;
nsplit = 1;

% S10pr
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2,inten_sum2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

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
    
    % K9me1
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

% K9me2
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K9me3
hno = 5;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);%#ok

[His.rt_ref(5),His.rt_ref(4)] = find_pair_new(top1_rt5,rts4,top1_rt4,inten_sum4,0);

% K9ac
hno = 6;
t1 = (His.rt_ref(1)+His.rt_ref(4))/2;
t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K14ac
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K9me1K14ac
hno = 8;
t1 = His.rt_ref(1)+delta;
if 0==His.rt_ref(2)
    t2 = His.rt_ref(1)+9;
else
    t2 = His.rt_ref(2)-delta;
end;
[rts8,top1_rt8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt8;
end;

% K9me2K14ac
hno = 9;
t1 = 6;
if 0==His.rt_ref(4)
    t2 = His.rt_ref(1)-4;
else
    t2 = His.rt_ref(4)-delta;
end;
[rts9,top1_rt9,inten_sum9] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K9me3K14ac
hno = 10;
t1 = 6;
if 0==His.rt_ref(4)
    t2 = His.rt_ref(1)-4;
else
    t2 = His.rt_ref(4)-delta;
end;
[rts10,top1_rt10,inten_sum10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 0==isempty(top1_rt10) && top1_rt10>0 && His.rt_ref(4)>0 && top1_rt10<His.rt_ref(4)
    [His.rt_ref(10),His.rt_ref(9)] = find_pair_new(top1_rt10,rts9,top1_rt9,inten_sum9,0);
else
    [His.rt_ref(9),His.rt_ref(10)] = find_pair_new(top1_rt9,rts10,top1_rt10,inten_sum10,1);
end;

% K9acK14ac
hno = 11;
if 0==His.rt_ref(4)
    t1 = His.rt_ref(1)-18;
else
    t1 = His.rt_ref(4)+delta;
end;
if 0==His.rt_ref(6)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(6)-delta;
end;
[rts11,top1_rt11] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts11)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt11;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% S10pr
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2,inten_sum2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

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
    
    % K9me1
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

% K9me2
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K9me3
hno = 5;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);%#ok

[His.rt_ref(5),His.rt_ref(4)] = find_pair_new(top1_rt5,rts4,top1_rt4,inten_sum4,0);

% K14ac
hno = 7;
t1 = (His.rt_ref(1)+His.rt_ref(4))/2;
t2 = His.rt_ref(1)-delta;
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% K9ac
hno = 6;
if 0==His.rt_ref(hno+1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno+1);
    d = His.rt_ref(hno+1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K9me1K14ac
hno = 8;
t1 = His.rt_ref(1)+delta;
if 0==His.rt_ref(2)
    t2 = His.rt_ref(1)+9;
else
    t2 = His.rt_ref(2)-delta;
end;
[rts8,top1_rt8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt8;
end;

% K9me2K14ac
hno = 9;
t1 = 6;
if 0==His.rt_ref(4)
    t2 = His.rt_ref(1)-4;
else
    t2 = His.rt_ref(4)-delta;
end;
[rts9,top1_rt9,inten_sum9] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K9me3K14ac
hno = 10;
t1 = 6;
if 0==His.rt_ref(4)
    t2 = His.rt_ref(1)-4;
else
    t2 = His.rt_ref(4)-delta;
end;
[rts10,top1_rt10,inten_sum10] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 0==isempty(top1_rt10) && top1_rt10>0 && His.rt_ref(4)>0 && top1_rt10<His.rt_ref(4)
    [His.rt_ref(10),His.rt_ref(9)] = find_pair_new(top1_rt10,rts9,top1_rt9,inten_sum9,0);
else
    [His.rt_ref(9),His.rt_ref(10)] = find_pair_new(top1_rt9,rts10,top1_rt10,inten_sum10,1);
end;

% K9acK14ac
hno = 11;
if 0==His.rt_ref(4)
    t1 = His.rt_ref(1)-18;
else
    t1 = His.rt_ref(4)+delta;
end;
if 0==His.rt_ref(6)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(6)-delta;
end;
[rts11,top1_rt11] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts11)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt11;
end;