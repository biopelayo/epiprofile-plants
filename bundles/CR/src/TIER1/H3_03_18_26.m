function H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_03_18_26';
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

His.pep_seq = 'KQLATKAAR';
His.mod_short = {'unmod';
    'K23me1';
    'K18me1';
    'K18me1K23me1';
    'K18ac';
    'K23ac';
    'K18acK23ac';
    'K18acT22pr';
    'T22ac';
    'T22pr'};
His.mod_type = {'0,pr;1,pr;6,pr;';
    '0,pr;1,pr;6,me1;';
    '0,pr;1,me1;6,pr;';
    '0,pr;1,me1;6,me1;';
    '0,pr;1,ac;6,pr;';
    '0,pr;1,pr;6,ac;';
    '0,pr;1,ac;6,ac;';
    '0,pr;1,ac;5,pr;6,pr;';
    '0,pr;1,pr;5,ac;6,pr;';
    '0,pr;1,pr;5,pr;6,pr;'};

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
%{
His.pep_mz = [1154.6892	577.8482	385.5679
    1168.7048	584.8561	390.2398
    1168.7048	584.8561	390.2398
    1182.7205	591.8639	394.9117
    1140.6735	570.8404	380.8960
    1140.6735	570.8404	380.8960
    1126.6579	563.8326	376.2241
    1196.6997	598.8535	399.5714
    1196.6997	598.8535	399.5714
    1210.7154	605.8613	404.2433];
%}
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [37.41
    39.76
    39.77
    40.96
    36.13
    36.14
    34.64
    40.04
    40.05
    41.27];
His.display = ones(length(His.mod_type),1);
His.display([8 9 10]) = 0;

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
            hno = 2;% K23me1
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2);%#ok
            hno = 5;% K18ac
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);%#ok
            if 0==isempty(top1_rt2) && 0==isempty(top1_rt5) && top1_rt2>top1_rt5
                if 1==isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt5+(top1_rt2-top1_rt5)*0.4;
                else
                    p = find(rts1>top1_rt5+(top1_rt2-top1_rt5)*0.2 & rts1<top1_rt2);
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

% K23me1/K18me1
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

% K18me1K23me1
hno = 4;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K18ac/K23ac
hno = 5;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;

% K18acK23ac
hno = 7;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% T22ac/K18acT22pr
if His.rt_ref(9)-His.rt_ref(8)>0.4
    for hno=8:9
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone10(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end;
    end;
else
    hno = 8;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end;
end;

% T22pr
hno = 10;
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

% K23me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
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
    
    % K18me1
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

% K18me1K23me1
hno = 4;
if 0==His.rt_ref(3)
    t1 = His.rt_ref(1)+delta;
else
    t1 = His.rt_ref(3)+delta;
end;
if 0==His.rt_ref(3)
    t2 = His.rt_ref(1)+22;
else
    t2 = His.rt_ref(3)+(His.rt_ref(3)-His.rt_ref(1))+2;
end;
[rts4,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt4;
end;

% K18ac
hno = 5;
t1 = His.rt_ref(1)-14;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% K23ac
hno = 6;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K18acK23ac
hno = 7;
t1 = His.rt_ref(1)-20;
if 0==His.rt_ref(5)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(5)-delta;
end;
[rts7,top1_rt7] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% T22ac
hno = 8;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts8,top1_rt8,inten_sum8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[tmp_sum,ix] = sort(inten_sum8,'descend');
tmp_rts = rts8(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
    His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts8)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt8;
    end;
    
    % K18acT22pr
    hno = 9;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

% T22pr
hno = 10;
if 0==His.rt_ref(8)
    t1 = His.rt_ref(1)+delta;
else
    t1 = His.rt_ref(8)+delta;
end;
t2 = His.rt_ref(1)+20;
[rts10,top1_rt10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts10)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt10;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K23me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts2,top1_rt2,inten_sum2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,0,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K18me1
hno = 3;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
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
            
            % K18me1
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

% K18me1K23me1
hno = 4;
if 0==His.rt_ref(3)
    t1 = His.rt_ref(1)+delta;
else
    t1 = His.rt_ref(3)+delta;
end;
if 0==His.rt_ref(3)
    t2 = His.rt_ref(1)+22;
else
    t2 = His.rt_ref(3)+(His.rt_ref(3)-His.rt_ref(1))+2;
end;
[rts4,top1_rt4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt4;
end;

% K23ac
hno = 6;
t1 = His.rt_ref(1)-14;
t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K18ac
hno = 5;
if 0==His.rt_ref(hno+1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno+1);
    d = His.rt_ref(hno+1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K18acK23ac
hno = 7;
t1 = His.rt_ref(1)-20;
if 0==His.rt_ref(5)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(5)-delta;
end;
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% T22ac
hno = 8;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts8,top1_rt8,inten_sum8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[tmp_sum,ix] = sort(inten_sum8,'descend');
tmp_rts = rts8(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
    His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts8)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt8;
    end;
    
    % K18acT22pr
    hno = 9;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

% T22pr
hno = 10;
if 0==His.rt_ref(8)
    t1 = His.rt_ref(1)+delta;
else
    t1 = His.rt_ref(8)+delta;
end;
t2 = His.rt_ref(1)+20;
[rts10,top1_rt10] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts10)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt10;
end;