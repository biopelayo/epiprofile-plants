function HH2AN_02m1_4_11(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'HH2A_02m1_4_11';
fprintf(1,'%s..',out_filename(2:end));
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% init
His = init_histone(cur_outpath,out_filename,special);

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

His.pep_seq = 'GKQGGKAR';
His.mod_short = {'unmod';
    'K5ac';
    'K9ac';
    'K5acK9ac';
    'K9me1';
    'K5me1'};
if 1==special.nsubtype
    His.mod_type = {'0,pr;2,pr;6,pr;';
        '0,pr;2,ac;6,pr;';
        '0,pr;2,pr;6,ac;';
        '0,pr;2,ac;6,ac;';
        '0,pr;2,pr;6,me1;';
        '0,pr;2,me1;6,pr;'};
else
    His.mod_type = {'0,pr;2,pr;6,pr;';
        '0,pr;2,hac;6,pr;';
        '0,pr;2,pr;6,hac;';
        '0,pr;2,hac;6,hac;';
        '0,pr;2,pr;6,hme1;';
        '0,pr;2,hme1;6,pr;'};
end;

His.pep_ch = repmat([1 2],length(His.mod_type),1);
%{
His.pep_mz = [969.5476	485.2774
    955.5320	478.2696
    955.5320	478.2696
    941.5163	471.2618
    983.5633	492.2853
    983.5633	492.2853];
%}
if 2==special.nsubtype
    His.pep_mz = calculate_pepmz(His);
else
    His.pep_mz = calculate_pepmzH(His);
end;
His.rt_ref = [24.23
    22.54
    22.67
    20.8
    26.33
    26.41];
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

% K5ac/K9ac
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;

% K5acK9ac
hno = 4;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K9me1/K5me1
if His.rt_ref(6)-His.rt_ref(5)>0.4
    for hno=5:6
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone10(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end;
    end;
else
    hno = 5;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end;
end;

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1;
nsplit = 1;

% K5ac
hno = 2;
t1 = His.rt_ref(1)-14;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K9ac
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K5acK9ac
hno = 4;
t1 = His.rt_ref(1)-20;
if 0==His.rt_ref(2)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(2)-delta;
end;
[rts4,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt4;
end;

% K9me1
hno = 5;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts5,top1_rt5,inten_sum5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[tmp_sum,ix] = sort(inten_sum5,'descend');
tmp_rts = rts5(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
    His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts5)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt5;
    end;
    
    % K5me1
    hno = 6;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K9ac
hno = 3;
t1 = His.rt_ref(1)-14;
t2 = His.rt_ref(1)-delta;
[rts3,top1_rt3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end;

% K5ac
hno = 2;
if 0==His.rt_ref(hno+1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno+1);
    d = His.rt_ref(hno+1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K5acK9ac
hno = 4;
t1 = His.rt_ref(1)-20;
if 0==His.rt_ref(2)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(2)-delta;
end;
[rts4,top1_rt4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt4;
end;

% K9me1
hno = 5;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts5,top1_rt5,inten_sum5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,0,t1,t2,nhmass);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% K5me1
hno = 6;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,0,t1,t2,nhmass);

if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

if 0==isempty(rts5) && 0==isempty(rts6)
    if top1_rt5==top1_rt6
        hno = 5;
        [tmp_sum,ix] = sort(inten_sum5,'descend');
        tmp_rts = rts5(ix);
        if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
            His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
            His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
        else
            old_t = His.rt_ref(hno);
            if 1==isempty(rts5)
                His.rt_ref(hno) = 0;
            else
                His.rt_ref(hno) = top1_rt5;
            end;
            
            % K5me1
            hno = 6;
            if 0==His.rt_ref(hno-1)
                His.rt_ref(hno) = 0;
            elseif old_t~= His.rt_ref(hno-1);
                d = His.rt_ref(hno-1) - old_t;
                His.rt_ref(hno) = His.rt_ref(hno) + d;
            end;
        end;
    elseif top1_rt5>top1_rt6 && abs(top1_rt5-top1_rt6)<2
        hno = 5;
        His.rt_ref(hno) = top1_rt6;
        His.rt_ref(hno+1) = top1_rt5;
    end;
end;