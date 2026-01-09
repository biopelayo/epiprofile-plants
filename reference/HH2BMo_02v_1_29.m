function HH2BMo_02v_1_29(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'HH2B_02v_1_29';
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
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

function His = init_histone()
%%

His.pep_seq = 'unmod';
His.mod_short = {'1C.PEPAKSAPAPKKGSKKAVTKAQKKDGKKR';
    '1H.PEPAKSAPAPKKGSKKALTKAQKKDGKKR';
    '2B.PDPAKSAPAPKKGSKKAVTKVQKKDGKKR';
    '1B.PEPSKSAPAPKKGSKKAISKAQKKDGKKR';
    '3B.PDPSKSAPAPKKGSKKAVTKAQKKDGKKR';
    '1M.PEPTKSAPAPKKGSKKAVTKAQKKDGKKR';
    '1P.PEPVKSVPAPKKGSKKAVTKAQKKDGKKR';
    '2E.PELAKSAPAPKKGSKKAVTKAQKKDGKKR'};
His.mod_short2 = {'1C/1K/1J';
    '1H';
    '2B';
    '1B';
    '3B';
    '1M';
    '1P';
    '2E'};
His.mod_type = {'0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;';
    '0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;';
    '0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;';
    '0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;';
    '0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;';
    '0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;';
    '0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;';
    '0,pr;5,pr;11,pr;12,pr;15,pr;16,pr;20,pr;23,pr;24,pr;27,pr;28,pr;'};

His.pep_ch = repmat([4 3],length(His.mod_type),1);
%{
His.pep_mz = [919.7816	1226.0397
    923.2855	1230.7116
    923.2855	1230.7116
    923.7804	1231.3714
    920.2764	1226.6995
    927.2843	1236.0433
    933.7973	1244.7273
    923.7895	1231.3835];
%}
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [37.89
    38.98
    39.41
    37.35
    37.89
    37.76
    38.89
    38.28];
His.display = ones(length(His.mod_type),1);

% main ch
main_ch = His.pep_ch(1,1);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);%#ok
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = [1 3];%1:npep;
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

% '1C/1K/1J'
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        p = strfind(His.mod_short{1},'.');
        [His.rt_ref(1),special.ndebug] = check_ref(special.raw_path,[His.mod_short{1}(p+1:end),His.mod_type{1}],His.rt_ref(1),special.ndebug);
    else
        nhmass = special.nhmass;
        p = strfind(His.mod_short{1},'.');
        His.rt_ref(1) = check_ref(special.raw_path,[His.mod_short{1}(p+1:end),His.mod_type{1}],His.rt_ref(1),special.ndebug);
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

% '1H'
% '2B'
% '1B'
% '3B'
% '1M'
% '1P'
% '2E'
for hno=2:3
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;
for hno=4:8
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

% '1H'
hno = 2;
t1 = His.rt_ref(1)-delta;
t2 = His.rt_ref(1)+5;
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
    
    % '2B'
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

% '1B'
hno = 4;
t1 = His.rt_ref(1)-5;
t2 = His.rt_ref(1)+delta;
[rts4,top1_rt4,inten_sum4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);%#ok
if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    id = find(rts4<His.rt_ref(1));
    if 0==isempty(id)
        [tmp,ix] = max(inten_sum4(id));%#ok
        His.rt_ref(hno) = rts4(id(ix));
    else
        His.rt_ref(hno) = 0;
    end;
end;

% '3B'
hno = 5;
t1 = His.rt_ref(1)-5;
t2 = His.rt_ref(1)+5;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% '1M'
hno = 6;
t1 = His.rt_ref(1)-5;
t2 = His.rt_ref(1)+5;
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% '1P'
hno = 7;
t1 = His.rt_ref(1)-5;
t2 = His.rt_ref(1)+5;
[rts7,top1_rt7] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% '2E'
hno = 8;
t1 = His.rt_ref(1)-delta;
t2 = His.rt_ref(1)+5;
[rts8,top1_rt8,inten_sum8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);%#ok
if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    id = find(rts8>His.rt_ref(1));
    if 0==isempty(id)
        [tmp,ix] = max(inten_sum8(id));%#ok
        His.rt_ref(hno) = rts8(id(ix));
    else
        His.rt_ref(hno) = 0;
    end;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% '1H'
hno = 2;
t1 = His.rt_ref(1)-delta;
t2 = His.rt_ref(1)+5;
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
    
    % '2B'
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;

% '1B'
hno = 4;
t1 = His.rt_ref(1)-5;
t2 = His.rt_ref(1)+delta;
[rts4,top1_rt4,inten_sum4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);%#ok
if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    id = find(rts4<His.rt_ref(1));
    if 0==isempty(id)
        [tmp,ix] = max(inten_sum4(id));%#ok
        His.rt_ref(hno) = rts4(id(ix));
    else
        His.rt_ref(hno) = 0;
    end;
end;

% '3B'
hno = 5;
t1 = His.rt_ref(1)-5;
t2 = His.rt_ref(1)+5;
[rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

% '1M'
hno = 6;
t1 = His.rt_ref(1)-5;
t2 = His.rt_ref(1)+5;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% '1P'
hno = 7;
t1 = His.rt_ref(1)-5;
t2 = His.rt_ref(1)+5;
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% '2E'
hno = 8;
t1 = His.rt_ref(1)-delta;
t2 = His.rt_ref(1)+5;
[rts8,top1_rt8,inten_sum8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);%#ok
if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    id = find(rts8>His.rt_ref(1));
    if 0==isempty(id)
        [tmp,ix] = max(inten_sum8(id));%#ok
        His.rt_ref(hno) = rts8(id(ix));
    else
        His.rt_ref(hno) = 0;
    end;
end;