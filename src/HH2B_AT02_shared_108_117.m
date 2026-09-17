function HH2B_AT02_shared_108_117(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'HH2B_AT02_shared_108_117';
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

His.pep_seq = 'YNKKPTITSR';
His.mod_short = {'unmod';
    'K3me1';
    'K3me2';
    'K3me3';
    'K3ac';
    'K4ac';
    'K3me1K4ac';
    'K3acK4ac'};
His.mod_type = {'0,pr;3,pr;4,pr;';
    '0,pr;3,me1;4,pr;';
    '0,pr;3,me2;4,pr;';
    '0,pr;3,me3;4,pr;';
    '0,pr;3,ac;4,pr;';
    '0,pr;3,pr;4,ac;';
    '0,pr;3,me1;4,ac;';
    '0,pr;3,ac;4,ac;'};

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [27.93
    30.60
    21.30
    21.20
    26.50
    26.60
    29.40
    24.90];
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
        if 2==special.ndebug
            hno = 1;% unmod
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2);
            hno = 2;% K3me1
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);%#ok
            [tmp_sum,ix] = sort(inten_sum1,'descend');
            tmp_rts = rts1(ix);
            if length(tmp_sum)>=3 && tmp_sum(3)>=tmp_sum(1)/30
                new_rts = sort(tmp_rts(1:3));
                if new_rts(2)-new_rts(1)>new_rts(3)-new_rts(2)
                    His.rt_ref(1) = new_rts(2);
                else
                    His.rt_ref(1) = new_rts(1);
                end;
            elseif length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<7
                His.rt_ref(1) = min([tmp_rts(2),tmp_rts(1)]);
            else
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

% K3me1
% K3me2
% K3me3
for hno=2:4
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

% K3ac/K4ac
hno = 5;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end;

% K3me1K4ac
hno = 7;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% K3acK4ac
hno = 8;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.5;
nsplit = 1;

% K3me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K3me2
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K3me3
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);%#ok

[His.rt_ref(4),His.rt_ref(3)] = find_pair_new(top1_rt4,rts3,top1_rt3,inten_sum3,0);

% K3ac
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

% K4ac
hno = 6;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1);
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K3me1K4ac
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

% K3acK4ac
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

delta = 0.5;
nsplit = 1;

% K3me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K3me2
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K3me3
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);%#ok

[His.rt_ref(4),His.rt_ref(3)] = find_pair_new(top1_rt4,rts3,top1_rt3,inten_sum3,0);

% K4ac
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

% K3ac
hno = 5;
if 0==His.rt_ref(hno+1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno+1);
    d = His.rt_ref(hno+1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end;

% K3me1K4ac
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

% K3acK4ac
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
[rts8,top1_rt8] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt8;
end;
