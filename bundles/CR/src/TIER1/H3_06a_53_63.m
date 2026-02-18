function H3_06a_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_06a_53_63';
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

His.pep_seq = 'RYQKSTELLIR';
His.mod_short = {'unmod';
    'S57pr';
    'Y54pr'};
His.mod_type = {'0,pr;4,pr;';
    '0,pr;4,pr;5,pr;';
    '0,pr;2,pr;4,pr;'};

His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);
%{
His.pep_mz = [1518.8639	759.9356	506.9595	380.4714
    1574.8901	787.9487	525.6349	394.4780
    1574.8901	787.9487	525.6349	394.4780];
%}
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [37
    41.58
    42];
His.display = zeros(length(His.mod_type),1);

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

% S57pr/Y54pr
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

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1;
nsplit = 0;

% S57pr
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts2,top1_rt2,inten_sum2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[tmp_sum,ix] = sort(inten_sum2,'descend');
tmp_rts = rts2(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<5
    His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end;
    
    % Y54pr
    hno = 3;
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
nsplit = 0;

% S57pr
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts2,top1_rt2,inten_sum2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[tmp_sum,ix] = sort(inten_sum2,'descend');
tmp_rts = rts2(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<5
    His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end;
    
    % Y54pr
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1);
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end;
end;