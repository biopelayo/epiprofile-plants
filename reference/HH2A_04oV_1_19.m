function HH2A_04oV_1_19(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'HH2A_04oV_1_19';
fprintf(1,'%s..',out_filename(2:end));
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

His.pep_seq = 'AGGKAGKDSGKAKAKAVSR';
His.mod_short = {'unmod';
    'K4ac';
    'K7ac';
    'K11ac';
    'K15ac';
    'K4acK7ac';
    'K4acK11ac';
    'K4acK15ac';
    'K7acK11ac';
    'K7acK15ac';
    'K11acK15ac';
    'K7acK11acK15ac';
    'K4acK11acK15ac';
    'K4acK7acK15ac';
    'K4acK7acK11ac';
    'K4acK7acK11acK15ac'};
His.mod_type = {'0,pr;4,pr;7,pr;11,pr;13,pr;15,pr;';
    '0,pr;4,ac;7,pr;11,pr;13,pr;15,pr;';
    '0,pr;4,pr;7,ac;11,pr;13,pr;15,pr;';
    '0,pr;4,pr;7,pr;11,ac;13,pr;15,pr;';
    '0,pr;4,pr;7,pr;11,pr;13,pr;15,ac;';
    '0,pr;4,ac;7,ac;11,pr;13,pr;15,pr;';
    '0,pr;4,ac;7,pr;11,ac;13,pr;15,pr;';
    '0,pr;4,ac;7,pr;11,pr;13,pr;15,ac;';
    '0,pr;4,pr;7,ac;11,ac;13,pr;15,pr;';
    '0,pr;4,pr;7,ac;11,pr;13,pr;15,ac;';
    '0,pr;4,pr;7,pr;11,ac;13,pr;15,ac;';
    '0,pr;4,pr;7,ac;11,ac;13,pr;15,ac;';
    '0,pr;4,ac;7,pr;11,ac;13,pr;15,ac;';
    '0,pr;4,ac;7,ac;11,pr;13,pr;15,ac;';
    '0,pr;4,ac;7,ac;11,ac;13,pr;15,pr;';
    '0,pr;4,ac;7,ac;11,ac;13,pr;15,ac;'};

His.pep_ch = repmat([2 3 4 5],length(His.mod_type),1);
%{
His.pep_mz = [1062.0946	708.3988	531.5509	425.4422
    1055.0868	703.7269	528.0470	422.6391
    1055.0868	703.7269	528.0470	422.6391
    1055.0868	703.7269	528.0470	422.6391
    1055.0868	703.7269	528.0470	422.6391
    1048.0789	699.0550	524.5431	419.8359
    1048.0789	699.0550	524.5431	419.8359
    1048.0789	699.0550	524.5431	419.8359
    1048.0789	699.0550	524.5431	419.8359
    1048.0789	699.0550	524.5431	419.8359
    1048.0789	699.0550	524.5431	419.8359
    1041.0711	694.3832	521.0392	417.0328
    1041.0711	694.3832	521.0392	417.0328
    1041.0711	694.3832	521.0392	417.0328
    1041.0711	694.3832	521.0392	417.0328
    1034.0633	689.7113	517.5353	414.2297];
%}
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [37.41
    36.57
    36.57
    36.57
    36.57
    35.70
    35.70
    35.70
    35.70
    35.70
    35.70
    34.86
    34.86
    34.86
    34.86
    34.06];
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

% K4ac/K7ac/K11ac/K15ac
%{
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
%}
%
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone4(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+3,1:ncharge) = cur_rts(1:4,:);
    pep_intens(hno:hno+3,1:ncharge) = cur_intens(1:4,:);
    mono_isointens(1:num_MS1,hno:hno+3) = cur_mono_isointens(:,1:4);
end;
%}

% K4acK7ac/K4acK11ac/K4acK15ac/K7acK11ac/K7acK15ac/K11acK15ac
%{
hno = 6;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
%}
%
hno = 6;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone6(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+5,1:ncharge) = cur_rts(1:6,:);
    pep_intens(hno:hno+5,1:ncharge) = cur_intens(1:6,:);
    mono_isointens(1:num_MS1,hno:hno+5) = cur_mono_isointens(:,1:6);
end;
%}

% K7acK11acK15ac/K4acK11acK15ac/K4acK7acK15ac/K4acK7acK11ac
%{
hno = 12;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
%}
%
hno = 12;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone4(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+3,1:ncharge) = cur_rts(1:4,:);
    pep_intens(hno:hno+3,1:ncharge) = cur_intens(1:4,:);
    mono_isointens(1:num_MS1,hno:hno+3) = cur_mono_isointens(:,1:4);
end;
%}

% K4acK7acK11acK15ac
hno = 16;
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

% K4ac
hno = 2;
t1 = His.rt_ref(1)-8;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K7ac/K11ac/K15ac
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:5) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:5) = His.rt_ref(hno-1);
end;

% K4acK7ac
hno = 6;
t1 = His.rt_ref(1)-13;
if 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K4acK11ac/K4acK15ac/K7acK11ac/K7acK15ac/K11acK15ac
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:11) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:11) = His.rt_ref(hno-1);
end;

% K7acK11acK15ac
hno = 12;
t1 = His.rt_ref(1)-18;
if 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts12,top1_rt12] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts12)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt12;
end;

% K4acK11acK15ac/K4acK7acK15ac/K4acK7acK11ac
hno = 13;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:15) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:15) = His.rt_ref(hno-1);
end;

% K4acK7acK11acK15ac
hno = 16;
t1 = His.rt_ref(1)-22;
if 0<His.rt_ref(12)
    t2 = His.rt_ref(12)-delta;
elseif 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts16,top1_rt16] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts16)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt16;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K4ac
hno = 2;
t1 = His.rt_ref(1)-8;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K7ac/K11ac/K15ac
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:5) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:5) = His.rt_ref(hno-1);
end;

% K4acK7ac
hno = 6;
t1 = His.rt_ref(1)-13;
if 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end;

% K4acK11ac/K4acK15ac/K7acK11ac/K7acK15ac/K11acK15ac
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:11) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:11) = His.rt_ref(hno-1);
end;

% K7acK11acK15ac
hno = 12;
t1 = His.rt_ref(1)-18;
if 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts12,top1_rt12] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts12)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt12;
end;

% K4acK11acK15ac/K4acK7acK15ac/K4acK7acK11ac
hno = 13;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:15) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:15) = His.rt_ref(hno-1);
end;

% K4acK7acK11acK15ac
hno = 16;
t1 = His.rt_ref(1)-22;
if 0<His.rt_ref(12)
    t2 = His.rt_ref(12)-delta;
elseif 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end;
[rts16,top1_rt16] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts16)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt16;
end;