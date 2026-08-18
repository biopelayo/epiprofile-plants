function H3_04v3a_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'H3_04v3a_27_40';
% fprintf(1,'%s..',out_filename);
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

His.pep_seq = 'KSAPSTGGVKKPHR';
His.mod_short = {'S28ph';
    'K27me1S28ph';
    'K27me2S28ph';
    'K27me3S28ph';
    'K27me2S28phK36me2';
    'K27me3S28phK36me2';
    'S28ac';
    'S28pr';
    'K27me2S28ac';
    'K27me3S28ac';
    'K27me3S28acK36me1'};
His.mod_type = {'0,pr;1,pr;2,ph;10,pr;11,pr;';
    '0,pr;1,me1;2,ph;10,pr;11,pr;';
    '0,pr;1,me2;2,ph;10,pr;11,pr;';
    '0,pr;1,me3;2,ph;10,pr;11,pr;';
    '0,pr;1,me2;2,ph;10,me2;11,pr;';
    '0,pr;1,me3;2,ph;10,me2;11,pr;';
    '0,pr;1,pr;2,ac;10,pr;11,pr;';
    '0,pr;1,pr;2,pr;10,pr;11,pr;';
    '0,pr;1,me2;2,ac;10,pr;11,pr;';
    '0,pr;1,me3;2,ac;10,pr;11,pr;';
    '0,pr;1,me3;2,ac;10,me1;11,pr;'};

His.pep_ch = repmat([2 3 4],length(His.mod_type),1);
%{
His.pep_mz = [877.4535	585.3047	439.2304
    884.4613	589.9766	442.7343
    863.4560	575.9731	432.2316
    870.4638	580.6450	435.7356
    
    
    858.4756	572.6528	429.7414
    865.4834	577.3247	433.2453
    844.4781	563.3212	422.7427
    851.4860	567.9931	426.2466
    858.4938	572.6649	429.7505];
%}
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [31.35
    33.20
    27.34
    27.31
    25
    25
    31.35
    33.20
    27.34
    27.31
    29.15];
His.display = ones(length(His.mod_type),1);
His.display(8) = 0;

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

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special)%#ok
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        nhmass = special.nhmass;
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end;
end;

% S28ph
% K27me1S28ph
% K27me2S28ph
% K27me3S28ph
% K27me2S28phK36me2
% K27me3S28phK36me2
% S28ac
% S28pr
% K27me2S28ac
% K27me3S28ac
% K27me3S28acK36me1
for hno=1:11
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone12(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.3;
nsplit = 1;

ref_rts = get_KSAPSTGGVKKPHR_rt(His);
if 1==isempty(ref_rts)
    return;
end;

% S28ph
hno = 1;
t1 = ref_rts(1)+delta;
t2 = ref_rts(1)+20;
[rts1,top1_rt1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts1)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt1;
end;

% K27me1S28ph
hno = 2;
if His.rt_ref(1)>0
    t1 = His.rt_ref(1)+delta;
elseif ref_rts(3)>0
    t1 = ref_rts(3)+delta;
else
    t1 = ref_rts(1)+delta;
end;
if His.rt_ref(1)>0
    t2 = His.rt_ref(1)+11;
elseif ref_rts(3)>0
    t2 = ref_rts(3)+18;
else
    t2 = ref_rts(1)+25;
end;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K27me2S28ph
hno = 3;
if 0==ref_rts(4)
    t1 = ref_rts(1)-35;
else
    t1 = ref_rts(4)+delta;
end;
if 0==ref_rts(9)
    t2 = ref_rts(1)-delta;
else
    t2 = ref_rts(9)-delta;
end;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K27me3S28ph
hno = 4;
if 0==ref_rts(4)
    t1 = ref_rts(1)-35;
else
    t1 = ref_rts(4)+delta;
end;
if 0==ref_rts(9)
    t2 = ref_rts(1)-delta;
else
    t2 = ref_rts(9)-delta;
end;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[His.rt_ref(3),His.rt_ref(4)] = find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3,rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K27me2S28phK36me2
hno = 5;
if 0==His.rt_ref(3)
    t1 = ref_rts(1)-45;
else
    t1 = His.rt_ref(3)-10;
end;
if 0==His.rt_ref(3)
    t2 = ref_rts(1)-10;
else
    t2 = His.rt_ref(3)-delta;
end;
[rts5,top1_rt5,inten_sum5,top1_inten_sum5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K27me3S28phK36me2
hno = 6;
if 0==His.rt_ref(3)
    t1 = ref_rts(1)-45;
else
    t1 = His.rt_ref(3)-10;
end;
if 0==His.rt_ref(3)
    t2 = ref_rts(1)-10;
else
    t2 = His.rt_ref(3)-delta;
end;
[rts6,top1_rt6,inten_sum6,top1_inten_sum6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[His.rt_ref(5),His.rt_ref(6)] = find_pair(rts5,top1_rt5,inten_sum5,top1_inten_sum5,rts6,top1_rt6,inten_sum6,top1_inten_sum6);

% S28ac
hno = 7;
t1 = ref_rts(1)+delta;
if ref_rts(3)>0
    t2 = ref_rts(3)+11;
else
    t2 = ref_rts(1)+18;
end;
[rts7,top1_rt7] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% S28pr
hno = 8;
if His.rt_ref(7)>0
    t1 = His.rt_ref(7)+delta;
elseif ref_rts(3)>0
    t1 = ref_rts(3)+delta;
else
    t1 = ref_rts(1)+delta;
end;
if His.rt_ref(7)>0
    t2 = His.rt_ref(7)+11;
elseif ref_rts(3)>0
    t2 = ref_rts(3)+18;
else
    t2 = ref_rts(1)+25;
end;
[rts8,top1_rt8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt8;
end;

% K27me2S28ac
hno = 9;
if ref_rts(4)>0
    t1 = ref_rts(4)+delta;
else
    t1 = delta;
end;
t2 = ref_rts(1)-delta;
[rts9,top1_rt9,inten_sum9,top1_inten_sum9] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K27me3S28ac
hno = 10;
if ref_rts(4)>0
    t1 = ref_rts(4)+delta;
else
    t1 = delta;
end;
t2 = ref_rts(1)-delta;
[rts10,top1_rt10,inten_sum10,top1_inten_sum10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[His.rt_ref(9),His.rt_ref(10)] = find_pair(rts9,top1_rt9,inten_sum9,top1_inten_sum9,rts10,top1_rt10,inten_sum10,top1_inten_sum10);

% K27me3S28acK36me1
hno = 11;
if His.rt_ref(10)>0
    t1 = His.rt_ref(10)+delta;
else
    t1 = delta;
end;
t2 = ref_rts(1)+delta;
[rts11,top1_rt11] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts11)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt11;
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.3;
nsplit = 1;

ref_rts = get_KSAPSTGGVKKPHR_rt(His);
if 1==isempty(ref_rts)
    return;
end;

% S28ph
hno = 1;
t1 = ref_rts(1)+delta;
t2 = ref_rts(1)+20;
[rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts1)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt1;
end;

% K27me1S28ph
hno = 2;
if His.rt_ref(1)>0
    t1 = His.rt_ref(1)+delta;
elseif ref_rts(3)>0
    t1 = ref_rts(3)+delta;
else
    t1 = ref_rts(1)+delta;
end;
if His.rt_ref(1)>0
    t2 = His.rt_ref(1)+11;
elseif ref_rts(3)>0
    t2 = ref_rts(3)+18;
else
    t2 = ref_rts(1)+25;
end;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K27me2S28ph
hno = 3;
if 0==ref_rts(4)
    t1 = ref_rts(1)-35;
else
    t1 = ref_rts(4)+delta;
end;
if 0==ref_rts(9)
    t2 = ref_rts(1)-delta;
else
    t2 = ref_rts(9)-delta;
end;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K27me3S28ph
hno = 4;
if 0==ref_rts(4)
    t1 = ref_rts(1)-35;
else
    t1 = ref_rts(4)+delta;
end;
if 0==ref_rts(9)
    t2 = ref_rts(1)-delta;
else
    t2 = ref_rts(9)-delta;
end;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[His.rt_ref(3),His.rt_ref(4)] = find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3,rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K27me2S28phK36me2
hno = 5;
if 0==His.rt_ref(3)
    t1 = ref_rts(1)-45;
else
    t1 = His.rt_ref(3)-10;
end;
if 0==His.rt_ref(3)
    t2 = ref_rts(1)-10;
else
    t2 = His.rt_ref(3)-delta;
end;
[rts5,top1_rt5,inten_sum5,top1_inten_sum5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K27me3S28phK36me2
hno = 6;
if 0==His.rt_ref(3)
    t1 = ref_rts(1)-45;
else
    t1 = His.rt_ref(3)-10;
end;
if 0==His.rt_ref(3)
    t2 = ref_rts(1)-10;
else
    t2 = His.rt_ref(3)-delta;
end;
[rts6,top1_rt6,inten_sum6,top1_inten_sum6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[His.rt_ref(5),His.rt_ref(6)] = find_pair(rts5,top1_rt5,inten_sum5,top1_inten_sum5,rts6,top1_rt6,inten_sum6,top1_inten_sum6);

% S28ac
hno = 7;
t1 = ref_rts(1)+delta;
if ref_rts(3)>0
    t2 = ref_rts(3)+11;
else
    t2 = ref_rts(1)+18;
end;
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts7)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt7;
end;

% S28pr
hno = 8;
if His.rt_ref(7)>0
    t1 = His.rt_ref(7)+delta;
elseif ref_rts(3)>0
    t1 = ref_rts(3)+delta;
else
    t1 = ref_rts(1)+delta;
end;
if His.rt_ref(7)>0
    t2 = His.rt_ref(7)+11;
elseif ref_rts(3)>0
    t2 = ref_rts(3)+18;
else
    t2 = ref_rts(1)+25;
end;
[rts8,top1_rt8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts8)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt8;
end;

% K27me2S28ac
hno = 9;
if ref_rts(4)>0
    t1 = ref_rts(4)+delta;
else
    t1 = delta;
end;
t2 = ref_rts(1)-delta;
[rts9,top1_rt9,inten_sum9,top1_inten_sum9] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K27me3S28ac
hno = 10;
if ref_rts(4)>0
    t1 = ref_rts(4)+delta;
else
    t1 = delta;
end;
t2 = ref_rts(1)-delta;
[rts10,top1_rt10,inten_sum10,top1_inten_sum10] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[His.rt_ref(9),His.rt_ref(10)] = find_pair(rts9,top1_rt9,inten_sum9,top1_inten_sum9,rts10,top1_rt10,inten_sum10,top1_inten_sum10);

% K27me3S28acK36me1
hno = 11;
if His.rt_ref(10)>0
    t1 = His.rt_ref(10)+delta;
else
    t1 = delta;
end;
t2 = ref_rts(1)+delta;
[rts11,top1_rt11] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts11)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt11;
end;

function ref_rts = get_KSAPSTGGVKKPHR_rt(His)
%%

out_file0 = fullfile(His.outpath,'H3_04v3_27_40.xls');
if 0~=exist(out_file0,'file')
    fp = fopen(out_file0,'r');
    str = fgetl(fp);
    while 0==feof(fp) && 0==strcmp(str,'[rt]')
        str = fgetl(fp);
    end;
    str = fgetl(fp);%#ok peptide
    no = 0;
    ref_rts = [];
    while 0==feof(fp)
        str = fgetl(fp);
        p = strfind(str,'	');
        c_rt = str2num( str(p(1)+1:p(2)-1) );%#ok
        no = no + 1;
        ref_rts(no) = c_rt;%#ok
    end;
    fclose(fp);
else
    ref_rts = [];
end;