function [rts,top1_rt,inten_sum,top1_inten_sum] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2)
%%

num_MS1 = size(MS1_index,1);
if t1>MS1_index(num_MS1,2) || t2<0 || t2<t1
    rts = [];
    top1_rt = [];
    inten_sum = [];
    top1_inten_sum = [];
    return;
end;

% get MS1 profile
c_mz = His.pep_mz(hno,1);
c_ch = His.pep_ch(hno,1);
% 2º error del codigo lo mismopara rts
% p = find( MS1_index(:,2)>=t1 );
% rt_i1 = p(1);
% pp = find( MS1_index(:,2)<=t2 );
% rt_i2 = pp(end);
rtcol = MS1_index(:,2);

p  = find(rtcol >= t1);
pp = find(rtcol <= t2);

if isempty(p) || isempty(pp)
    rts = [];
    top1_rt = [];
    inten_sum = [];
    top1_inten_sum = [];
    return;
end

rt_i1 = p(1);
rt_i2 = pp(end);
if isempty(p) || isempty(pp)
    fprintf('get_rts window empty: hno=%d t1=%g t2=%g  RT=[%g,%g]\n', ...
        hno, t1, t2, min(rtcol), max(rtcol));
    return;
end

if rt_i2 < rt_i1
    rts = [];
    top1_rt = [];
    inten_sum = [];
    top1_inten_sum = [];
    return;
end

c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
if ptol==100
    ptol = 10;
end;
if ptol>100 && c_ch>=3
    nC13 = 1;
else
    nC13 = 0;
end;
[c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
c_mono_isointens = c_ref_isointens(:,2);
if 1==nsplit
    [nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens);%#ok
else
    [nt,nb,top1_idx,inten_sum] = GetTopBottom11(c_mono_isointens);%#ok
end;
rts = c_isorts(nt);
top1_rt = c_isorts(nt(top1_idx));
top1_inten_sum = inten_sum(top1_idx);

if top1_rt<4
    rts = [];
    top1_rt = [];
    inten_sum = [];
    top1_inten_sum = [];
    return;
end;