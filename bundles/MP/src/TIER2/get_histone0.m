function [cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special)
%%

[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([1,ncharge]);
cur_intens = zeros([1,ncharge]);
num_MS1 = size(MS1_index,1);
end_rt = MS1_index(num_MS1,2);
if end_rt>90
    gradient = 2;
else
    gradient = 1;
end;

rtcol = MS1_index(:,2);

% get MS1 profile
delta = 1;
c_mz = His.pep_mz(hno,1);
c_ch = His.pep_ch(hno,1);
c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
if ptol==100
    ptol = 10;
end;
if ptol>100 && c_ch>=3
    nC13 = 1;
else
    nC13 = 0;
end;

% --- Ventana inicial para el charge principal ---
if His.rt_unmod_orig==His.rt_ref(1) && 1~=special.ndebug
    rt_i1 = 1;
    rt_i2 = num_MS1;
else
    rt1 = His.rt_ref(hno)-5;
    rt2 = His.rt_ref(hno)+5;

    p  = find(rtcol >= rt1);
    pp = find(rtcol <= rt2);

    if isempty(p) || isempty(pp) || p(1)>pp(end)
        % No hay scans en ventana: devolver no-detectado sin crashear
        % fprintf('get_histone0: ventana vacía (init) hno=%d rt1=%g rt2=%g RT=[%g,%g]\n', hno, rt1, rt2, min(rtcol), max(rtcol));
        cur_mono_isointens = zeros(num_MS1,1);
        return;
    end

    rt_i1 = p(1);
    rt_i2 = pp(end);
end

[c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
j = 2;
c_mono_isointens = c_ref_isointens(:,j);
cur_mono_isointens = c_mono_isointens;

% get rt and area
[nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens);

if 1==special.ndebug || -1==special.ndebug
    if 1==special.ndebug
        llimit = -delta;
        rlimit = delta;
    else
        llimit = -2;
        if 1==isfield(His,'outfile') && 1==ismember(His.outfile,{'H3_02_9_17','H3_02a_9_17'})
            rlimit = 3;
        else
            rlimit = 2;
        end;
    end;
    ref_rt = His.rt_ref(hno);
    flag = zeros([1,length(nt)]);
    for i=1:length(nt)
        if llimit<=c_isorts(nt(i))-ref_rt && c_isorts(nt(i))-ref_rt<=rlimit
            flag(i) = 1;
        end;
    end;
    x = find(flag==1);
    if 0==isempty(x)
        if -1==special.ndebug && 1==isfield(His,'outfile') && 1==ismember(His.outfile,{'H3_02_9_17','H3_02a_9_17'}) && length(x)>=2
            [tmp,ix] = sort(inten_sum(x),'descend');%#ok
            ix = ix(1:2);
            [tmp,xx] = min( c_isorts(nt(x(ix))) );%#ok
            id = ix(xx);
        else
            [tmp,id] = max(inten_sum(x));%#ok
        end;
        top1_idx = x(id);
        cur_pos = nt(top1_idx);
        cur_rts(1) = c_isorts(cur_pos);
        cur_intens(1) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
    else
        cur_rts(1) = ref_rt;
        cur_intens(1) = 0;
    end;
else
    if His.rt_unmod_orig~=His.rt_ref(1) || ( 1==isfield(His,'outfile') && 1==strcmp(His.outfile,'H3_04v3_27_40') )
        if His.rt_unmod_orig~=His.rt_ref(1)
            ref_rt = His.rt_ref(1);
        elseif 1==isfield(His,'outfile') && 1==strcmp(His.outfile,'H3_04v3_27_40')
            ref_rt = get_H33_27_40_unmod_rt(His,gradient);
        end;
        flag = zeros([1,length(nt)]);
        for i=1:length(nt)
            if abs( c_isorts(nt(i))-ref_rt )<=delta
                flag(i) = 1;
            end;
        end;
        x = find(flag==1);
        if 0==isempty(x)
            [tmp,id] = max(inten_sum(x));%#ok
            top1_idx = x(id);
            cur_pos = nt(top1_idx);
            cur_rts(1) = c_isorts(cur_pos);
            cur_intens(1) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
        else
            cur_rts(1) = ref_rt;
            cur_intens(1) = 0;
        end;
    else
        if 0==isempty(nt) && c_isorts(nt(top1_idx))>4
            cur_pos = nt(top1_idx);
            cur_rts(1) = c_isorts(cur_pos);
            cur_intens(1) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
        else
            cur_rts(1) = His.rt_ref(hno);
            cur_intens(1) = 0;
        end;
    end;
end;

% --- Ventana para charges 2..ncharge basada en cur_rts(1) ---
rt1 = cur_rts(1)-delta;
rt2 = cur_rts(1)+delta;

p  = find(rtcol >= rt1);
pp = find(rtcol <= rt2);

if isempty(p) || isempty(pp) || p(1)>pp(end)
    % No hay scans en la ventana alrededor del RT estimado: dejamos charges>1 en 0 y salimos
    % fprintf('get_histone0: ventana vacía (charges) hno=%d rt1=%g rt2=%g cur_rt=%g RT=[%g,%g]\n', hno, rt1, rt2, cur_rts(1), min(rtcol), max(rtcol));
    return;
end

rt_i1 = p(1);
rt_i2 = pp(end);

for jno=2:ncharge
    % get MS1 profile
    c_mz = His.pep_mz(hno,jno);
    c_ch = His.pep_ch(hno,jno);
    c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
    if ptol>100 && c_ch>=3
        nC13 = 1;
    else
        nC13 = 0;
    end;
    [c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
    j = 2;
    c_mono_isointens = c_ref_isointens(:,j);

    % get rt and area
    [nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens);%#ok
    flag = zeros([1,length(nt)]);
    for i=1:length(nt)
        if abs( c_isorts(nt(i))-cur_rts(1) )<=delta
            flag(i) = 1;
        end;
    end;
    x = find(flag==1);
    if 0==isempty(x)
        [tmp,id] = max(inten_sum(x));%#ok
        top1_idx = x(id);
        cur_pos = nt(top1_idx);
        cur_rts(jno) = c_isorts(cur_pos);
        cur_intens(jno) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
    else
        cur_rts(jno) = cur_rts(1);
        cur_intens(jno) = 0;
    end;
end;

function ref_rt = get_H33_27_40_unmod_rt(His,gradient)
%%

out_file0 = fullfile(His.outpath,'H3_04_27_40.xls');
if 0~=exist(out_file0,'file')
    fp = fopen(out_file0,'r');
    str = fgetl(fp);
    while 0==feof(fp) && 0==strcmp(str,'[rt]')
        str = fgetl(fp);
    end;
    str = fgetl(fp);%#ok
    str = fgetl(fp);
    fclose(fp);
    p = strfind(str,'	');
    rt0 = str2num( str(p(1)+1:p(2)-1) );%#ok
    if 1==gradient
        ref_rt = rt0-0.3;
    else
        ref_rt = rt0-1.6;
    end;
else
    ref_rt = His.rt_ref(1);
end;
