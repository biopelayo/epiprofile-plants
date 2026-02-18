function [cur_rts,cur_intens,cur_mono_isointens] = get_histone10(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
%%

[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([1,ncharge]);
cur_intens = zeros([1,ncharge]);

delta = 2;
p = find( MS1_index(:,2)>=His.rt_ref(hno)-delta );
rt_i1 = p(1);
pp = find( MS1_index(:,2)<=His.rt_ref(hno)+delta );
rt_i2 = pp(end);

if ptol==100
    ptol = 10;
end;

for jno=1:ncharge
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
    if 1==jno
        cur_mono_isointens = c_mono_isointens;
    end;

    % get rt and area
    [nt,nb,top1_idx,inten_sum] = GetTopBottom11(c_mono_isointens);%#ok
    flag = zeros([1,length(nt)]);
    for i=1:length(nt)
        if abs( c_isorts(nt(i))-His.rt_ref(hno) )<=0.5
            flag(i) = 1;
        end;
    end;
    x = find(flag==1);
    if 0==isempty(x)
        [tmp,id] = min(abs(c_isorts(nt(x))-His.rt_ref(hno)));%#ok
        %[tmp,id] = max(inten_sum(x));%#ok
        top1_idx = x(id);
        cur_pos = nt(top1_idx);
        cur_rts(jno) = c_isorts(cur_pos);
        if length(x)>=2
            nb = [nb(x(1)) nb(x(end)+1)];
        end;
        cur_intens(jno) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
    else
        cur_rts(jno) = His.rt_ref(hno);
        cur_intens(jno) = 0;
    end;
end;