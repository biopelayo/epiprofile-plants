function [top1_rt,top1_inten_sum] = get_rts0(MS1_index,MS1_peaks,ptol,His,hno)
%%

unitdiff = 1.0032;
num_MS1 = size(MS1_index,1);

if ptol==100
    ptol = 10;
end;
if ptol>100 && His.pep_ch(hno)>=3
    nC13 = 1;
else
    nC13 = 0;
end;

str = {'TKQTAR','KSTGGKAPR','KSAPATGGVKKPHR','KVLR'};
if 1==ismember(His.pep_seq{hno,1},str)
    c_ch = His.pep_ch(hno);
    c_mz = His.pep_mz(hno)+14.01565/c_ch;
    c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
    [c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,1:num_MS1);
    c_mono_isointens = c_ref_isointens(:,2);
    [nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens);%#ok
    top1_rt = c_isorts(nt(top1_idx));
    if 1==isempty(top1_rt) || top1_rt<4
        top1_rt = 0;
    end;
    if 0==top1_rt
        new_pos = num_MS1;
    else
        p = find( MS1_index(:,2)<=top1_rt );
        new_pos = p(end);
    end;
else
    new_pos = num_MS1;
end;

% get MS1 profile
c_mz = His.pep_mz(hno);
c_ch = His.pep_ch(hno);

c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
[c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,1:new_pos);
c_mono_isointens = c_ref_isointens(:,2);
[nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens);%#ok
top1_rt = c_isorts(nt(top1_idx));
top1_inten_sum = inten_sum(top1_idx);

if 1==isempty(top1_rt) || top1_rt<4
    top1_rt = 0;
    top1_inten_sum = 0;
end;