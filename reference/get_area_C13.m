function pep_area = get_area_C13(c_isorts,c_ref_isointens,nb,cur_pos,cur_mz,cur_ch,MS1_index,MS1_peaks,unitdiff,ptol)
%%

% init
pep_area = 0;
c_mono_isointens = c_ref_isointens(:,2);

nw = c_mono_isointens;
% nw = smooth(nw,3);
nw(find(nw<0.015*max(nw))) = 0;%#ok
if 0==nw(cur_pos)
    nw = c_mono_isointens;
end;

bflag = 0;
for ino=1:length(nb)-1
    i1 = nb(ino);
    i2 = nb(ino+1);
    while i1<=length(nw) && 0==nw(i1)
        i1 = i1+1;
    end;
    if i1>nb(ino)
        i1 = i1-1;
    end;
    while i2>=1 && 0==nw(i2)
        i2 = i2-1;
    end;
    if i2<nb(ino+1)
        i2 = i2+1;
    end;
    IX = i1:i2;
    nmid = find( IX==cur_pos );%#ok
    if 0==isempty( nmid )
        bflag = 1;
        break;
    end;
end;
if 0==bflag
    return;
end;

% n_rt, n_inten
n_rt = c_isorts(IX);
n_inten = c_mono_isointens(IX);
n_inten = smooth(n_inten,3);

% area
nOverlap = 1; % 1: yes, 0: no
if 1==nOverlap
    pep_area = judgeOverlap2(cur_mz,cur_ch,IX,cur_pos,MS1_index,MS1_peaks,unitdiff,ptol);
else
    xx = n_rt(1):0.005:n_rt(end);
    yy = spline(n_rt,n_inten,xx);
    pep_area = sum(abs(yy));
end;
if pep_area<0
    pep_area = 0;
end;

% Sum the isotopic clusters
nSumIso = 1; % 1: yes, 0: no
if 1==nSumIso
    IPV = GetIPV();
    M = cur_mz*cur_ch-cur_ch*1.007276;
    tIPV = IPV(floor(M),1:5);
    fCali = sum(tIPV);
else
    fCali = 1;
end;
pep_area = pep_area*fCali;

function pep_area = judgeOverlap2(cur_mz,cur_ch,IX,cur_pos,MS1_index,MS1_peaks,unitdiff,ptol)
%%

% mz, inten
num_MS1 = size(MS1_index,1);
index = [1;MS1_index(1:num_MS1,3)];
I = index(cur_pos):index(cur_pos+1)-1;
mz = MS1_peaks(I,1);
inten = MS1_peaks(I,2);

% cur_mz as the center
tmp_isomzs = [cur_mz-9*unitdiff/cur_ch cur_mz-8*unitdiff/cur_ch cur_mz-7*unitdiff/cur_ch cur_mz-6*unitdiff/cur_ch cur_mz-5*unitdiff/cur_ch cur_mz-4*unitdiff/cur_ch cur_mz-3*unitdiff/cur_ch cur_mz-2*unitdiff/cur_ch cur_mz-unitdiff/cur_ch cur_mz];
nmz = length(tmp_isomzs);
tmp_intens = repmat(0,[1,nmz]);
for jno=1:nmz
    c_mz = tmp_isomzs(jno);
    c_ptol = ptol*c_mz*1e-6;
    left = c_mz - c_ptol;
    right = c_mz + c_ptol;
    pos = find( mz>=left & mz<=right );
    if 0==isempty(pos)
        tmp_intens(jno) = max(inten(pos));
    end;
end;

IPV = GetIPV();
M = cur_mz*cur_ch-cur_ch*1.007276;
tIPV = IPV(floor(M),1:5);

% judge overlap
thr = 5;
t1 = tIPV(2);
a23 = tmp_intens(3)/(eps+tmp_intens(2));
a45 = tmp_intens(5)/(eps+tmp_intens(4));
a67 = tmp_intens(7)/(eps+tmp_intens(6));
a89 = tmp_intens(9)/(eps+tmp_intens(8));
r_at23 = a23/t1;
r_at45 = a45/t1;
r_at67 = a67/t1;
r_at89 = a89/t1;

%{
c12 = tmp_intens(1)<tmp_intens(2);
c34 = tmp_intens(3)<tmp_intens(4);
c56 = tmp_intens(5)<tmp_intens(6);
c78 = tmp_intens(7)<tmp_intens(8);
%}
cc23 = (r_at23>=1/thr && r_at23<=thr);
cc45 = (r_at45>=1/thr && r_at45<=thr);
cc67 = (r_at67>=1/thr && r_at67<=thr);
cc89 = (r_at89>=1/thr && r_at89<=thr);

if cc23 && cc45 && cc67 && cc89% && c12
    c_ref_isomzs = [cur_mz-8*unitdiff/cur_ch cur_mz-6*unitdiff/cur_ch cur_mz-4*unitdiff/cur_ch cur_mz-2*unitdiff/cur_ch cur_mz];
    idx = [2 4 6 8 10];
elseif cc45 && cc67 && cc89% && c34
    c_ref_isomzs = [cur_mz-6*unitdiff/cur_ch cur_mz-4*unitdiff/cur_ch cur_mz-2*unitdiff/cur_ch cur_mz];
    idx = [4 6 8 10];
elseif cc67 && cc89% && c56
    c_ref_isomzs = [cur_mz-4*unitdiff/cur_ch cur_mz-2*unitdiff/cur_ch cur_mz];
    idx = [6 8 10];
elseif cc89% && c78
    c_ref_isomzs = [cur_mz-2*unitdiff/cur_ch cur_mz];
    idx = [8 10];
else
    c_ref_isomzs = cur_mz;
    idx = 10;
end;

% get the final area
[c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,cur_ch,ptol,1,IX);
n_rt = c_isorts(IX);
xx = n_rt(1):0.005:n_rt(end);

nmz = length(c_ref_isomzs);
tmp_areas = repmat(0,[1,nmz]);
for jno=1:nmz
    n_inten = c_ref_isointens(IX,jno);
    n_inten = smooth(n_inten,3);
    yy = spline(n_rt,n_inten,xx);
    tmp_areas(jno) = sum(abs(yy));
end;

Factor = tmp_areas(end)/(eps+tmp_intens(end));
for jno=1:nmz
    tmp_areas(jno) = min(tmp_areas(jno),Factor*(eps+tmp_intens(idx(jno))));
end;

pep_area = tmp_areas(1);
if 1==nmz
    return;
end;
for jno=2:nmz
    id = idx(jno)-idx(jno-1)+1;
    pep_area = tmp_areas(jno)-pep_area*tIPV(id);
end;