function [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special)
%%

[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([2,ncharge]);
cur_intens = zeros([2,ncharge]);
num_MS1 = size(MS1_index,1);
cur_mono_isointens = zeros([num_MS1,2]);

[h_rts,h_intens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if 0==h_intens(1)
    return;
end;

[s_rts,ratio,cur_mono_isointens] = get_ratio_2iso(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,h_rts,special);
cur_rts(1,1:ncharge) = repmat(s_rts(1),[1,ncharge]);
cur_rts(2,1:ncharge) = repmat(s_rts(2),[1,ncharge]);
cur_intens(1,1:ncharge) = h_intens*ratio;
cur_intens(2,1:ncharge) = h_intens*(1-ratio);

function [s_rts,ratio,cur_mono_isointens] = get_ratio_2iso(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,h_rts,special)
%%

npart = 2;
s_rts = zeros([1,npart]);
ratio = 0;
num_MS1 = size(MS1_index,1);
cur_mono_isointens = zeros([num_MS1,npart]);

delta = 0.5;
p = find( MS1_index(:,2)>=His.rt_ref(hno)-delta );
rt_i1 = p(1);
pp = find( MS1_index(:,2)<=His.rt_ref(hno)+delta );
rt_i2 = pp(end);

% c_mono_isointens
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
[c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
j = 2;
c_mono_isointens = c_ref_isointens(:,j);

% rt1,rt2
IX = rt_i1:rt_i2;
% [nt,nb] = GetTopBottom(c_mono_isointens);%#ok
% p = find(c_isorts<=h_rts(1));
% c_ms1pos = p(end);
% [localmax_rt,localmax_inten,IX] = GetLocal(c_ms1pos,c_isorts,c_mono_isointens,nb);%#ok
rt1 = c_isorts(IX(1));
rt2 = c_isorts(IX(end));
if 2==special.nDAmode
    premzs = unique(MS2_index(:,4));
    [tmp,ii] = min( abs(premzs-c_mz) );%#ok
    target = premzs(ii);
    num_MS2 = size(MS2_index,1);
    flag = zeros([num_MS2,1]);
    p = find( MS2_index(:,2)>=rt1 );
    if 1==isempty(p)
        return;
    end;
    i1 = p(1);
    pp = find( MS2_index(:,2)<=rt2 );
    if 1==isempty(pp)
        return;
    end;
    i2 = pp(end);
    for i=i1:i2
        cen_mz = MS2_index(i,4);
        if 0==cen_mz-target
            flag(i) = 1;
        end;
    end;
    
    % check MS2
    ms2pos = find(flag==1);
    if 1==isempty(ms2pos)
        return;
    end;
    IX = zeros([length(ms2pos),1]);
    for i=1:length(ms2pos)
        IX(i) = find(MS1_index(:,1)==MS2_index(ms2pos(i),1));
    end;
    rt1 = c_isorts(IX(1));
    rt2 = c_isorts(IX(end));
end;

% match MS2
[ms2pos,ms2ratios] = MatchMS2(MS2_index,MS2_peaks,c_mz,c_ch,ptol,unitdiff,Mods,His,hno,rt1,rt2,special);

% ratios for profile
nlen = length(IX);
ratio1 = zeros([nlen,1]);
if length(ms2pos)<nlen/2
    if nlen<=3
        return;
    end;
    for i=1:nlen
        ratio1(i) = cos((i-1)*pi/(2*nlen));% cosine
    end;
    if His.rt_ref(hno)>His.rt_ref(hno+1)
        ratio1 = rot90(rot90(ratio1));
    end;
    ratio1 = ratio1/10;
    ratio2 = 1-ratio1;
else
    ms1scans = MS2_index(ms2pos,1);
    for i=1:nlen
        c_scan = MS1_index(IX(i),1);
        [tf,loc] = ismember(c_scan,ms1scans);
        if 1==tf
            ratio1(i) = ms2ratios(loc);
        end;
    end;
    if His.rt_ref(hno)<=His.rt_ref(hno+1)
        [tmp,x1] = max(ratio1);
        nX1 = x1:nlen;
        x2 = find(ratio1(nX1)<0.1);
        if tmp>0.1 && 0==isempty(x2)
            nX2 = nX1(x2(1)):nlen;
            ratio1(nX2) = 0.001;
        end;
    else
        [tmp,x1] = max(ratio1);
        nX1 = 1:x1;
        x2 = find(ratio1(nX1)<0.1);
        if tmp>0.1 && 0==isempty(x2)
            nX2 = 1:nX1(x2(end));
            ratio1(nX2) = 0.001;
        end;
    end;
    ratio1 = smooth(ratio1,3);
    ratio2 = 1-ratio1;
end;

n_rt = c_isorts(IX);
n_inten1 = c_mono_isointens(IX);
n_inten2 = n_inten1.*ratio1;
n_inten3 = n_inten1.*ratio2;

new_inten2 = smooth(n_inten2,3);
new_inten3 = smooth(n_inten3,3);
new_inten2(find(new_inten2<0.2*max(new_inten2))) = 0;%#ok
new_inten3(find(new_inten3<0.2*max(new_inten3))) = 0;%#ok
[nt2,nb2,top1_idx2] = GetTopBottom11(new_inten2);%#ok
[nt3,nb3,top1_idx3] = GetTopBottom11(new_inten3);%#ok
if 1==isempty(nt2)
    top1 = His.rt_ref(hno);
else
    top1 = n_rt(nt2(top1_idx2));
end;
if 1==isempty(nt3)
    top2 = His.rt_ref(hno);
else
    top2 = n_rt(nt3(top1_idx3));
end;

xx = n_rt(1):0.005:n_rt(end);
yy = spline(n_rt,n_inten2,xx);
area1 = sum(abs(yy));

xx = n_rt(1):0.005:n_rt(end);
yy = spline(n_rt,n_inten3,xx);
area2 = sum(abs(yy));

s_rts = [top1 top2];
ratio = area1/(eps+area1+area2);
if 2~=special.nDAmode
    cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]);n_inten2;zeros([num_MS1-IX(end),1])];
    cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]);n_inten3;zeros([num_MS1-IX(end),1])];
else
    xx = c_isorts(IX(1):IX(end));
    yy = spline(n_rt,n_inten2,xx);
    cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
    yy = spline(n_rt,n_inten3,xx);
    cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
end;
%
set(gcf,'visible','off');
out_file1 = fullfile(His.outpath,['Iso_',His.outfile,'_',His.mod_short{hno},'_',His.mod_short{hno+1},'.pdf']);
plot(n_rt,n_inten1,'linestyle','-','linewidth',2,'color','k');
hold on;
plot(n_rt,n_inten2,'linestyle','--','linewidth',2,'color','r');
plot(n_rt,n_inten3,'linestyle','-.','linewidth',2,'color','b');
xlabel('time (min)');
ylabel('intensity');
legend('experiment',His.mod_short{hno},His.mod_short{hno+1});
if area1/(eps+area1+area2)<0.03
    r1 = floor(1000*area1/(eps+area1+area2))*0.1;
else
    r1 = floor(100*area1/(eps+area1+area2));
end;
r2 = 100-r1;
title([His.mod_short{hno},'/',His.mod_short{hno+1},':',num2str(r1),'%:',num2str(r2),'%']);
print('-dpdf',out_file1);
close();
%}

function [ms2pos,ms2ratios] = MatchMS2(MS2_index,MS2_peaks,c_mz,c_ch,ptol,unitdiff,Mods,His,hno,rt1,rt2,special)
%%

% nhmass
if 4==special.nsource && (0~=special.nsubtype && 2~=special.nsubtype)
    nhmass = 1;
else
    nhmass = 0;
end;

% get precursors in MS1 profile
num_MS2 = size(MS2_index,1);
if 1==special.nDAmode
    % DDA
    sets = [0 1 2];
    mzs = c_mz + sets*unitdiff/c_ch;
    flag = zeros([num_MS2,1]);
    p = find( MS2_index(:,2)>=rt1 );
    if 1==isempty(p)
        ms2pos = [];
        ms2ratios = [];
        return;
    end;
    i1 = p(1);
    pp = find( MS2_index(:,2)<=rt2 );
    if 1==isempty(pp)
        ms2pos = [];
        ms2ratios = [];
        return;
    end;
    i2 = pp(end);
    for i=i1:i2
        cen_mz = MS2_index(i,4);
        ix = find(abs(mzs-cen_mz)<ptol*cen_mz*1e-6);%#ok
        if 0==isempty(ix)
            flag(i) = 1;
        end;
    end;

    % check MS2
    ms2pos = find(flag==1);
    if 1==isempty(ms2pos)
        ms2ratios = [];
        return;
    end;
elseif 2==special.nDAmode
    % DIA
    premzs = unique(MS2_index(:,4));
    [tmp,ii] = min( abs(premzs-c_mz) );%#ok
    target = premzs(ii);
    flag = zeros([num_MS2,1]);
    p = find( MS2_index(:,2)>=rt1 );
    if 1==isempty(p)
        ms2pos = [];
        ms2ratios = [];
        return;
    end;
    i1 = p(1);
    pp = find( MS2_index(:,2)<=rt2 );
    if 1==isempty(pp)
        ms2pos = [];
        ms2ratios = [];
        return;
    end;
    i2 = pp(end);
    for i=i1:i2
        cen_mz = MS2_index(i,4);
        if 0==cen_mz-target
            flag(i) = 1;
        end;
    end;

    % check MS2
    ms2pos = find(flag==1);
    if 1==isempty(ms2pos)
        ms2ratios = [];
        return;
    end;
else
    ms2pos = [];
    ms2ratios = [];
    return;
end;

instruments = MS2_index(ms2pos,6);% MS2dirs = {'CIDIT','CIDFT','ETDIT','ETDFT','HCDIT','HCDFT'};
if 1==length(unique(instruments))
    % ActiveType, tol
    c_instrument = instruments(1);
    if 3==c_instrument || 4==c_instrument
        ActiveType = 'ETD';
    else
        ActiveType = 'CID';
    end;
    if 1==mod(c_instrument,2)
        tol = 0.4;
    else
        tol = 0.02;
    end;

    % K1,K2
    if 1==nhmass
        [K1,K2] = get_key_ionsH(His,hno,hno+1,Mods,ActiveType);
    else
        [K1,K2] = get_key_ions(His,hno,hno+1,Mods,ActiveType);
    end;
end;

index = [1;MS2_index(1:num_MS2,7)];
ms2ratios = zeros([1,length(ms2pos)]);
for i=1:length(ms2pos)
    cno = ms2pos(i);
%     if 2==special.nDAmode
%         newpos = cno-1:cno+1;
%     else
%         newpos = cno;
%     end;
    newpos = cno;
    for pno = newpos
        if pno<1 || pno>num_MS2
            continue;
        end;
        if 1<length(unique(instruments))
            % ActiveType, tol
            c_instrument = MS2_index(pno,6);% MS2dirs = {'CIDIT','CIDFT','ETDIT','ETDFT','HCDIT','HCDFT'};
            if 3==c_instrument || 4==c_instrument
                ActiveType = 'ETD';
            else
                ActiveType = 'CID';
            end;
            if 1==mod(c_instrument,2)
                tol = 0.4;
            else
                tol = 0.02;
            end;
        end;

        if 1<length(unique(instruments))
            % K1,K2
            if 1==nhmass
                [K1,K2] = get_key_ionsH(His,hno,hno+1,Mods,ActiveType);
            else
                [K1,K2] = get_key_ions(His,hno,hno+1,Mods,ActiveType);
            end;
        end;

        % mz, inten
        IX = index(pno):index(pno+1)-1;
        mz = MS2_peaks(IX,1);
        inten = MS2_peaks(IX,2);

        % match key ions
        intens1 = zeros([1,length(K1)]);
        intens2 = zeros([1,length(K1)]);
        for j=1:length(K1)
            ix1 = find(abs(mz-K1(j))<=tol);
            ix2 = find(abs(mz-K2(j))<=tol);
            [tmp,x1] = min(abs(mz(ix1)-K1(j)));%#ok
            [tmp,x2] = min(abs(mz(ix2)-K2(j)));%#ok
            %[tmp,x1] = max(inten(ix1));
            %[tmp,x2] = max(inten(ix2));
            if 0==isempty(ix1)
                intens1(j) = inten(ix1(x1));
            end;
            if 0==isempty(ix2)
                intens2(j) = inten(ix2(x2));
            end;
        end;
        %r1 = intens1./(eps+intens1+intens2);
        %ratio1 = median(r1);
        ratio1 = sum(intens1)/(eps+sum(intens1)+sum(intens2));
        if ms2ratios(i)<ratio1
            ms2ratios(i) = ratio1;
        end;
    end;
end;