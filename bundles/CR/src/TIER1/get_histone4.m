function [cur_rts,cur_intens,cur_mono_isointens] = get_histone4(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special)
%%

npart = 4;
[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([npart,ncharge]);
cur_intens = zeros([npart,ncharge]);
num_MS1 = size(MS1_index,1);
cur_mono_isointens = zeros([num_MS1,npart]);

[h_rts,h_intens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if 0==h_intens(1)
    return;
end;

[s_rts,ratio,cur_mono_isointens] = get_ratio_4iso(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,h_rts,special);
cur_rts(1,1:ncharge) = repmat(s_rts(1),[1,ncharge]);
cur_rts(2,1:ncharge) = repmat(s_rts(2),[1,ncharge]);
cur_rts(3,1:ncharge) = repmat(s_rts(3),[1,ncharge]);
cur_rts(4,1:ncharge) = repmat(s_rts(4),[1,ncharge]);
cur_intens(1,1:ncharge) = h_intens*ratio(1);
cur_intens(2,1:ncharge) = h_intens*ratio(2);
cur_intens(3,1:ncharge) = h_intens*ratio(3);
cur_intens(4,1:ncharge) = h_intens*ratio(4);

function [s_rts,ratio,cur_mono_isointens] = get_ratio_4iso(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,h_rts,special)%#ok
%%

npart = 4;
s_rts = zeros([1,npart]);
ratio = zeros([1,npart]);
num_MS1 = size(MS1_index,1);
cur_mono_isointens = zeros([num_MS1,npart]);

delta = 0.5;
p = find( MS1_index(:,2)>=His.rt_ref(hno)-delta );
rt_i1 = p(1);
pp = find( MS1_index(:,2)<=His.rt_ref(hno+3)+delta );
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
% [nt,nb] = GetTopBottom(c_mono_isointens);
% p = find(c_isorts<=h_rts(1));
% c_ms1pos = p(end);
% [localmax_rt,localmax_inten,IX] = GetLocal(c_ms1pos,c_isorts,c_mono_isointens,nb);
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
if length(ms2pos)<=3
    if nlen<=3
        return;
    end;
    s_rts(1:npart) = [His.rt_ref(hno) His.rt_ref(hno+1) His.rt_ref(hno+2) His.rt_ref(hno+3)];
    ratio(1:npart) = repmat(1/npart,[1,npart]);
    n_inten1 = c_mono_isointens(IX)/npart;
    if 2~=special.nDAmode
        cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]);n_inten1;zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]);n_inten1;zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,3) = [zeros([IX(1)-1,1]);n_inten1;zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,4) = [zeros([IX(1)-1,1]);n_inten1;zeros([num_MS1-IX(end),1])];
    else
        xx = c_isorts(IX(1):IX(end));
        yy = spline(c_isorts(IX),n_inten1,xx);
        cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,3) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,4) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
    end;
    return;
end;

ratio1 = zeros([nlen,npart]);
ms1scans = MS2_index(ms2pos,1);
for i=1:nlen
    c_scan = MS1_index(IX(i),1);
    [tf,loc] = ismember(c_scan,ms1scans);
    if 1==tf
        ratio1(i,1:npart) = ms2ratios(1:npart,loc)';
    end;
end;
ratio1(:,1) = smooth(ratio1(:,1),3);
ratio1(:,2) = smooth(ratio1(:,2),3);
ratio1(:,3) = smooth(ratio1(:,3),3);
ratio1(:,4) = smooth(ratio1(:,4),3);
sumrow = sum(ratio1,2);
ratio1(:,1) = ratio1(:,1)./(eps+sumrow);
ratio1(:,2) = ratio1(:,2)./(eps+sumrow);
ratio1(:,3) = ratio1(:,3)./(eps+sumrow);
ratio1(:,4) = ratio1(:,4)./(eps+sumrow);

n_rt = c_isorts(IX);
n_inten1 = c_mono_isointens(IX);
n_inten2 = n_inten1.*ratio1(:,1);
n_inten3 = n_inten1.*ratio1(:,2);
n_inten4 = n_inten1.*ratio1(:,3);
n_inten5 = n_inten1.*ratio1(:,4);

[tmp,x1] = max(n_inten2);%#ok
[tmp,x2] = max(n_inten3);%#ok
[tmp,x3] = max(n_inten4);%#ok
[tmp,x4] = max(n_inten5);%#ok

xx = n_rt(1):0.005:n_rt(end);
yy = spline(n_rt,n_inten2,xx);
area1 = sum(abs(yy));

xx = n_rt(1):0.005:n_rt(end);
yy = spline(n_rt,n_inten3,xx);
area2 = sum(abs(yy));

xx = n_rt(1):0.005:n_rt(end);
yy = spline(n_rt,n_inten4,xx);
area3 = sum(abs(yy));

xx = n_rt(1):0.005:n_rt(end);
yy = spline(n_rt,n_inten5,xx);
area4 = sum(abs(yy));
sumarea = eps+area1+area2+area3+area4;

s_rts(1:npart) = [n_rt(x1) n_rt(x2) n_rt(x3) n_rt(x4)];
ratio(1:npart) = [area1/sumarea area2/sumarea area3/sumarea area4/sumarea];
if 2~=special.nDAmode
    cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]);n_inten2;zeros([num_MS1-IX(end),1])];
    cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]);n_inten3;zeros([num_MS1-IX(end),1])];
    cur_mono_isointens(:,3) = [zeros([IX(1)-1,1]);n_inten4;zeros([num_MS1-IX(end),1])];
    cur_mono_isointens(:,4) = [zeros([IX(1)-1,1]);n_inten5;zeros([num_MS1-IX(end),1])];
else
    xx = c_isorts(IX(1):IX(end));
    yy = spline(n_rt,n_inten2,xx);
    cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
    yy = spline(n_rt,n_inten3,xx);
    cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
    yy = spline(n_rt,n_inten4,xx);
    cur_mono_isointens(:,3) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
    yy = spline(n_rt,n_inten5,xx);
    cur_mono_isointens(:,4) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
end;

%
set(gcf,'visible','off');
xx = strfind(His.mod_short{hno},'ac');
out_file1 = fullfile(His.outpath,['Iso_',His.outfile,'_',num2str(length(xx)),'ac','_',His.mod_short{hno},'.pdf']);
plot(n_rt,n_inten1,'linestyle','-','linewidth',2,'color','k');
hold on;
plot(n_rt,n_inten2,'linestyle','--','linewidth',2,'color','r');
plot(n_rt,n_inten3,'linestyle','-.','linewidth',2,'color','b');
plot(n_rt,n_inten4,'linestyle','--','linewidth',2,'color','g');
plot(n_rt,n_inten5,'linestyle','-.','linewidth',2,'color','m');
xlabel('time (min)');
ylabel('intensity');
legend('experiment',His.mod_short{hno},His.mod_short{hno+1},His.mod_short{hno+2},His.mod_short{hno+3});
r1 = floor(100*ratio(1));
r2 = floor(100*ratio(2));
r3 = floor(100*ratio(3));
r4 = 100-r1-r2-r3;
title([His.mod_short{hno},'/',His.mod_short{hno+1},'/',His.mod_short{hno+2},'/',His.mod_short{hno+3},':',num2str(r1),'%:',num2str(r2),'%:',num2str(r3),'%:',num2str(r4),'%']);
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

npart = 4;
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
        [K11,K12] = get_key_ionsH(His,hno,hno+1,Mods,ActiveType);
        [K21,K22] = get_key_ionsH(His,hno+1,hno+2,Mods,ActiveType);
        [K31,K32] = get_key_ionsH(His,hno+2,hno+3,Mods,ActiveType);
    else
        [K11,K12] = get_key_ions(His,hno,hno+1,Mods,ActiveType);
        [K21,K22] = get_key_ions(His,hno+1,hno+2,Mods,ActiveType);
        [K31,K32] = get_key_ions(His,hno+2,hno+3,Mods,ActiveType);
    end;
end;

index = [1;MS2_index(1:num_MS2,7)];
ms2ratios = zeros([npart,length(ms2pos)]);
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
                [K11,K12] = get_key_ionsH(His,hno,hno+1,Mods,ActiveType);
                [K21,K22] = get_key_ionsH(His,hno+1,hno+2,Mods,ActiveType);
                [K31,K32] = get_key_ionsH(His,hno+2,hno+3,Mods,ActiveType);
            else
                [K11,K12] = get_key_ions(His,hno,hno+1,Mods,ActiveType);
                [K21,K22] = get_key_ions(His,hno+1,hno+2,Mods,ActiveType);
                [K31,K32] = get_key_ions(His,hno+2,hno+3,Mods,ActiveType);
            end;
        end;

        % mz, inten
        IX = index(pno):index(pno+1)-1;
        mz = MS2_peaks(IX,1);
        inten = MS2_peaks(IX,2);

        % match key ions
        a = get_ratio(mz,inten,tol,K11,K12);
        ab = get_ratio(mz,inten,tol,K21,K22);
        abc = get_ratio(mz,inten,tol,K31,K32);
        d = 1-abc;
        b = max([ab-a 0]);
        c = 1-a-b-d;
        if a+d>1
            d = 1-a;
            b = 0;
            c = 0;
        elseif a+d+b>1
            b = 1-a-d;
            c = 0;
        end;
        X = [a;b;c;d];
        if length(find(ms2ratios(1:npart,i)==0))>length(find(X==0))
            ms2ratios(1:npart,i) = X;
        end;
    end;
end;

% no fragments
x = find(ms2ratios(npart,:)<1);%#ok
if 1==isempty(x)
    ms2pos = [];
    ms2ratios = [];
else
    xx = find(ms2ratios(npart,:)==1);
    if 0==isempty(xx)
        for ino=1:length(xx)
            ms2ratios(1:npart,xx(ino)) = [0 0 0.01 0.99]';
        end;
    end;
end;

function ratio1 = get_ratio(mz,inten,tol,K1,K2)
%%

intens1 = zeros([1,length(K1)]);
intens2 = zeros([1,length(K1)]);
for j=1:length(K1)
    ix1 = find(abs(mz-K1(j))<=tol);
    ix2 = find(abs(mz-K2(j))<=tol);
    [tmp,x1] = min(abs(mz(ix1)-K1(j)));%#ok
    [tmp,x2] = min(abs(mz(ix2)-K2(j)));%#ok
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