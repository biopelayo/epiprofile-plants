function [rts,top1_rt,inten_sum,top1_inten_sum] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass)
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

p = find( MS1_index(:,2)>=t1 );
rt_i1 = p(1);
pp = find( MS1_index(:,2)<=t2 );
rt_i2 = pp(end);

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
    [nt,nb,top1_idx,inten_sum] = GetTopBottom11(c_mono_isointens);
else
    [nt,nb,top1_idx,inten_sum] = GetTopBottom11(c_mono_isointens);
end;
rts = c_isorts(nt);
% top1_rt = c_isorts(nt(top1_idx));
% top1_inten_sum = inten_sum(top1_idx);

%+
Mods = GetMods();
flen = 0;
% w = smooth(c_mono_isointens,3);
w = c_mono_isointens;
w(find(w<0.015*max(w))) = 0;%#ok
similarity = zeros([length(nt),1]);
for ino=1:length(nt)
    i1 = nb(ino);
    i2 = nb(ino+1);
    while i1<=length(w) && 0==w(i1)
        i1 = i1+1;
    end;
    if i1>nb(ino)
        i1 = i1-1;
    end;
    while i2>=1 && 0==w(i2)
        i2 = i2-1;
    end;
    if i2<nb(ino+1)
        i2 = i2+1;
    end;
    IX = i1:i2;
    rt1 = c_isorts(IX(1));
    rt2 = c_isorts(IX(end));
    
    [ms2pos,ms2rts,ms2intens,posn,posc] = MatchMS2(MS2_index,MS2_peaks,Mods,His,hno,rt1,rt2,nhmass);
    if flen<length(posn)
        flen = length(posn);
    end;
    if 1==isempty(ms2pos)
        continue;
    end;
    
    nstep = median(ms2rts(ms2pos(2:end))-ms2rts(ms2pos(1:end-1)))*3.5;
    nst = c_isorts(nt(ino))-nstep;
    ntm = c_isorts(nt(ino))+nstep;
    np1 = find(ms2rts(ms2pos)>=nst);
    np2 = find(ms2rts(ms2pos)<=ntm);
    X = np1(1):np2(end);
    XX = zeros([1,length(X)]);
    for jno=1:length(X)
        XX(jno) = find(MS1_index(:,1)==MS2_index(ms2pos(X(jno)),1));
    end;
    p_intens = c_mono_isointens(XX);
    p_intens = smooth(p_intens,3);
    [sim_thr,ex_intens] = get_thr(p_intens);
    p_intens = p_intens+ex_intens;
    for kno=1:length(posn)
        f_intens = ms2intens(ms2pos(X),kno);
        f_intens = smooth(f_intens,3)+ex_intens;
        e_sim = sum(p_intens.*f_intens)/sqrt( sum(p_intens.*p_intens)*sum(f_intens.*f_intens) );
        if e_sim>sim_thr
            similarity(ino) = similarity(ino) + e_sim;
        end;
    end;
    
    for kno=1:length(posc)
        qno = kno+length(posn);
        f_intens = ms2intens(ms2pos(X),qno);
        f_intens = smooth(f_intens,3)+ex_intens;
        e_sim = sum(p_intens.*f_intens)/sqrt( sum(p_intens.*p_intens)*sum(f_intens.*f_intens) );
        if e_sim>sim_thr
            similarity(ino) = similarity(ino) + e_sim;
        end;
    end;
end;

tmp = max(ceil(similarity));
if tmp>flen/2
%     if (ceil(similarity(top1_idx))==tmp-1 && tmp>1) || (ceil(similarity(top1_idx))==tmp-2 && tmp>2)
%         similarity(top1_idx) = tmp;
%     end;
    ii = find(ceil(similarity)==tmp);
    [t,id] = max(inten_sum(ii));%#ok
    top1_idx = ii(id);
    top1_rt = c_isorts(nt(top1_idx));
    top1_inten_sum = inten_sum(top1_idx);
else
    top1_rt = c_isorts(nt(top1_idx));
    top1_inten_sum = inten_sum(top1_idx)/1e4;
end;
%+

if top1_rt<4
    rts = [];
    top1_rt = [];
    inten_sum = [];
    top1_inten_sum = [];
    return;
end;

function [sim_thr,ex_intens] = get_thr(p_intens0)
%%

x(1:length(p_intens0),1) = (ones(1,length(p_intens0)))'*eps;
x(1:length(p_intens0),2) = (1:length(p_intens0))'*eps;
x(1:length(p_intens0),3) = (length(p_intens0):-1:1)'*eps;

p_intens = p_intens0 + x(:,1);
f_intens = x(:,1);
sim_thrs(1) = sum(p_intens.*f_intens)/sqrt( sum(p_intens.*p_intens)*sum(f_intens.*f_intens) );

p_intens = p_intens0 + x(:,2);
f_intens = x(:,2);
sim_thrs(2) = sum(p_intens.*f_intens)/sqrt( sum(p_intens.*p_intens)*sum(f_intens.*f_intens) );

p_intens = p_intens0 + x(:,3);
f_intens = x(:,3);
sim_thrs(3) = sum(p_intens.*f_intens)/sqrt( sum(p_intens.*p_intens)*sum(f_intens.*f_intens) );

[sim_thr,i] = min(sim_thrs);
ex_intens = x(:,i);

function [ms2pos,ms2rts,ms2intens,posn,posc,ActiveType,K1] = MatchMS2(MS2_index,MS2_peaks,Mods,His,hno,rt1,rt2,nhmass)
%%

% get precursors in MS1 profile
num_MS2 = size(MS2_index,1);
c_mz = His.pep_mz(hno,1);
premzs = unique(MS2_index(:,4));
[tmp,ii] = min( abs(premzs-c_mz) );%#ok
target = premzs(ii);
flag = zeros([num_MS2,1]);
p = find( MS2_index(:,2)>=rt1 );
pp = find( MS2_index(:,2)<=rt2 );
if 1==isempty(p) || 1==isempty(pp)
    ms2pos = [];
    ms2rts = [];
    ms2intens = [];
    posn = [];
    posc = [];
    ActiveType = [];
    K1 = [];
    return;
end;
i1 = p(1);
i2 = pp(end);
for i=i1:i2
    cen_mz = MS2_index(i,4);
    if 0==cen_mz-target
        flag(i) = 1;
    end;
end;
ms2pos = find(flag==1);
if length(ms2pos)<=4
    ms2pos = [];
    ms2rts = [];
    ms2intens = [];
    posn = [];
    posc = [];
    ActiveType = [];
    K1 = [];
    return;
end;

% check MS2
ms2rts = MS2_index(:,2);

instruments = MS2_index(ms2pos,6);% MS2dirs = {'CIDIT','CIDFT','ETDIT','ETDFT','HCDIT','HCDFT'};
% if 1==length(unique(instruments))
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
    [K1,posn,posc] = get_key_ions1H(His,hno,Mods,ActiveType);
else
    [K1,posn,posc] = get_key_ions1(His,hno,Mods,ActiveType);
end;
% end;

index = [1;MS2_index(1:num_MS2,7)];
ms2intens = zeros([num_MS2,length(K1)]);
for i=1:length(ms2pos)
    cno = ms2pos(i);
    for pno = cno;%cno-1:cno+1
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
                [K1,posn,posc] = get_key_ions1H(His,hno,Mods,ActiveType);
            else
                [K1,posn,posc] = get_key_ions1(His,hno,Mods,ActiveType);
            end;
        end;
        
        % mz, inten
        IX = index(pno):index(pno+1)-1;
        mz = MS2_peaks(IX,1);
        inten = MS2_peaks(IX,2);
        I = find(inten>=3*MS2_index(pno,8));%+
        mz = mz(I);
        inten = inten(I);%+
        
        % match key ions
        for j=1:length(K1)
            ix1 = find(abs(mz-K1(j))<=tol);
            [tmp,x1] = min(abs(mz(ix1)-K1(j)));%#ok
            %[tmp,x1] = max(inten(ix1));
            if 0==isempty(ix1) && ms2intens(cno,j)<inten(ix1(x1))
                ms2intens(cno,j) = inten(ix1(x1));
            end;
        end;
    end;
end;