function [K1,posn,posc] = get_key_ions1(His,pno,Mods,ActiveType)
%%

c_seq = His.pep_seq;
if 1==strcmp(c_seq,'unmod')
    c_seq = His.mod_short{pno};
    p = strfind(c_seq,'.');
    if 0==isempty(p)
        c_seq = c_seq(p(end)+1:end);
    end;
end;
c_mod1 = His.mod_type{pno};

% get_mod_postype
modpos1 = get_mod_postype(c_seq,c_mod1,Mods);
peplen = length(c_seq);
ix = find(modpos1>0);
if 1==isempty(ix) || modpos1(ix(end))-modpos1(ix(1))<=3
    t1 = 1;
    t2 = peplen;
else
    t1 = min( modpos1(ix) );
    t2 = max( modpos1(ix) );
end;
nlen = 2*(t2-t1);
K1 = zeros([1,nlen]);

posn = t1:t2-1;
posc = (peplen-t2+1):(peplen-t1);

% get_theo_mz
[theo_mz1,theo_tp1] = get_theo_mz(c_seq,c_mod1,1,Mods,ActiveType);

n1 = find(theo_tp1<2000);
tpn1 = floor( (theo_tp1(n1)-1000)/10 );
c1 = find(theo_tp1>2000);
tpc1 = floor( (theo_tp1(c1)-2000)/10 );

for ino=1:length(posn)
    [tf1,loc1] = ismember(posn(ino),tpn1);
    if 1==tf1
        K1(ino) = theo_mz1(n1(loc1));
    end;
end;

for ino=1:length(posc)
    [tf1,loc1] = ismember(posc(ino),tpc1);
    if 1==tf1
        K1(ino+length(posn)) = theo_mz1(c1(loc1));
    end;
end;