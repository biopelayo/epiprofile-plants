function [K1,K2] = get_key_ions2(His,pno,qno,Mods,ActiveType)
%%

c_seq = His.pep_seq;
c_mod1 = His.mod_type{pno};
c_mod2 = His.mod_type{qno};

% get_mod_postype
[modpos1,modtype1] = get_mod_postype(c_seq,c_mod1,Mods);
[modpos2,modtype2] = get_mod_postype(c_seq,c_mod2,Mods);

newtype1 = modpos1*10000+modtype1;
newtype2 = modpos2*10000+modtype2;
comtype = intersect(newtype1,newtype2);
uniontype = union(newtype1,newtype2);
lefttype = setdiff(uniontype,comtype );
t1 = floor(min(lefttype)/10000);
if 0==t1
    t1 = 1;
end;
t2 = floor(max(lefttype)/10000);
nlen = 2*(t2-t1);
K1 = repmat(0,[1,nlen]);
K2 = repmat(0,[1,nlen]);

peplen = length(c_seq);
posn = t1:t2-1;
posc = (peplen-t2+1):(peplen-t1);

% get_theo_mz
[theo_mz1,theo_tp1] = get_theo_mz(c_seq,c_mod1,1,Mods,ActiveType);
[theo_mz2,theo_tp2] = get_theo_mz(c_seq,c_mod2,1,Mods,ActiveType);

n1 = find(theo_tp1<2000);
tpn1 = floor( (theo_tp1(n1)-1000)/10 );
c1 = find(theo_tp1>2000);
tpc1 = floor( (theo_tp1(c1)-2000)/10 );

n2 = find(theo_tp2<2000);
tpn2 = floor( (theo_tp2(n2)-1000)/10 );
c2 = find(theo_tp2>2000);
tpc2 = floor( (theo_tp2(c2)-2000)/10 );

%{
for ino=1:length(posn)
    [tf1,loc1] = ismember(posn(ino),tpn1);
    [tf2,loc2] = ismember(posn(ino),tpn2);
    if 1==tf1 && 1==tf2
        K1(ino) = theo_mz1(n1(loc1));
        K2(ino) = theo_mz2(n2(loc2));
    end;
end;
%}

for ino=1:length(posc)
    [tf1,loc1] = ismember(posc(ino),tpc1);
    [tf2,loc2] = ismember(posc(ino),tpc2);
    if 1==tf1 && 1==tf2
        K1(ino+length(posn)) = theo_mz1(c1(loc1));
        K2(ino+length(posn)) = theo_mz2(c2(loc2));
    end;
end;

x = find(K1>0);
K1 = K1(x);
K2 = K2(x);