function [localmax_rt,localmax_inten,IX] = GetLocal(cur_ms1pos,isorts,c_ref_isointens,nb)
%%

nmz = size(c_ref_isointens,2);
localmax_rt = repmat(0,[nmz,1]);
localmax_inten = repmat(0,[nmz,1]);
IX = [];

nw = c_ref_isointens(:,1);
% nw = smooth(nw,3);
nw(find(nw<0.015*max(nw))) = 0;%#ok
if 0==nw(cur_ms1pos)
    nw = c_ref_isointens(:,1);
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
    c_pos = find( IX==cur_ms1pos );%#ok
    if 0==isempty( c_pos )
        bflag = 1;
        %{
        c_max = max(w(IX));
        [t,b] = JudgeLocalmaxmin(w(IX));
        if length(b)>2
            for j=2:length(b)-1
                mid = IX(b(j));
                if nw(mid)<0.2*c_max
                    s1 = sum( w(i1:mid) );
                    s2 = sum( w(mid:i2) );
                    m1 = min(s1,s2);
                    m2 = max(s1,s2);
                    if m1>0.2*m2
                        if cur_ms1pos<mid
                            IX = i1:mid;
                        else
                            IX = mid:i2;
                        end;
                        break;
                    end;
                end;
            end;
        end;
        %}
        break;
    end;
end;
if 0==bflag
    return;
end;

for j=1:nmz
    c_isointens = c_ref_isointens(:,j);
    [localmax_inten(j),x] = max(c_isointens(IX));
    localmax_rt(j) = isorts(IX(x));
end;