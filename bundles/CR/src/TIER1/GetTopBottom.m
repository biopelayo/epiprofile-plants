function [nt,nb,top1_idx,inten_sum] = GetTopBottom(w)
%% JudgeLocalmaxmin

% w = smooth(w,3);
w(find(w<0.015*max(w))) = 0;%#ok

ix = find(w==0);
if length(w)==length(ix)
    nt = [];
    nb = [];
    top1_idx = [];
    inten_sum = [];
    return;
end;

x1 = ix(1:end-1);
x2 = ix(2:end);
dx = x2-x1;
idx = find(dx>1);

if 1==isempty(idx)
    nb = [1;length(w)];
else
    nb = [ix(idx);length(w)];
end;
% split non-zero peaks
for i=1:length(nb)-1
    i1 = nb(i);
    i2 = nb(i+1);
    IX = i1:i2;

    c_max = max(w(IX));
    [t,b] = JudgeLocalmaxmin(w(IX));%#ok
    if length(b)>2
        for j=2:length(b)-1
            nmid = IX(b(j));
            if w(nmid)<0.2*c_max
                s1 = sum( w(i1:nmid) );
                s2 = sum( w(nmid:i2) );
                m1 = min(s1,s2);
                m2 = max(s1,s2);
                if m1>0.05*m2
                    nb = [nb;nmid];%#ok
                    break;
                end;
            end;
        end;
    end;
end;
nb = sort(nb);
% split non-zero peaks
for i=1:length(nb)-1
    i1 = nb(i);
    i2 = nb(i+1);
    IX = i1:i2;

    c_max = max(w(IX));
    [t,b] = JudgeLocalmaxmin(w(IX));%#ok
    if length(b)>2
        for j=2:length(b)-1
            nmid = IX(b(j));
            if w(nmid)<0.2*c_max
                s1 = sum( w(i1:nmid) );
                s2 = sum( w(nmid:i2) );
                m1 = min(s1,s2);
                m2 = max(s1,s2);
                if m1>0.05*m2
                    nb = [nb;nmid];%#ok
                    break;
                end;
            end;
        end;
    end;
end;
nb = sort(nb);

nt = zeros([length(nb)-1,1]);
inten_sum = zeros([length(nb)-1,1]);
for i=1:length(nb)-1
    i1 = nb(i);
    i2 = nb(i+1);
    while i1<=length(w) && 0==w(i1)
        i1 = i1+1;
    end;
    if i1>nb(i)
        i1 = i1-1;
    end;
    while i2>=1 && 0==w(i2)
        i2 = i2-1;
    end;
    if i2<nb(i+1)
        i2 = i2+1;
    end;
    IX = i1:i2;
    %[tmp,id] = max(w(IX));%#ok
    %nt(i) = IX(id);
    nt(i) = floor(IX*w(IX)/sum(w(IX)));

    nb(i) = i1;
    nb(i+1) = i2;
    inten_sum(i) = sum(w(nb(i):nb(i+1)));
end;
[tmp,top1_idx] = max(inten_sum);%#ok