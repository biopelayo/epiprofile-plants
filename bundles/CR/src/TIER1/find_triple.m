function xt = find_triple(rts1,top1_rt1,rts2,rts3,inten_sum1,inten_sum2,inten_sum3)
%%

if 1==isempty(top1_rt1)
    xt = zeros([3,3]);
    return;
end;

% [tmp_sum,ix] = sort(inten_sum1,'descend');
% tmp_rts = rts1(ix);
% if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<12
%     top1_rt1 = min([tmp_rts(2),tmp_rts(1)]);
% end;

xt = find_triple_once(rts1,top1_rt1,rts2,rts3,inten_sum1,inten_sum2,inten_sum3);

function xt = find_triple_once(rts1,top1_rt1,rts2,rts3,inten_sum1,inten_sum2,inten_sum3)
%%

xt = zeros([3,3]);

% 1st col
xt(1,1) = top1_rt1;

id = find( rts2>top1_rt1-1.2 & rts2<top1_rt1+0.7 );
if 0==isempty(id)
    [tmp,ix] = max(inten_sum2(id));%#ok
    xt(2,1) = rts2(id(ix));
    %top1_rt1 = min([xt(1,1) xt(2,1)]);% smaller
end;

% 2nd col
id = find(rts1>=top1_rt1+0.7);
if 0==isempty(id)
    %[tmp,ix] = max(inten_sum1(id));
    %xt(1,2) = rts1(id(ix));
    [tmp_sum,ix] = sort(inten_sum1(id),'descend');
    tmp_rts = rts1(id(ix));
    if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/10
        new_tmp_rt = min([tmp_rts(2),tmp_rts(1)]);
        if 0==isempty(find( abs(rts2-new_tmp_rt)<0.7 ))%#ok
            xt(1,2) = new_tmp_rt;
        else
            xt(1,2) = tmp_rts(1);
        end;
    else
        xt(1,2) = tmp_rts(1);
    end;
end;

if xt(1,2)>0
    id = find(rts2>=top1_rt1+0.7 & rts2<=xt(1,2)+1.2);
else
    id = find(rts2>=top1_rt1+0.7 & rts2<=top1_rt1+1.9);
end;
if 0==isempty(id)
    [tmp_sum,ix] = sort(inten_sum2(id),'descend');
    tmp_rts = rts2(id(ix));
    if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/10 && abs(tmp_rts(2)-tmp_rts(1))<1.3
        xt(2,2) = min([tmp_rts(2),tmp_rts(1)]);
        xt(1,3) = max([tmp_rts(2),tmp_rts(1)]);
    else
        xt(2,2) = tmp_rts(1);
        xt(1,3) = tmp_rts(1)+0.02;
    end;
end;

if xt(1,2)>0
    id = find(rts3>=top1_rt1+0.7 & rts3<=xt(1,2)+1.2);
else
    id = find(rts3>=top1_rt1+0.7 & rts3<=top1_rt1+1.9);
end;
if 0==isempty(id)
    [tmp,ix] = max(inten_sum3(id));%#ok
    xt(3,2) = rts3(id(ix));
end;

id = find(abs(rts1-(xt(2,2)+xt(3,2))/2)<0.7);
if 0==isempty(id)
    [tmp,ix] = min(abs(rts1(id)-(xt(2,2)+xt(3,2))/2));%#ok
    xt(1,2) = rts1(id(ix));
end;

if (0==xt(1,2) && abs(xt(2,2)-xt(3,2))>0.7) || 0==xt(2,2) || 0==xt(3,2)
    return;
end;

% 3rd col
id = find(rts2>=xt(1,3)+1);
if 0==isempty(id)
    [tmp,ix] = max(inten_sum2(id));%#ok
    xt(2,3) = rts2(id(ix));
end;

id = find(rts3>=xt(1,3)+1 & rts3<=xt(2,3)+1.2);
if 0==isempty(id)
    %[tmp,ix] = max(inten_sum3(id));%#ok
    i = find(abs(rts3(id)-xt(2,3))<=0.5);
    if 1==isempty(i)
        xt(3,3) = xt(2,3);
    else
        [tmp,ix] = max(inten_sum3(id(i)));%#ok
        xt(3,3) = rts3(id(i(ix)));
    end;
end;