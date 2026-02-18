function [t3,t4] = find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3,rts4,top1_rt4,inten_sum4,top1_inten_sum4)
%%

if 0==isempty(rts3) && 0==isempty(rts4)
    if top1_rt4>top1_rt3-2 && top1_rt4<top1_rt3+0.5
        t3 = top1_rt3;
        t4 = top1_rt4;
    elseif top1_inten_sum3>top1_inten_sum4
        t3 = top1_rt3;
        id = find( rts4>top1_rt3-2 & rts4<top1_rt3+0.5 );
        if 0==isempty(id)
            [tmp,ix] = max(inten_sum4(id));%#ok
            t4 = rts4(id(ix));
        else
            t4 = top1_rt3-0.5;
        end;
    else
        id = find( rts3>top1_rt4-0.5 & rts3<top1_rt4+2 );
        if 0==isempty(id)
            [tmp,ix] = max(inten_sum3(id));%#ok
            t3 = rts3(id(ix));
        else
            t3 = top1_rt4+0.5;
        end;
        t4 = top1_rt4;
    end;
elseif 0==isempty(rts3)
    t3 = top1_rt3;
    t4 = top1_rt3-0.5;
elseif 0==isempty(rts4)
    t3 = top1_rt4+0.5;
    t4 = top1_rt4;
else
    t3 = 0;
    t4 = 0;
end;