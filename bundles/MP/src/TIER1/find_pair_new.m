function [t3,t4] = find_pair_new(top1_rt3,rts4,top1_rt4,inten_sum4,rt3_up)
%%

if 1==rt3_up
    if 0==isempty(top1_rt3) && 0==isempty(top1_rt4)
        t3 = top1_rt3;
        id = find( rts4>top1_rt3-2 & rts4<top1_rt3+0.5 );
        if 0==isempty(id)
            [tmp,ix] = max(inten_sum4(id));%#ok
            t4 = rts4(id(ix));
        else
            t4 = top1_rt3-0.5;
        end;
    elseif 0==isempty(top1_rt3)
        t3 = top1_rt3;
        t4 = top1_rt3-0.5;
    elseif 0==isempty(top1_rt4)
        t3 = top1_rt4+0.5;
        t4 = top1_rt4;
    else
        t3 = 0;
        t4 = 0;
    end;
else
    if 0==isempty(top1_rt3) && 0==isempty(top1_rt4)
        t3 = top1_rt3;
        id = find( rts4>top1_rt3-0.5 & rts4<top1_rt3+2 );
        if 0==isempty(id)
            [tmp,ix] = max(inten_sum4(id));%#ok
            t4 = rts4(id(ix));
        else
            t4 = top1_rt3+0.5;
        end;
    elseif 0==isempty(top1_rt3)
        t3 = top1_rt3;
        t4 = top1_rt3+0.5;
    elseif 0==isempty(top1_rt4)
        t3 = top1_rt4-0.5;
        t4 = top1_rt4;
    else
        t3 = 0;
        t4 = 0;
    end;
end;