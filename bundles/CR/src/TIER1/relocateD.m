function His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His)
% debug

if ptol==100
    ptol = 10;
end;
delta = 1;

for hno = 2:length(His.rt_ref)
    t1 = His.rt_ref(hno)-delta;
    t2 = His.rt_ref(hno)+delta;
    [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);
    
    if 1==isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end;
end;