function [isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,cur_ch,ptol,nC13,PosRange)
%%
% MS1_index (scan, retention time, peak number, baseline), MS1_peaks (m/z, intensity)

nmz = length(c_ref_isomzs);
num_MS1 = size(MS1_index,1);
index = [1;MS1_index(1:num_MS1,3)];

isorts = MS1_index(1:num_MS1,2);
c_ref_isointens = zeros([num_MS1,nmz]);

if 1==nC13
    for ino = PosRange
        IX = index(ino):index(ino+1)-1;
        mz = MS1_peaks(IX,1);
        inten = MS1_peaks(IX,2);

        tmp_intens = zeros([1,nmz]);
        for jno=1:nmz
            c_mz = c_ref_isomzs(jno);
            c_ptol = min([ptol*c_mz*1e-6,0.3]);
            left = c_mz - c_ptol;
            right = c_mz + c_ptol;
            pos = find( mz>=left & mz<=right );
            if 0==isempty(pos)
                tmp_intens(jno) = max(inten(pos));
            end;
        end;
        c_ref_isointens(ino,1:nmz) = tmp_intens;
    end;
else
    thr = 5;
    IPV = GetIPV();
    M = c_ref_isomzs(2)*cur_ch-cur_ch*1.007276;
    tIPV = IPV(floor(M),1:5);
    t1 = tIPV(2);

    for ino = PosRange
        IX = index(ino):index(ino+1)-1;
        mz = MS1_peaks(IX,1);
        inten = MS1_peaks(IX,2);

        flags = zeros([1,length(mz)]);
        for jno=1:nmz
            c_mz = c_ref_isomzs(jno);
            c_ptol = min([ptol*c_mz*1e-6,0.3]);
            left = c_mz - c_ptol;
            right = c_mz + c_ptol;
            pos = find( mz>=left & mz<=right );
            if 0==isempty(pos)
                flags(pos) = jno;
            end;
        end;

        tmp_intens = zeros([1,nmz]);
        x1 = find(flags==1);
        x2 = find(flags==2);
        x3 = find(flags==3);
        x4 = find(flags==4);
        if 0==isempty(x2) && 0==isempty(x3)
            array23 = zeros([length(x2)*length(x3),1]);
            for pno = 1:length(x2)
                for qno = 1:length(x3)
                    array23((pno-1)*length(x3)+qno,1) = abs( log( (inten(x3(qno))/inten(x2(pno)))/t1 ) );
                end;
            end;
            [minlogR23,idx] = min(array23);
            qno = mod(idx,length(x3));
            if 0==qno
                qno = length(x3);
            end;
            pno = (idx-qno)/length(x3)+1;

            if (1==isempty(x1) || min(inten(x1))<inten(x2(pno))) && minlogR23<=log(thr)
                if 0==isempty(x1)
                    tmp_intens(1) = min(inten(x1));
                end;
                tmp_intens(2) = inten(x2(pno));
                tmp_intens(3) = inten(x3(qno));
                if 0==isempty(x4)
                    array34 = zeros([length(x4),1]);
                    for rno = 1:length(x4)
                        array34(rno,1) = abs( log( (inten(x4(rno))/inten(x3(qno)))/t1 ) );
                    end;
                    [minlogR34,rno] = min(array34);
                    if 2==nC13 || minlogR34<=log(thr)
                        tmp_intens(4) = inten(x4(rno));
                    end;
                end;
            end;
        end;
        c_ref_isointens(ino,1:nmz) = tmp_intens;
    end;
end;