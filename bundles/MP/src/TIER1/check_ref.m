function [rt_ref,ndebug] = check_ref(raw_path,new_seq,rt_ref,ndebug)

mat_file = fullfile(raw_path,'histone_layouts','0_ref_info.mat');
if (0==ndebug || -1==ndebug) && 0~=exist(mat_file,'file')
    load(mat_file);
    new_godel = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    ix = find(AllUnHis.seq_godel==new_godel);
    if 1==isempty(ix)
        return;
    elseif 1==length(ix)
        rt_ref = AllUnHis.rt_ref(ix);
        ndebug = -1;
    else
        nflag = 0;
        for i=1:length(ix)
            cur_seq = [AllUnHis.pep_seq{ix(i),1},AllUnHis.mod_type{ix(i),1}];
            if 1==strcmp(cur_seq,new_seq)
                nflag = 1;
                break;
            end;
        end;
        if 1==nflag
            rt_ref = AllUnHis.rt_ref(ix(i));
            ndebug = -1;
        end;
    end;
end;