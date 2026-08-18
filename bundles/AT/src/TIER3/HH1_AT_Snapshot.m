function HH1_AT_Snapshot(cur_outpath)
%%

% H1.1 (P26568, 274 aa, without M)
H1 = 'SEVEIENAATIEGNTAADAPVTDAAVEKKPAAKGRKTKNVKEVKEKKTVAAAPKKRTVSSHPTYEEMIKDAIVTLKERTGSSQYAIQKFIEEKRKELPPTFRKLLLLNLKRLVASGKLVKVKASFKLPSASAKASSPKAAAEKSAPAKKKPATVAVTKAKRKVAAASKAKKTIAVKPKTAAAKKVTAKAKAKPVPRATAAATKRKAVDAKPKAKARPAKAAKTAKVTSPAKKAVAATKKVATVATKKKTPVKKVVKPKTVKSPAKRASSRVKK';

% get_pos_modi
poses = [];
modis = {};

% H1.1: KLLLLNLKR (pos 104-112)
out_filename = 'HH1_AT01_H11_104_112';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H1.1: TGSSQYAIQKFIEEKR (pos 80-95)
out_filename = 'HH1_AT02_H11_80_95';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H1.2: KLLLVNLKR (pos 104-112)
out_filename = 'HH1_AT03_H12_104_112';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H1.2: TGSSQYAIQKFIEEKHKSLPPTFR (pos 80-103)
out_filename = 'HH1_AT04_H12_80_103';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% get unique
[poses,I] = sort(poses,'ascend');
modis = modis(I);

flag = repmat(1,[length(poses),1]);
for ino=1:length(poses)-1
    if 0==flag(ino)
        continue;
    end;
    p_i = poses(ino);
    m_i = modis{ino};
    for jno=ino+1:length(poses)
        p_j = poses(jno);
        m_j = modis{jno};
        if 0~=p_j-p_i
            break;
        end;
        if 1==strcmp(m_j,m_i)
            flag(jno) = 0;
        end;
    end;
end;
II = find(flag==1);
poses = poses(II);
modis = modis(II);

% output
snapshotfile = fullfile(fileparts(cur_outpath),'HH1_AT_Snapshot.xls');
fp = fopen(snapshotfile,'w');
if -1==fp
    fprintf('can not open:%s\n',snapshotfile);
    return;
end;
for ino=1:length(H1)
    fprintf(fp,'%s',H1(ino));
    III = find(poses==ino);
    if 0==isempty(III)
        for jno=1:length(III)
            fprintf(fp,'\t%s',modis{III(jno)});
        end;
    end;
    fprintf(fp,'\n');
end;
fclose(fp);

function [poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis)
%%

out_file = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file,'file')
    p = strfind(out_filename,'_');
    if length(p)>=3
        start_pt = str2double(out_filename(p(end-1)+1:p(end)-1))-1;
    else
        start_pt = 0;
    end;

    load(out_file);% His, auc
    nlen = length(His.mod_type);
    for ino=1:nlen
        if 0==auc(ino,2)
            continue;
        end;
        cur_mod_type = [';',His.mod_type{ino}];
        p1 = strfind(cur_mod_type,';');
        p2 = strfind(cur_mod_type,',');
        for jno=1:length(p2)
            cur_pose = str2double( cur_mod_type(p1(jno)+1:p2(jno)-1) );
            cur_modi = cur_mod_type(p2(jno)+1:p1(jno+1)-1);
            if 0==ismember(cur_modi,{'pr','ox'})
                poses(end+1) = cur_pose+start_pt;%#ok
                modis{end+1,1} = cur_modi;%#ok
            end;
        end;
    end;
end;
