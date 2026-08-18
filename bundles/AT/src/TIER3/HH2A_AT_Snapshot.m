function HH2A_AT_Snapshot(cur_outpath)
%%
% H2A AT Snapshot — aggregates all AT H2A modules
% Repaired: 2026-04-16, Cleaned: 2026-05-24, Sistema Pelamovic
% - Replaced generic/human refs with actual AT module out_filenames
% - [Retired 2026-05-24] AT09 (KGKYAER) - peptido fantasma
%   (Tribunal de Codigo Capa B + Art. 16 FN-1)

% H2A canonical (A0A178WDN2, 132 aa, without M)
H2A = 'AGRGKTLGSGSAKKATTRSSSKAGLQFPVGRIHRFLKKGKYAERVGAGAPVYLAAVLEYLAAEVLELAGNAARRDNKKTRIHRIHIQLAVRNDEELSKLLGDVTIANGGVMPNIHNLLLPKKTGASKPSAEDD';

% get_pos_modi
poses = [];
modis = {};

% H2A canonical shared: SSKAGLQFPVGR (pos 20-31)
out_filename = 'HH2A_AT01_can_20_31';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2A canonical 1: GKTLGSGSAKKATTR (pos 5-19)
out_filename = 'HH2A_AT02_can1_5_19';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2A.Z: AGIQFPVGR (pos 33-41)
out_filename = 'HH2A_AT03_Z_33_41';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2A.X: GKPKATKSVSR (pos 15-25)
out_filename = 'HH2A_AT04_X_15_25';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2A.Xa: FLKSGKYAER (pos 41-50)
out_filename = 'HH2A_AT05_Xa_41_50';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2A.W12: SGGGPKKKPVSR (pos 16-27)
out_filename = 'HH2A_AT06_W12_16_27';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2A.W12: SVKSGLQFPVGR (pos 28-39)
out_filename = 'HH2A_AT07_W12_28_39';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2A.W6: MESTGKVKKAFGGR (pos 1-14, N-term)
out_filename = 'HH2A_AT08_W6_1_14';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% [Retirado 2026-05-24] HH2A_AT09 (KGKYAER) - peptido fantasma.

% H2A conserved: DNKKTR (pos 72-77, shared human/AT)
out_filename = 'HH2A_06m1_72_77';
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
snapshotfile = fullfile(fileparts(cur_outpath),'HH2A_AT_Snapshot.xls');
fp = fopen(snapshotfile,'w');
if -1==fp
    fprintf('can not open:%s\n',snapshotfile);
    return;
end;
for ino=1:length(H2A)
    fprintf(fp,'%s',H2A(ino));
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
