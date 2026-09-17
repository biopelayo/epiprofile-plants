function HH2B_AT_Snapshot(cur_outpath)
%%
% H2B AT Snapshot — aggregates all AT H2B modules
% Repaired: 2026-04-16, Sistema Pelamovic
% - Fixed ref HH2B_03u_1_100 (was missing, now uses HH2BMo_03u_1_100)
% - Replaced generic HH2B_AT_disc with actual module names
% - Added AT08 (LARYNKK, adapted from human LAHYNKR)

% H2B.6 (O23629, 150 aa, without M)
H2B = 'APRAEKKPAEKKPAAEKPVEEKSKAEKAPAEKKPKAGKKLPKEAGAGGDKKKKMKKKSVETYKIYIFKVLKQVHPDIGISSKAMGIMNSFINDIFEKLASESSKKLARYNKKPTITSREIQTAVRLVLPGELAKHAVSEGTKAVTKFTSS';

% get_pos_modi
poses = [];
modis = {};

% H2B conserved unmod panel (EIQTAVR, EVQTAVR, etc.)
out_filename = 'HH2BMo_03u_1_100';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2B AT shared: EIQTAVR (pos 94-100)
out_filename = 'HH2B_AT01_shared_94_100';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2B AT shared: YNKKPTITSR (pos 108-117)
out_filename = 'HH2B_AT02_shared_108_117';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2B AT shared: LVLPGELAKHAVSEGTK (pos 118-134)
out_filename = 'HH2B_AT03_shared_118_134';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2B AT discriminating peptides (variant-specific N-term regions)
out_filename = 'HH2B_AT04_67_disc';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_AT05_14_disc';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_AT06_10_disc';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_AT07_311_disc';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2B AT shared: LARYNKK (pos 80-86, adapted from human LAHYNKR)
out_filename = 'HH2B_AT08_shared_80_86';
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
snapshotfile = fullfile(fileparts(cur_outpath),'HH2B_AT_Snapshot.xls');
fp = fopen(snapshotfile,'w');
if -1==fp
    fprintf('can not open:%s\n',snapshotfile);
    return;
end;
for ino=1:length(H2B)
    fprintf(fp,'%s',H2B(ino));
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
