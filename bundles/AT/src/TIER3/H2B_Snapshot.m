function H2B_Snapshot(cur_outpath)
%%
% H2B_Snapshot — Aggregates PTM detections across plant H2B variant modules.
%
%  TIER: T3 (plant-specific)
%  Currently all AT H2B modules quantify unmod only. This snapshot is
%  structural: it maps detected peptides onto the canonical H2B sequence
%  and will capture any unexpected PTMs if present.
%  Modules: HH2B_01oA..HH2B_13oA (H2B.1), HH2B_14oA..HH2B_26oA (H2B.11)

% H2B canonical (AT_H2B_1 / HTB1, At1g07790)
H2B = 'MAAAEKKPAEKTAAERPVEENKAAPAEKAAAAKKPESGDEAGDKNVETYKIYIFKVLKQVHPDIGISSKAMGIMNSFINDIFEKLAQESSKYNKKPTITSREIQTAVRLVLPGELAKHAVSEGTKAVTKFTSS';

% get_pos_modi
poses = [];
modis = {};

% H2B.1 modules
out_filename = 'HH2B_01oA_8_12';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_02oA_14_24';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_03oA_29_33';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_04oA_45_49';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_05oA_57_62';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_06oA_63_67';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_07oA_71_81';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_08oA_82_96';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_09oA_97_103';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_10oA_110_116';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_11oA_117_123';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_12oA_124_132';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_13oA_133_140';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% H2B.11 modules
out_filename = 'HH2B_14oA_8_12';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_15oA_14_24';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_16oA_29_33';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_17oA_45_49';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_18oA_57_62';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_19oA_63_67';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_20oA_71_81';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_21oA_82_96';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_22oA_97_103';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_23oA_110_116';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_24oA_117_123';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_25oA_124_132';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'HH2B_26oA_133_140';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% get unique
if isempty(poses)
    return;
end;
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
snapshotfile = fullfile(fileparts(cur_outpath),'H2B_Snapshot.xls');
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
    start_pt = str2double(out_filename(p(2)+1:p(3)-1))-1;% subtract 1

    load(out_file);% His, auc
    nlen = length(His.mod_type);
    for ino=1:nlen
        if 0==auc(ino,2)% area
            continue;
        end;
        cur_mod_type = [';',His.mod_type{ino}];
        p1 = strfind(cur_mod_type,';');
        p2 = strfind(cur_mod_type,',');
        for jno=1:length(p2)
            cur_pose = str2double( cur_mod_type(p1(jno)+1:p2(jno)-1) );
            cur_modi = cur_mod_type(p2(jno)+1:p1(jno+1)-1);
            if 0==ismember(cur_modi,{'pr','ox','pic','acD3'})
                poses(end+1) = cur_pose+start_pt;%#ok
                modis{end+1,1} = cur_modi;%#ok
            end;
        end;
    end;
end;
