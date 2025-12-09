function H3_Snapshot(cur_outpath)
%%
% H3_Snapshot
% -------------------------------------------------------------------------
% Purpose:
%   Consolidate site-level PTM evidence across multiple H3 panels into a
%   single tab-delimited snapshot file ('H3_Snapshot.xls'). For each panel:
%     - Load <panel>.mat (expects 'His' and 'auc')
%     - For each quantified row (auc(:,2) ~= 0), parse the compact mod_type
%       string "idx,mod;idx,mod;..." to collect (position, modification)
%       pairs, ignoring 'pr' and 'ox'.
%     - Convert panel-relative positions to absolute coordinates in H3
%       using the <start> index derived from the panel filename.
%   Finally:
%     - Sort by position, remove duplicate (pos,mod) entries
%     - Emit one line per residue of full-length H3 with any mods found
%       listed as additional tab-separated columns.
%
% Inputs:
%   cur_outpath : directory containing the panel .mat files
%
% Outputs:
%   <parent_of_cur_outpath>/H3_Snapshot.xls
%
% Notes:
%   - The file extension ".xls" is a plain tab-delimited text file for
%     convenience; it is not a native Excel binary format.
%   - Panel filename pattern must be H3_<block>_<start>_<end>
%     so that <start> can be extracted to align relative positions.
% -------------------------------------------------------------------------

% ---------------------------- H3 full sequence ---------------------------
% Full H3 protein sequence used to emit one line per residue.
H3 = 'ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA';

% -------------------------- aggregate containers -------------------------
poses = [];  % absolute positions accross panels
modis = {};  % corresponding modification labels

% ------------------------------ panel list -------------------------------
% For each panel, we call get_pos_modi which appends positions and mods.
out_filename = 'H3_01_3_8';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_02_9_17';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_02a_9_17';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_02b_9_17';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_03_18_26';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_04_27_40';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_04a_27_40';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_06_54_63';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_07_73_83';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H3_08_117_128';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% --------------------------- sort and unique -----------------------------
% Sort by absolute position
[poses,I] = sort(poses,'ascend');
modis = modis(I);

% Remove duplicates for identical positions (same mod label)
flag = repmat(1,[length(poses),1]);
for ino=1:length(poses)-1
    if 0==flag(ino), continue; end
    p_i = poses(ino);
    m_i = modis{ino};
    for jno=ino+1:length(poses)
        p_j = poses(jno);
        m_j = modis{jno};
        if 0~=p_j-p_i
            % Positions differ; no further duplicates at this position
            break;
        end
        if 1==strcmp(m_j,m_i)
            % Exact duplicate (position + same mod label)
            flag(jno) = 0;
        end
    end
end
II    = find(flag==1);
poses = poses(II);
modis = modis(II);

% ------------------------------- write out --------------------------------
snapshotfile = fullfile(fileparts(cur_outpath),'H3_Snapshot.xls');
fp = fopen(snapshotfile,'w');
if -1==fp
    fprintf('can not open:%s\n',snapshotfile);
    return;
end

% Emit one line per residue (AA + any mods at that absolute index)
for ino=1:length(H3)
    fprintf(fp,'%s',H3(ino));
    III = find(poses==ino);
    if 0==isempty(III)
        for jno=1:length(III)
            fprintf(fp,'\t%s',modis{III(jno)});
        end
    end
    fprintf(fp,'\n');
end
fclose(fp);

end % H3_Snapshot

% =========================================================================
function [poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis)
%%
% get_pos_modi
% -------------------------------------------------------------------------
% Purpose:
%   For a given panel <out_filename>, append absolute (pos, mod) entries
%   to the running lists, derived from the panel's .mat results.
%
% Inputs:
%   cur_outpath  : directory containing <out_filename>.mat
%   out_filename : basename of the panel (H3_<block>_<start>_<end>)
%   poses, modis : current aggregates to update
%
% Outputs:
%   poses, modis : updated aggregates with entries from this panel
%
% Requirements in the .mat:
%   - 'His' with cell array 'mod_type'
%   - 'auc' where auc(ino,2) is the area/intensity used to filter rows
%
% Notes:
%   - Ignores 'pr' and 'ox' labels (derivatization / oxidation).
%   - Converts panel-relative indices to absolute H3 coordinates using the
%     <start> extracted from out_filename (minus 1).
% -------------------------------------------------------------------------

out_file = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file,'file')
    % Extract <start> index from the filename pattern H3_<block>_<start>_<end>
    p = strfind(out_filename,'_');
    start_pt = str2double(out_filename(p(2)+1:p(3)-1))-1; % subtract 1 (core convention)

    % Load required variables: His, auc
    load(out_file); % expects variables 'His' and 'auc'

    nlen = length(His.mod_type);
    for ino=1:nlen
        % Skip rows with zero area (no quantitative evidence)
        if 0==auc(ino,2)
            continue;
        end

        % Prepend a semicolon to simplify parsing "idx,mod;" segments
        cur_mod_type = [';',His.mod_type{ino}];
        p1 = strfind(cur_mod_type,';'); % semicolon positions
        p2 = strfind(cur_mod_type,','); % comma positions

        % Iterate over segments (each "idx,mod;" pair)
        for jno=1:length(p2)
            % Extract "idx" and "mod" between delimiters
            cur_pose = str2double(cur_mod_type(p1(jno)+1 : p2(jno)-1));
            cur_modi = cur_mod_type(p2(jno)+1 : p1(jno+1)-1);

            % Ignore derivatization ('pr') and oxidation ('ox')
            if 0==ismember(cur_modi,{'pr','ox'})
                poses(end+1)     = cur_pose + start_pt; %#ok<AGROW>
                modis{end+1,1}   = cur_modi;            %#ok<AGROW>
            end
        end
    end
end

end % get_pos_modi
