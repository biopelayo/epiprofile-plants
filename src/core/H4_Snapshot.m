function H4_Snapshot(cur_outpath)
%%
% =========================================================================
% H4_Snapshot — Build a position-wise PTM snapshot across H4 protein
% -------------------------------------------------------------------------
% PURPOSE
%   Aggregate biologically relevant PTM labels detected in multiple H4
%   peptide panels into a single position-indexed table for the *entire*
%   H4 sequence. The output is a tab-separated file (with .xls extension
%   for Excel-friendliness) listing, per residue, any PTM labels observed
%   at that absolute protein position across the input panels.
%
% INPUT
%   cur_outpath : folder containing per-panel .mat outputs, each expected
%                 to define 'His' (with 'mod_type') and 'auc' (areas).
%
% OUTPUT (side effect)
%   <parent_of_cur_outpath>/H4_Snapshot.xls — TSV where each line has:
%     - column 1: the single-letter amino acid for that H4 position
%     - subsequent columns: PTM labels (excluding 'pr' and 'ox') detected
%                           for that absolute position (one label per cell)
%
% ASSUMPTIONS
%   - The .mat files follow the naming pattern 'H4_XX_start_end.mat'.
%   - Column 2 of 'auc' is used as presence proxy (auc(state,2) > 0).
%   - 'His.mod_type' uses the "index,label;" grammar where index 0 refers
%     to peptide N-terminus, whereas 1..L refer to residue indices within
%     the peptide window. Protein-absolute positions are obtained by adding
%     (start-1) parsed from the file name.
%   - Labels 'pr' (propionyl) and 'ox' (oxidation) are considered chemical
%     artifacts or derivatization marks and thus are *ignored* in the
%     snapshot to keep biological PTMs only.
% =========================================================================

% ----------------------------
% H4 canonical sequence used for indexing in this analysis context.
% Each printed line in the snapshot will begin with this amino acid.
% ----------------------------
H4 = 'SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG';

% ----------------------------
% Accumulators for absolute positions and PTM labels
% 'poses' holds protein-absolute indices (1-based w.r.t. H4 above)
% 'modis' holds the corresponding PTM label strings (e.g., 'ac','me1',...)
% ----------------------------
poses = [];
modis = {};

% ----------------------------
% Pull PTMs from selected peptide panels (extend as needed)
% NOTE: Each get_pos_modi call will append entries to poses/modis if the
%       panel's .mat exists and reports non-zero auc(:,2).
% ----------------------------
out_filename = 'H4_01_4_17';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H4_02_20_23';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

out_filename = 'H4_04_40_45';
[poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis);

% ----------------------------
% Deduplicate (position,label) pairs
% Strategy:
%   1) sort by absolute position
%   2) mark duplicates where both position and label match
% ----------------------------
[poses,I] = sort(poses,'ascend');
modis = modis(I);

flag = repmat(1,[length(poses),1]); % 1=keep, 0=drop duplicate
for ino=1:length(poses)-1
    if 0==flag(ino)
        continue;
    end
    p_i = poses(ino);
    m_i = modis{ino};
    for jno=ino+1:length(poses)
        p_j = poses(jno);
        m_j = modis{jno};
        if 0~=p_j-p_i
            break;              % positions differ; no more duplicates ahead
        end
        if 1==strcmp(m_j,m_i)   % same position AND same label -> drop later
            flag(jno) = 0;
        end
    end
end
II = find(flag==1);
poses = poses(II);
modis = modis(II);

% ----------------------------
% Write snapshot to disk (TSV but .xls extension for convenience)
% One line per residue of H4; if PTMs exist at that position, print each
% PTM label in its own tab-separated cell after the amino acid.
% ----------------------------
snapshotfile = fullfile(fileparts(cur_outpath),'H4_Snapshot.xls');
fp = fopen(snapshotfile,'w');
if -1==fp
    fprintf('can not open:%s\n',snapshotfile);
    return;
end
for ino=1:length(H4)
    fprintf(fp,'%s',H4(ino));         % print the amino acid letter
    III = find(poses==ino);           % all PTM entries for this position
    if 0==isempty(III)
        for jno=1:length(III)
            fprintf(fp,'\t%s',modis{III(jno)}); % each label in its own cell
        end
    end
    fprintf(fp,'\n');                  % next residue/line
end
fclose(fp);

end % H4_Snapshot


% =========================================================================
% get_pos_modi — read one panel's .mat, collect absolute-position PTMs
% -------------------------------------------------------------------------
% For each detected state (auc(state,2) > 0), parse His.mod_type{state}
% entries of the form ';idx,label;idx,label;...' and:
%   - obtain the peptide-relative index 'idx' (0 = peptide N-terminus)
%   - map to protein-absolute index by adding (start-1) parsed from the
%     file name 'H4_XX_start_end'
%   - keep the label if it is NOT 'pr' or 'ox'
%   - append to poses/modis accumulators
% =========================================================================
function [poses,modis] = get_pos_modi(cur_outpath,out_filename,poses,modis)
%%
out_file = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file,'file')
    % ------------------------------------------------------------
    % Parse peptide window start from file name
    % Expected pattern: 'H4_XX_start_end' (underscores at fixed places)
    % start_pt is (start-1) to convert peptide-relative indices to
    % protein-absolute indices.
    % ------------------------------------------------------------
    p = strfind(out_filename,'_');
    start_pt = str2double(out_filename(p(2)+1:p(3)-1))-1; % subtract 1

    % ------------------------------------------------------------
    % Load per-panel results; expected to contain:
    %   - His.mod_type : cell array of PTM strings ('i,label;...;')
    %   - auc          : areas matrix (presence when auc(state,2) > 0)
    % ------------------------------------------------------------
    load(out_file); % loads variables from the .mat into workspace (His, auc)

    nlen = length(His.mod_type);
    for ino=1:nlen
        if 0==auc(ino,2) % area==0 -> considered absent -> skip
            continue;
        end
        % --------------------------------------------------------
        % Normalize the string so we can safely locate ';' and ','
        % We prepend a leading ';' to simplify parsing the first token.
        % --------------------------------------------------------
        cur_mod_type = [';',His.mod_type{ino}];
        p1 = strfind(cur_mod_type,';'); % token boundaries
        p2 = strfind(cur_mod_type,','); % index/label separators

        % Iterate over all "index,label" pairs in this state string
        for jno=1:length(p2)
            % Extract peptide-relative index (as number)
            cur_pose = str2double( cur_mod_type(p1(jno)+1:p2(jno)-1) );
            % Extract label between comma and next semicolon
            cur_modi = cur_mod_type(p2(jno)+1:p1(jno+1)-1);

            % Skip derivatization/artefact labels ('pr') and generic ox ('ox')
            if 0==ismember(cur_modi,{'pr','ox'})
                % Map to protein-absolute index: add start offset
                poses(end+1) = cur_pose+start_pt; %#ok
                modis{end+1,1} = cur_modi;        %#ok
            end
        end
    end
end

end % get_pos_modi
