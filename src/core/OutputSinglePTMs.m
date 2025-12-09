function OutputSinglePTMs(layout_path,raw_names)
%% ========================================================================
%  OutputSinglePTMs
%  ------------------------------------------------------------------------
%  PURPOSE (high level, plain English, non-proteomics audience):
%    This function scans a set of per-sample result files produced by a
%    previous histone-proteomics workflow (EpiProfile-like). For each sample,
%    it reads peptide-level results and aggregates quantitative signals for a
%    list of "single PTMs" (post-translational modifications, e.g. H3K4me1).
%
%    In simpler words: given a folder with many small .mat files (one per
%    peptide window/region), it sums up the “ratio” signal for specific
%    histone marks (targets) across all relevant peptides, and then outputs
%    a final summary table (one row per PTM, one column per sample).
%
%  WHAT THE FUNCTION OUTPUTS:
%    1) <layout_path>/histone_ratios_single_PTMs.mat
%         - contains:
%             targets : cell array of PTM names (e.g., 'H3K4me1', ...)
%             sratios : matrix [numTargets x numSamples] with the aggregated
%                       "ratio" per target and sample.
%    2) <layout_path>/histone_ratios_single_PTMs.xls
%         - a tab-separated text file (".xls" extension for convenience)
%           with header (sample IDs) and each row = PTM target.
%
%  KEY ASSUMPTIONS ABOUT THE .mat FILES BEING READ:
%    - For each sample i, there is a directory:
%         <layout_path>/<twoDigitIndex>_<raw_name>/detail/
%      Example: layout_path = '/project/results'
%               raw_names{1} = 'SampleA'
%               Folder: '/project/results/01_SampleA/detail/'
%
%    - Inside each /detail/ directory there are many .mat files. Each .mat
%      file corresponds to a peptide “window/region” and must contain at least:
%         * His : a struct with fields used here:
%                 - pep_mz     : [nPeptides x 1] precursor m/z per peptide
%                 - mod_short  : {nPeptides x 1} short labels of mods
%                                (e.g. 'K4me1.R8-K14' or similar)
%         * auc : [nPeptides x 3] numeric matrix with three columns:
%                 (1) RT (retention time), (2) Area (peak area), (3) Ratio
%                 NOTE: this function only aggregates the 3rd column (Ratio).
%
%  WHAT "RATIO" MEANS HERE (conceptual):
%    - In EpiProfile-style histone workflows, each peptide may appear in
%      different modified forms (e.g., unmodified, mono-, di-, tri-methylated).
%      A "ratio" typically represents the fraction of a specific modified
%      form relative to the sum of all forms of that peptide (normalization).
%      This function *sums those ratios across relevant peptides* to get a
%      target-level summary per sample. (This is a practical heuristic used
%      in many pipelines; interpret with caution.)
%
%  HOW TARGET MATCHING WORKS (simple rule-based string matching):
%    - The 'targets' list contains PTM names like 'H3K4me1'.
%    - For each target:
%        * If it begins with 'H31' or 'H33', the code limits the search to
%          peptides from specific sequence windows (e.g. 'H3_27_40' vs 'H33_27_40')
%          to distinguish histone variants H3.1 vs H3.3.
%        * Otherwise, it searches across all peptides.
%        * It then matches modification substrings (e.g., 'K4me1', 'K9ac')
%          inside each peptide label and sums their "ratio".
%
%  INPUTS:
%    layout_path : (char) root directory containing numbered sample folders.
%    raw_names   : (cellstr) sample names, e.g., {'S1','S2','S3'}
%
%  OUTPUT FILES:
%    - histone_ratios_single_PTMs.mat  (MATLAB .mat file)
%    - histone_ratios_single_PTMs.xls  (tab-separated table with .xls suffix)
%
%  IMPORTANT LIMITATIONS:
%    - This function presumes a very specific folder layout and variable
%      names inside each .mat file (His, auc). If your pipeline produces a
%      different structure, adapt accordingly.
%    - PTM detection uses substring matching; if naming conventions differ
%      (e.g., extra characters), adjust the rules accordingly.
%    - The ".xls" is *not* a real Excel binary file; it is a tab-separated
%      text file that Excel can open.
%
%  ------------------------------------------------------------------------
%  AUTHORSHIP / CONTEXT:
%    Original logic adapted from EpiProfile-style outputs. This is a pure
%    aggregation step; it does not perform peak picking nor identification.
%
%  ========================================================================


%% nos, peptides, and info
% -------------------------------------------------------------------------
% Initialize processing for the FIRST sample to discover how many peptides
% ('npeps') we will have when stacking all peptide entries across all the
% .mat files in its /detail/ folder. We only look at the first sample here
% to size arrays consistently for all samples.
i = 1;
cur_rawname = raw_names{i};

% Build a two-digit prefix "01", "02", ... for folder naming.
if i<10
    prefix = ['0',num2str(i)];
else
    prefix = num2str(i);
end;
cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');

% ------------------------------
% Count total number of peptide rows (npeps) contributed by all .mat files
% in the current /detail/ directory. We will later preallocate arrays with
% this size to avoid dynamic resizing.
npeps = 0;
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);

for j=1:length(matfiles)
    % Load each peptide-window .mat file (must contain 'His' and 'auc').
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1);  % loads variables into workspace: expects 'His', 'auc'
    
    % cnp = number of peptide entries in this window/file
    cnp = size(His.pep_mz,1);
    
    % Derive a human-readable peptide file name (without .mat)
    pepname = matfiles(j).name(1:end-4);
    % Some files may start with 'HH' (historical quirk); normalize to drop
    % one 'H' so that downstream string parsing stays consistent.
    if 1==strcmp(pepname(1:2),'HH')
        pepname = pepname(2:end);
    end;
    
    % Only count into 'npeps' when histone is H3 or H4.
    % (This reflects the downstream targets list which focuses on H3/H4 marks.)
    if 1==strcmp(pepname(1:2),'H3')
        npeps = npeps + cnp;
    elseif 1==strcmp(pepname(1:2),'H4')
        npeps = npeps + cnp;
    end;
end;

% -------------------------------------------------------------------------
% Build a cell array of peptide labels, one per peptide row, concatenating
% the position window + short modification name. This label will be later
% used to filter and sum ratios for each target (e.g., contains 'K4me1').
m = 0;
peptides = repmat({''},[npeps,1]);  % preallocate cell array of empty strings
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);

for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1);
    
    pepname = matfiles(j).name(1:end-4);
    if 1==strcmp(pepname(1:2),'HH')
        pepname = pepname(2:end);
    end;
    
    % 'pepname' encodes multiple parts separated by underscores. We extract
    % the "version" and build a standardized "histone_pos" label that
    % encodes histone type and positional window (e.g., 'H3_27_40').
    p = strfind(pepname,'_');
    version = pepname(p(1)+1:p(2)-1);
    
    % bvar = 1 indicates a "variant" naming branch (3-char version ending in 'v')
    bvar = 0;
    if length(version)>=4
        % Special-case for older naming 'H3_04v3...' → take char #4 or #4:end
        if 1==strcmp(pepname(1:7),'H3_04v3')
            version = version(4);
        else
            version = version(4:end);
        end;
        % Construct 'histone_pos' = <histone><version>_<start>_<end>
        histone_pos = [pepname(1:p(1)-1),version,'_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    elseif length(version)==3 && version(3)=='v'
        % Variant branch (example 'H2B' variants)
        bvar = 1;
        histone_pos = ['_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    else
        % Default: no explicit variant characters
        histone_pos = [pepname(1:p(1)-1),'_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    end;
    
    % Number of peptide entries in this file/window
    nlen = size(His.pep_mz,1);
    
    % Fill peptide labels for the current block [m+1 ... m+nlen]
    for ino = m+1:m+nlen
        if 1==bvar
            % For variant branch, strip the dot suffix from mod_short
            % (mod_short looks like 'Kxxme1.something'); keep left part.
            x = strfind(His.mod_short{ino-m},'.');
            if 1==strcmp(pepname(1:3),'H2B')
                % Keep explicit 'H2B' prefix if present.
                peptides{ino,1} = ['H2B',His.mod_short{ino-m}(1:x(1)-1),histone_pos];
            else
                peptides{ino,1} = [His.mod_short{ino-m}(1:x(1)-1),histone_pos];
            end;
        else
            % Standard branch: "<histone_pos> <mod_short>"
            peptides{ino,1} = [histone_pos,' ',His.mod_short{ino-m}];
        end;
    end;
    
    % Advance the running index 'm' by block size
    m = m+nlen;
    if m==npeps
        break;
    end;
end;

% -------------------------------------------------------------------------
% Build a big numeric matrix 'info' stacking, for each sample, the AUC info
% (RT, Area, Ratio) column-wise. The shape of 'info' is:
%   [npeps x (length(raw_names)*3)]
% such that for sample i, its 3 columns occupy [(i-1)*3+1 : (i-1)*3+3].
info = zeros([npeps,length(raw_names)*3]);

for i=1:length(raw_names)
    cur_rawname = raw_names{i};
    if i<10
        prefix = ['0',num2str(i)];
    else
        prefix = num2str(i);
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');
    
    m = 0;
    % We must iterate the SAME 'matfiles' order as above to maintain the
    % same peptide stacking order across samples (critical for alignment).
    for j = 1:length(matfiles)
        matfile1 = fullfile(cur_outpath,matfiles(j).name);
        load(matfile1);  % expects 'auc' with [RT, Area, Ratio]
        
        nlen = size(His.pep_mz,1);
        % Place the current block into the 3 columns for sample i
        info(m+1:m+nlen,(i-1)*3+1:(i-1)*3+3) = auc;
        m = m+nlen;
        if m==npeps
            break;
        end;
    end;
end;

% -------------------------------------------------------------------------
% Reorder columns into [all RT | all Area | all Ratio] groups to simplify
% later indexing. After this step, 'info' becomes:
%   info = [RT1..RTn | Area1..Arean | Ratio1..Ration], each block size = nSamples
info0 = zeros([npeps,length(raw_names)*3]);
for jno=1:length(raw_names)
    info0(:,jno) = info(:,(jno-1)*3+1);  % RT of sample jno
end;
for jno=1:length(raw_names)
    info0(:,length(raw_names)+jno) = info(:,(jno-1)*3+2);  % Area of sample jno
end;
for jno=1:length(raw_names)
    info0(:,2*length(raw_names)+jno) = info(:,(jno-1)*3+3);  % Ratio of sample jno
end;
info = info0;

% Convenience view: keep only the "Ratio" block (last third of columns).
% Shape: [npeps x nSamples]
info_ratio = info(1:npeps,2*length(raw_names)+1:3*length(raw_names));


%% output
% -------------------------------------------------------------------------
% Define the list of "single PTM" targets to summarize. Each string is
% interpreted with simple substring rules (see matching logic below).
targets = {'H3K4me1';
    'H3K4me2';
    'H3K4me3';
    'H3K4ac';
    'H3K9me1';
    'H3K9me2';
    'H3K9me3';
    'H3K9ac';
    'H3S10ph';
    'H3K14ac';
    'H3K18me1';
    'H3K18ac';
    'H3K23me1';
    'H3K23ac';
    'H31K27me1';
    'H31K27me2';
    'H31K27me3';
    'H31K27ac';
    'H31K36me1';
    'H31K36me2';
    'H31K36me3';
    'H33K27me1';
    'H33K27me2';
    'H33K27me3';
    'H33K27ac';
    'H33K36me1';
    'H33K36me2';
    'H33K36me3';
    'H3K56me1';
    'H3K56me2';
    'H3K56me3';
    'H3K56ac';
    'H3K79me1';
    'H3K79me2';
    'H3K79me3';
    'H3K79ac';
    'H3K122ac';
    'H4K5ac';
    'H4K8ac';
    'H4K12ac';
    'H4K16ac';
    'H4K20me1';
    'H4K20me2';
    'H4K20me3';
    'H4K20ac'};

% sratios will accumulate, for each target (row) and sample (column),
% the SUM of peptide-level ratios that match that target.
sratios = zeros([length(targets),length(raw_names)]);

for t=1:length(targets)
    c_pep = targets{t};  % current target string, e.g. 'H3K4me1'
    
    % -----------------------
    % Select a subset of peptide rows by histone variant when needed.
    % 'H31...' → restrict to peptides labeled like 'H3_27_40' (H3.1 window)
    % 'H33...' → restrict to peptides labeled like 'H33_27_40' (H3.3 window)
    if 1==strcmp(c_pep(1:3),'H31')
        X1 = [];
        no = 0;
        for p=1:npeps
            if 0==isempty(strfind(peptides{p},'H3_27_40'))
                no = no + 1;
                X1(no) = p;%#ok
            end;
        end;
        % For matching, we will search only the mod substring (drop 'H31')
        c_mod = c_pep(4:end);
        
    elseif 1==strcmp(c_pep(1:3),'H33')
        X1 = [];
        no = 0;
        for p=1:npeps
            if 0==isempty(strfind(peptides{p},'H33_27_40'))
                no = no + 1;
                X1(no) = p;%#ok
            end;
        end;
        c_mod = c_pep(4:end);
        
    else
        % Otherwise, consider all peptide rows (no variant restriction)
        X1 = 1:npeps;
        % Strip the 'H3' or 'H4' prefix so we match by modification only,
        % e.g. 'K4me1', 'K9ac', 'K79me2', etc.
        c_mod = c_pep(3:end);
    end;
    
    % -----------------------
    % Among the selected peptides X1, keep only those whose label contains
    % the modification substring (e.g., 'K4me1').
    X2 = [];
    no = 0;
    for x=1:length(X1)
        if 0==isempty(strfind(peptides{X1(x)},c_mod))
            no = no + 1;
            X2(no) = X1(x);%#ok
        end;
    end;
    
    % If no peptide matched this target, skip.
    if 0==no
        continue;
    end;
    
    % Sum ratios across all matched peptide rows, for each sample column.
    for x=1:length(X2)
        sratios(t,:) = sratios(t,:) + info_ratio(X2(x),:);
    end;
end;

% -------------------------------------------------------------------------
% Save the numeric result in a MATLAB file for downstream scripts.
mat_file = fullfile(layout_path,'histone_ratios_single_PTMs.mat');
save(mat_file,'targets','sratios');

% Also write a tab-separated table with '.xls' extension for convenient
% opening in Excel/LibreOffice. First column = target; then one column per
% sample, in the order given by 'raw_names'.
inten_file = fullfile(layout_path,'histone_ratios_single_PTMs.xls');
fp = fopen(inten_file,'w');
if -1==fp
    disp(['can not open: ',inten_file]);
    return;
end;

% Header line: tab + "index,name" per sample (kept exactly as original)
for i=1:length(raw_names)
    fprintf(fp,'\t%d,%s',i,raw_names{i});
end;
fprintf(fp,'\r\n');

% Body: each row is a target; then one value per sample (sratios)
for t=1:length(targets)
    fprintf(fp,'%s',targets{t});
    for jno=1:length(raw_names)
        fprintf(fp,'\t%f',sratios(t,jno));
    end;
    fprintf(fp,'\r\n');
end;

fclose(fp);
% ========================================================================
% End of function
% ========================================================================
