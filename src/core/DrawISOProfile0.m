function DrawISOProfile0(raw_path,raw_names,ptol,special)
%%
% DRAWISOPROFILE0  Build an initial reference (rt_ref) for histone peptides by mining MS1 data.
%
% Syntax:
%   DrawISOProfile0(raw_path, raw_names, ptol, special)
%
% Purpose:
%   This routine computes **reference retention times** (rt_ref) for all histone peptide entries
%   defined by `init_histone0(special)`. It scans MS1 indices/peaks from a *subset* of RAW runs
%   (top-3 selection heuristic) and, for each peptide, derives rt_ref as the **median** of the
%   observed RTs across those selected runs. The result is saved in:
%       <raw_path>/histone_layouts/0_ref_info.mat
%
% High-level workflow:
%   1) Safety checks:
%        - Require fewer than 100 runs and special.ndebug == 0; otherwise **return** (skip).
%        - Ensure '<raw_path>/histone_layouts' exists (create if needed).
%        - If '0_ref_info.mat' already exists, **return** to avoid overwriting an existing reference.
%   2) Initialize histone universe: AllUnHis = init_histone0(special).
%   3) **Run selection (top-3)**:
%        - If number of runs > 3, score each run by MS1 signal on 4 anchor entries (ipos=[1 2 3 4]).
%        - Compute ranks per anchor (higher intensity → higher rank), multiply ranks → a composite score.
%        - Pick the best run q; choose the anchor p with max intensity in q.
%        - Find other runs with similar RT (|ΔRT|<2) at anchor p; if ≥3, select left/middle/right
%          around the median RT with the highest intensities; else use all available (range = i).
%        - If ≤3 runs total, just use them all.
%   4) For each selected run and each peptide (length(AllUnHis.pep_mz)):
%        - Call get_rts0(MS1_index,MS1_peaks,ptol,AllUnHis,ino) → (RT, intensity).
%        - Store RTs and intensities; record selected_raws.
%   5) For each peptide, set AllUnHis.rt_ref = median(observed RTs across selected runs).
%   6) Save AllUnHis into '0_ref_info.mat'.
%
% Inputs:
%   - raw_path : char/string. Base path containing subfolders 'MS1' (with *_MS1scans.mat, *_MS1peaks.mat)
%                and (to be created/used) 'histone_layouts'.
%   - raw_names: cellstr (N×1). Basenames of RAW files (without extension); used to locate MS1 files.
%   - ptol     : double. Mass tolerance (ppm) for MS1 peak matching inside get_rts0.
%   - special  : struct. At least must contain 'ndebug':
%                  special.ndebug == 0  → proceed
%                  otherwise            → early return (skip).
%
% Outputs:
%   - None (side-effect): creates '<raw_path>/histone_layouts/0_ref_info.mat' with:
%       AllUnHis.rt_ref        : (numPeptides × 1) reference RTs (median across selected runs)
%       AllUnHis.selected_raws : (≤3 × 1) list of chosen runs
%       (and whatever fields init_histone0 populated, e.g., pep_mz, etc.)
%
% Dependencies (project-specific):
%   - init_histone0(special)       → returns AllUnHis with peptide definitions (fields: pep_mz, ...).
%   - get_rts0(MS1_index,MS1_peaks,ptol,AllUnHis,idx)
%       → for peptide index 'idx', returns (rt, intensity) from the given run.
%
% Data dependencies (per run):
%   - <raw_path>/MS1/<rawname>_MS1scans.mat  → defines 'MS1_index'
%   - <raw_path>/MS1/<rawname>_MS1peaks.mat  → defines 'MS1_peaks'
%
% File I/O behavior:
%   - If '<raw_path>/histone_layouts/0_ref_info.mat' **already exists**, the function **returns early**
%     to avoid overwriting an existing reference file.
%
% Notes / caveats:
%   - The selection heuristic uses 4 anchor peptide entries: ipos = [1 2 3 4], commented as:
%       H3_01_3_8, H3_02_9_17, H3_03_18_26, H3_04_27_40  (domain-specific anchors).
%   - Rank logic: higher intensity → higher rank; composite score = product of ranks across anchors.
%     The best run is the one with the **largest** composite score (sorted 'descend').
%   - RT window for similarity is fixed at **2 units** (likely minutes; domain-specific).
%   - For each peptide, rt_ref is set to the **median** RT across the selected runs (robust to outliers).
%   - This function assumes MS1 .mat files exist and contain variables named exactly 'MS1_index' and 'MS1_peaks'.
%

% check
if ~(length(raw_names)<100 && 0==special.ndebug)
    return;
end;
% Guard clause: only proceed when fewer than 100 runs AND ndebug==0.
% Otherwise, skip building a fresh reference (perhaps to avoid heavy work or when debugging).

layout_path = fullfile(raw_path,'histone_layouts');
if 0==exist(layout_path,'dir') && 0==mkdir(layout_path)
    fprintf(1,'can not create: %s\n',layout_path);
    return;
end;
% Ensure the target folder exists; if creation fails, abort.

mat_file = fullfile(layout_path,'0_ref_info.mat');
if 0~=exist(mat_file,'file')
    return;
end;
% If the reference file already exists, **do nothing** (prevent overwrite).

fprintf(1,'get a reference...\n');
AllUnHis = init_histone0(special);
% Initialize histone peptide definitions (m/z, sequences, etc.). Expected to define fields used below:
%   - pep_mz (vector/list length = number of peptide entries)
%   - possibly metadata for anchors referenced by ipos.

% get top 3 raw files
ntop = 3;
range = 1:length(raw_names);
if length(raw_names)>ntop
    ipos = [1 2 3 4];% H3_01_3_8, H3_02_9_17, H3_03_18_26, H3_04_27_40
    rts = zeros(length(ipos),length(raw_names));
    intens = zeros(length(ipos),length(raw_names));
    ranks = zeros(length(ipos),length(raw_names));
    score = ones(1,length(raw_names));
    for jno=1:length(raw_names)
        cur_rawname = raw_names{jno};
        MS1_scanfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1scans.mat']);
        MS1_peakfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1peaks.mat']);
        load(MS1_scanfile);% MS1_index
        load(MS1_peakfile);% MS1_peaks
        for ino=1:length(ipos)
            [rts(ino,jno),intens(ino,jno)] = get_rts0(MS1_index,MS1_peaks,ptol,AllUnHis,ipos(ino));
        end;
    end;
    % For each run jno and each anchor ino, collect RT and intensity via get_rts0.

    for ino=1:length(ipos)
        [tmp,xx] = sort(intens(ino,:));%#ok
        [tmp,ranks(ino,:)] = sort(xx);%#ok
        % Ranking trick:
        %   - 'xx' are indices of runs sorted by ascending intensity.
        %   - Sorting 'xx' gives, for each run index, its position in that ascending order.
        %   → Larger intensity ⇒ larger rank value.
        ii = find(intens(ino,:)==0);
        if 0==isempty(ii)
            ranks(ino,ii) = 1;
        end;
        % Zero intensity is forced to minimal rank (1), penalizing runs with missing signal.
    end;

    for jno=1:length(raw_names)
        for ino=1:length(ipos)
            score(jno) = score(jno)*ranks(ino,jno);
        end;
    end;
    % Composite score per run: product of ranks across the 4 anchors (higher is better).

    [tmp,ix] = sort(score,'descend');%#ok
    q = ix(1);
    % Best run index 'q' with the highest composite score.

    [tmp,p] = max(intens(:,q));%#ok
    % Choose the anchor 'p' with the strongest intensity in the best run q.

    [tmp,i] = find( abs(rts(p,:)-rts(p,q))<2 );%#ok
    % Collect runs whose anchor-p RT is within ±2 (time units) of that in run q.

    if length(i)<ntop
        range = i;
        % If fewer than 3 similar runs, use all of them.
    else
        [tmp,x] = min( abs(rts(p,i)-median(rts(p,i))) );%#ok
        % Pick the run closest to the median RT among the similar set as the "center" (i(x)).

        [tmp,l] = find(rts(p,i)<rts(p,i(x)));%#ok
        [tmp,r] = find(rts(p,i)>rts(p,i(x)));%#ok
        % Partition the neighbors into "left" (earlier RT) and "right" (later RT).

        [tmp,ll] = max(intens(p,i(l)));%#ok
        [tmp,rr] = max(intens(p,i(r)));%#ok
        % From left and right partitions, pick the runs with the highest anchor-p intensity.

        range = [i(l(ll)) i(x) i(r(rr))];
        % Final 3-run selection: [best-left, center (median-RT), best-right].
    end;
end;

% get rt_ref
info_rts = zeros(length(AllUnHis.pep_mz),length(range));
info_intens = zeros(length(AllUnHis.pep_mz),length(range));
for jno=1:length(range)
    cur_rawname = raw_names{range(jno)};
    MS1_scanfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1scans.mat']);
    MS1_peakfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1peaks.mat']);
    load(MS1_scanfile);% MS1_index
    load(MS1_peakfile);% MS1_peaks
    for ino=1:length(AllUnHis.pep_mz)
        [info_rts(ino,jno),info_intens(ino,jno)] = get_rts0(MS1_index,MS1_peaks,ptol,AllUnHis,ino);
    end;
    AllUnHis.selected_raws{jno,1} = cur_rawname;
end;
% For each selected run: load MS1 structures; for each peptide index 'ino':
%   - Query RT and intensity with get_rts0.
% Record the names of the selected runs in AllUnHis.selected_raws.

for ino=1:length(AllUnHis.pep_mz)
    %[tmp,ix] = max(info_intens(ino,:));%#ok
    %AllUnHis.rt_ref(ino,1) = info_rts(ino,ix);
    AllUnHis.rt_ref(ino,1) = median(info_rts(ino,:));
end;
% Reference RT per peptide: use the **median** across the selected runs.
% (There is a commented alternative: pick the RT from the run where the peptide has max intensity.)

save(mat_file,'AllUnHis');
% Persist the initialized and populated AllUnHis (including rt_ref and selected_raws)
% into '<raw_path>/histone_layouts/0_ref_info.mat'.
