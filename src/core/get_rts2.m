function [rts, top1_rt, inten_sum, top1_inten_sum] = get_rts2( ...
    MS1_index, MS1_peaks, MS2_index, MS2_peaks, ...
    ptol, unitdiff, His, hno, nsplit, t1, t2, nhmass)
%% get_rts2 — Joint MS1+MS2 peak picking and RT selection for a target histone peptide
%
% PURPOSE (high level, non-expert friendly)
%   This function finds the retention time(s) (RT) where a *specific* peptide
%   (given by 'hno' inside the histone library 'His') is most likely eluting.
%   It does so by:
%     1) Extracting the peptide’s MS1 isotopic profile (XIC) in a given RT window.
%     2) Splitting the RT window into candidate peaks (elution features) and
%        quantifying their MS1 intensity.
%     3) Checking whether *diagnostic MS2 fragments* supporting this peptide are
%        present around each candidate apex.
%     4) Scoring the match between MS1 shape and MS2 fragments (cosine similarity)
%        and choosing the best-supported apex (top1).
%
% WHAT YOU GET OUT
%   rts            -> vector of candidate RTs (one per detected peak apex)
%   top1_rt        -> the single best RT (top-1 candidate) after combining MS1+MS2
%   inten_sum      -> MS1 intensity sum per candidate (same length as rts)
%   top1_inten_sum -> MS1 intensity sum for the chosen top-1 peak
%
% INPUTS (formats and meaning)
%   MS1_index : [N_MS1_scans x ?] numeric
%       Column 1: MS1 scan identifier (integer).
%       Column 2: retention time for the scan (e.g., minutes).  <-- used here
%       Other columns may exist but are not used directly in this function.
%
%   MS1_peaks : [sum(ms1_peaks_per_scan) x 2] numeric
%       Concatenated MS1 centroid peaks for all MS1 scans.
%       Column 1: m/z ; Column 2: intensity.
%       (The mapping “which rows belong to which scan” is handled by GetProfiles.)
%
%   MS2_index : [N_MS2_scans x ?] numeric
%       Column 1: MS2 scan identifier (integer).
%       Column 2: retention time of the MS2 scan (same time unit as MS1).  <-- used
%       Column 4: precursor m/z of the MS2 scan.                             <-- used
%       Column 6: instrument code (odd=ion trap; even=FT; 3/4=ETD family).   <-- used
%       Column 7: index pointer into MS2_peaks (row start for this scan).    <-- used
%       Column 8: noise/threshold estimate for this MS2 scan (used as 3×).   <-- used
%       (Other columns may exist; these conventions follow EpiProfile2 data.)
%
%   MS2_peaks : [sum(ms2_peaks_per_scan) x 2] numeric
%       Concatenated MS2 centroid peaks for all MS2 scans.
%       Column 1: m/z ; Column 2: intensity.
%
%   ptol : numeric
%       MS1 ppm tolerance (or a code) used by GetProfiles. If ptol==100, we
%       internally set it to 10. If ptol>100 and charge>=3, we allow ^13C handling.
%       (This mirrors EpiProfile2 heuristics.)
%
%   unitdiff : numeric
%       Mass difference per 13C (≈ 1.00335 Da). Used to place isotopic m/z traces:
%       [M-1, M, M+1, M+2] as [c_mz−unitdiff/z, c_mz, c_mz+unitdiff/z, c_mz+2*unitdiff/z].
%
%   His : struct
%       Histone peptide library with at least:
%         - pep_mz(:,1): monoisotopic m/z per peptide
%         - pep_ch(:,1): precursor charge per peptide
%
%   hno : integer
%       Index of the peptide of interest inside 'His'.
%
%   nsplit : integer
%       Peak splitting mode. (In this code path both branches call the same
%       helper GetTopBottom11; retained for compatibility.)
%
%   t1, t2 : numeric
%       Retention time search window [t1, t2]. Must intersect the MS1 RT range.
%       (Units must match MS1_index(:,2) and MS2_index(:,2).)
%
%   nhmass : integer (0/1)
%       Switch to choose the key-ion generator: get_key_ions1 (0) vs get_key_ions1H (1).
%       This relates to whether an N-terminal hydrogen mass offset is considered.
%
% ASSUMPTIONS / DEPENDENCIES
%   Requires the following helper functions (from EpiProfile core):
%     - GetProfiles(...)      : Extract multi-isotope MS1 traces in an RT index range.
%     - GetTopBottom11(...)   : Segment the MS1 trace into peak regions and pick apex.
%     - GetMods(...)          : Return modification catalog used to build diagnostic ions.
%     - get_key_ions1(...)    : Compute theoretical MS2 diagnostic ions (standard).
%     - get_key_ions1H(...)   : Same but with hydrogen-offset logic.
%   These are not defined below (only MatchMS2 and get_thr are), and must be available
%   in your MATLAB path.
%
% CAVEATS
%   - We assume exact precursor m/z matching to select MS2 scans (see MatchMS2).
%     If your MS2_index has small rounding differences, you may need to adapt that.
%   - If top1_rt < 4 (very early elution, often solvent front) we return empties.
%   - Time units: this code treats MS1/ MS2 RTs consistently but does not enforce units.
%     Ensure everything (t1, t2, MS1_index(:,2), MS2_index(:,2)) uses the same unit.

    %-----------%
    % 1) GUARDS %
    %-----------%

    num_MS1 = size(MS1_index, 1);

    % Fast exit if the requested window does not overlap the MS1 run.
    if t1 > MS1_index(num_MS1, 2) || t2 < 0 || t2 < t1
        rts = []; top1_rt = []; inten_sum = []; top1_inten_sum = [];
        return;
    end

    %------------------------------------------%
    % 2) Resolve peptide m/z and isotopic grid %
    %------------------------------------------%

    % Target peptide monoisotopic m/z and charge from the histone library.
    c_mz = His.pep_mz(hno, 1);
    c_ch = His.pep_ch(hno, 1);

    % Convert the requested RT limits into *row indices* within MS1_index.
    p     = find(MS1_index(:, 2) >= t1);
    rt_i1 = p(1);
    pp    = find(MS1_index(:, 2) <= t2);
    rt_i2 = pp(end);

    % Build the reference isotopic m/z channels we will extract in MS1:
    % M-1, M, M+1, M+2 (spaced by unitdiff/charge).
    c_ref_isomzs = [ ...
        c_mz - unitdiff / c_ch, ...
        c_mz, ...
        c_mz + unitdiff / c_ch, ...
        c_mz + 2 * unitdiff / c_ch];

    % Normalize/interpret 'ptol' the same way as upstream code:
    if ptol == 100
        ptol = 10;
    end
    if ptol > 100 && c_ch >= 3
        nC13 = 1;   % enable 13C handling for high charges if ptol looks "special"
    else
        nC13 = 0;
    end

    %-----------------------------------------%
    % 3) Extract MS1 isotope traces (GetProfiles)
    %-----------------------------------------%

    % GetProfiles returns:
    %   c_isorts          -> RT (one per MS1 scan in the slice [rt_i1:rt_i2])
    %   c_ref_isointens   -> [nScans x 4] intensity matrix for the 4 isotope channels
    [c_isorts, c_ref_isointens] = GetProfiles( ...
        MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2);

    % We will work on the *monoisotopic* trace (2nd column corresponds to M).
    c_mono_isointens = c_ref_isointens(:, 2);

    %-----------------------------------------------------%
    % 4) Split the RT slice into peak regions and get apex
    %    (GetTopBottom11 returns nt, nb, top1_idx, inten_sum)
    %-----------------------------------------------------%

    if nsplit == 1
        [nt, nb, top1_idx, inten_sum] = GetTopBottom11(c_mono_isointens);
    else
        % Same helper currently; retained for compatibility with original code.
        [nt, nb, top1_idx, inten_sum] = GetTopBottom11(c_mono_isointens);
    end

    % Candidate RTs (apex positions mapped to RT units).
    rts = c_isorts(nt);

    %--------------------------------------------------------%
    % 5) Combine MS1 shape with MS2 fragment evidence (score)
    %--------------------------------------------------------%

    % Get the modification table used to build theoretical fragment ions.
    Mods = GetMods();

    % 'flen' will track the number of expected diagnostic ion families (posn+posc).
    flen = 0;

    % Optionally smooth the monoisotopic trace (disabled here per original code).
    % w = smooth(c_mono_isointens, 3);
    w = c_mono_isointens;

    % Zero out very small intensities (below 1.5% of max) to denoise tails.
    w(find(w < 0.015 * max(w))) = 0; %#ok<FNDSB>

    % Accumulator of similarity scores (per candidate peak).
    similarity = zeros([length(nt), 1]);

    % For each candidate elution region:
    for ino = 1:length(nt)
        % Boundaries for this region (nb has length length(nt)+1)
        i1 = nb(ino);
        i2 = nb(ino + 1);

        % Shrink boundaries to the first/last *non-zero* points in 'w'
        while i1 <= length(w) && w(i1) == 0, i1 = i1 + 1; end
        if i1 > nb(ino), i1 = i1 - 1; end
        while i2 >= 1 && w(i2) == 0, i2 = i2 - 1; end
        if i2 < nb(ino + 1), i2 = i2 + 1; end

        IX  = i1:i2;
        rt1 = c_isorts(IX(1));   % RT start of the local region
        rt2 = c_isorts(IX(end)); % RT end   of the local region

        % Pull MS2 scans that match this peptide within [rt1, rt2]
        [ms2pos, ms2rts, ms2intens, posn, posc] = ...
            MatchMS2(MS2_index, MS2_peaks, Mods, His, hno, rt1, rt2, nhmass);

        % Track the expected number of diagnostic ion groups
        if flen < length(posn), flen = length(posn); end

        % If there are no supporting MS2 scans, skip scoring for this candidate
        if isempty(ms2pos), continue; end

        % Build a time window around the candidate apex based on the
        % *median spacing* of MS2 RTs (robust to variable duty cycles)
        nstep = median(ms2rts(ms2pos(2:end)) - ms2rts(ms2pos(1:end-1))) * 3.5;
        nst   = c_isorts(nt(ino)) - nstep;
        ntm   = c_isorts(nt(ino)) + nstep;

        % Select MS2 scans around the apex
        np1 = find(ms2rts(ms2pos) >= nst);
        np2 = find(ms2rts(ms2pos) <= ntm);
        X   = np1(1):np2(end);

        % Map those MS2 scans to the *closest* MS1 scans to pull the MS1 profile
        XX = zeros([1, length(X)]);
        for jno = 1:length(X)
            % MS2_index(ms2pos(X(jno)), 1)  -> MS2 scan id
            % Find the MS1 row whose scan id matches that MS2 scan id
            XX(jno) = find(MS1_index(:, 1) == MS2_index(ms2pos(X(jno)), 1));
        end

        % Extract the MS1 monoisotopic intensities at those aligned indices (then smooth)
        p_intens = c_mono_isointens(XX);
        p_intens = smooth(p_intens, 3);

        % Compute a dynamic similarity threshold and a tiny “baseline” vector
        % to avoid zero-division and overly optimistic cosine values.
        [sim_thr, ex_intens] = get_thr(p_intens);
        p_intens = p_intens + ex_intens;

        % Compare MS1 shape to each family of MS2 diagnostic ions (posn: neutral type)
        for kno = 1:length(posn)
            f_intens = ms2intens(ms2pos(X), kno);
            f_intens = smooth(f_intens, 3) + ex_intens;

            % Cosine similarity between MS1 and MS2 shapes
            e_sim = sum(p_intens .* f_intens) / ...
                sqrt(sum(p_intens .* p_intens) * sum(f_intens .* f_intens));

            % Only accumulate evidence if similarity beats the local threshold
            if e_sim > sim_thr
                similarity(ino) = similarity(ino) + e_sim;
            end
        end

        % Same for 'posc' (e.g., complementary / charge-dependent ion set)
        for kno = 1:length(posc)
            qno     = kno + length(posn);
            f_intens = ms2intens(ms2pos(X), qno);
            f_intens = smooth(f_intens, 3) + ex_intens;

            e_sim = sum(p_intens .* f_intens) / ...
                sqrt(sum(p_intens .* p_intens) * sum(f_intens .* f_intens));

            if e_sim > sim_thr
                similarity(ino) = similarity(ino) + e_sim;
            end
        end
    end

    %-----------------------------------------------%
    % 6) Choose the top-1 apex, using MS2 similarity
    %-----------------------------------------------%

    tmp = max(ceil(similarity));

    if tmp > flen / 2
        % If we have strong MS2 support (tmp large relative to expected ion families),
        % pick among the candidates achieving that max similarity the one with the
        % highest MS1 intensity sum.
        ii = find(ceil(similarity) == tmp);
        [~, id] = max(inten_sum(ii)); %#ok<ASGLU>
        top1_idx = ii(id);
        top1_rt = c_isorts(nt(top1_idx));
        top1_inten_sum = inten_sum(top1_idx);
    else
        % Otherwise, keep the apex chosen by GetTopBottom11 but “soften” the intensity.
        top1_rt = c_isorts(nt(top1_idx));
        top1_inten_sum = inten_sum(top1_idx) / 1e4;
    end

    %----------------------------------------------%
    % 7) Final guard: discard very early elutions  %
    %----------------------------------------------%

    if top1_rt < 4
        % Heuristic: < ~4 minutes is often solvent front / unreliable region.
        rts = []; top1_rt = []; inten_sum = []; top1_inten_sum = [];
        return;
    end
end

%===========================%
% Internal helper functions %
%===========================%

function [sim_thr, ex_intens] = get_thr(p_intens0)
%% get_thr — Build a conservative similarity threshold vs. “null” shapes
%
% WHY:
%   When comparing MS1 and MS2 shapes via cosine similarity, zeros or very
%   small vectors can inflate the score. We add tiny structured baselines (x)
%   and compute three cosine values against simple patterns; we then take the
%   MINIMUM as a per-window threshold. That is, a candidate must beat this
%   baseline to be considered “similar enough”.

    % Build three tiny (epsilon-level) patterns:
    x = zeros(length(p_intens0), 3);
    x(:, 1) = (ones(length(p_intens0), 1))          * eps; % flat
    x(:, 2) = ( (1:length(p_intens0))' )            * eps; % rising
    x(:, 3) = ( (length(p_intens0):-1:1)' )         * eps; % falling

    sim_thrs = zeros(1, 3);

    % Cosine vs pattern 1 (flat)
    p_intens = p_intens0 + x(:, 1);
    f_intens = x(:, 1);
    sim_thrs(1) = sum(p_intens .* f_intens) / ...
        sqrt(sum(p_intens .* p_intens) * sum(f_intens .* f_intens));

    % Cosine vs pattern 2 (rising)
    p_intens = p_intens0 + x(:, 2);
    f_intens = x(:, 2);
    sim_thrs(2) = sum(p_intens .* f_intens) / ...
        sqrt(sum(p_intens .* p_intens) * sum(f_intens .* f_intens));

    % Cosine vs pattern 3 (falling)
    p_intens = p_intens0 + x(:, 3);
    f_intens = x(:, 3);
    sim_thrs(3) = sum(p_intens .* f_intens) / ...
        sqrt(sum(p_intens .* p_intens) * sum(f_intens .* f_intens));

    % Final threshold is the *minimum* of the three (most conservative).
    [sim_thr, i] = min(sim_thrs);
    ex_intens    = x(:, i);   % the baseline we used (added to both series)
end

function [ms2pos, ms2rts, ms2intens, posn, posc, ActiveType, K1] = ...
    MatchMS2(MS2_index, MS2_peaks, Mods, His, hno, rt1, rt2, nhmass)
%% MatchMS2 — Collect MS2 scans for the target precursor and extract key-ion intensities
%
% WHAT THIS DOES
%   1) Finds MS2 scans whose precursor m/z exactly equals the target peptide’s m/z.
%   2) Keeps only those scans within the RT region [rt1, rt2].
%   3) Decides the fragmentation “ActiveType” from instrument code (ETD vs CID/HCD).
%   4) Builds the list of theoretical diagnostic ions (K1) for that peptide,
%      using either get_key_ions1 or get_key_ions1H (based on 'nhmass').
%   5) For each selected MS2 scan, finds the nearest measured peaks within a
%      mass tolerance (0.4 for ion trap, 0.02 for FT) and records their intensities.
%
% RETURNS
%   ms2pos     -> indices (rows in MS2_index) of selected scans
%   ms2rts     -> vector of all MS2 RTs (same length as MS2_index(:,2))
%   ms2intens  -> [N_MS2 x N_key_ions] intensities for the selected scans
%   posn,posc  -> partition of diagnostic ion indices (two families)
%   ActiveType -> 'ETD' or 'CID' (CID covers CID/HCD here)
%   K1         -> vector of theoretical m/z values for diagnostic ions

    %--------------------------------------------%
    % 1) Select MS2 scans for this precursor m/z %
    %--------------------------------------------%

    num_MS2 = size(MS2_index, 1);
    c_mz    = His.pep_mz(hno, 1);

    premzs = unique(MS2_index(:, 4));              % unique precursor m/z in the run
    [~, ii] = min(abs(premzs - c_mz));             % pick closest to the target m/z
    target  = premzs(ii);

    % Restrict to RT window first
    p  = find(MS2_index(:, 2) >= rt1);
    pp = find(MS2_index(:, 2) <= rt2);
    if isempty(p) || isempty(pp)
        ms2pos=[]; ms2rts=[]; ms2intens=[]; posn=[]; posc=[]; ActiveType=[]; K1=[];
        return
    end
    i1 = p(1); i2 = pp(end);

    % Flag MS2 rows whose precursor m/z matches the chosen target *exactly*.
    % NOTE: if your data has rounding errors, consider using a ppm window here.
    flag = zeros([num_MS2, 1]);
    for i = i1:i2
        if MS2_index(i, 4) == target
            flag(i) = 1;
        end
    end
    ms2pos = find(flag == 1);

    % If too few scans are available, abort MS2-based scoring.
    if length(ms2pos) <= 4
        ms2pos=[]; ms2rts=[]; ms2intens=[]; posn=[]; posc=[]; ActiveType=[]; K1=[];
        return
    end

    %-----------------------------%
    % 2) Determine fragmentation  %
    %-----------------------------%

    ms2rts     = MS2_index(:, 2);
    instruments = MS2_index(ms2pos, 6); % instrument codes per scan

    % Pick ActiveType and m/z tolerance (trap vs FT)
    c_instrument = instruments(1);
    if c_instrument == 3 || c_instrument == 4
        ActiveType = 'ETD';
    else
        ActiveType = 'CID';
    end
    if mod(c_instrument, 2) == 1
        tol = 0.4;   % ion trap-like
    else
        tol = 0.02;  % high-res FT-like
    end

    %----------------------------------------%
    % 3) Build diagnostic ion list (K1)      %
    %----------------------------------------%

    if nhmass == 1
        [K1, posn, posc] = get_key_ions1H(His, hno, Mods, ActiveType);
    else
        [K1, posn, posc] = get_key_ions1(His, hno, Mods, ActiveType);
    end

    %----------------------------------------%
    % 4) Gather intensities for those ions   %
    %----------------------------------------%

    % Row pointer into MS2_peaks: index(k) : index(k+1)-1 are peaks for scan k.
    index      = [1; MS2_index(1:num_MS2, 7)];
    ms2intens  = zeros([num_MS2, length(K1)]);

    for i = 1:length(ms2pos)
        cno = ms2pos(i);

        % (Optionally consider neighboring scans: here only the exact row 'cno')
        for pno = cno
            if pno < 1 || pno > num_MS2, continue; end

            % If instruments vary across scans, adapt ActiveType and tol on the fly
            if length(unique(instruments)) > 1
                c_instrument = MS2_index(pno, 6);
                if c_instrument == 3 || c_instrument == 4
                    ActiveType = 'ETD';
                else
                    ActiveType = 'CID';
                end
                if mod(c_instrument, 2) == 1
                    tol = 0.4;
                else
                    tol = 0.02;
                end

                % Also rebuild diagnostic ions if ActiveType changes
                if nhmass == 1
                    [K1, posn, posc] = get_key_ions1H(His, hno, Mods, ActiveType);
                else
                    [K1, posn, posc] = get_key_ions1(His, hno, Mods, ActiveType);
                end
            end

            % Pull the peaks for this MS2 scan
            IX    = index(pno):index(pno + 1) - 1;
            mz    = MS2_peaks(IX, 1);
            inten = MS2_peaks(IX, 2);

            % Basic noise filtering: keep peaks >= 3 × (scan’s noise estimate)
            I = find(inten >= 3 * MS2_index(pno, 8));
            mz = mz(I); inten = inten(I);

            % For each diagnostic ion m/z in K1, match the closest peak within tol
            for j = 1:length(K1)
                ix1 = find(abs(mz - K1(j)) <= tol);
                if ~isempty(ix1)
                    % pick the nearest in m/z (alternatively, max intensity)
                    [~, x1] = min(abs(mz(ix1) - K1(j)));
                    ms2intens(cno, j) = max(ms2intens(cno, j), inten(ix1(x1)));
                end
            end
        end
    end
end
