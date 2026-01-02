function [cur_rts,cur_intens,cur_mono_isointens] = get_histone6( ...
    MS1_index, MS1_peaks, MS2_index, MS2_peaks, ...
    ptol, unitdiff, Mods, His, hno, special)
%GET_HISTONE6
% -------------------------------------------------------------------------
% Purpose (non-proteomics wording):
%   Decompose the MS1 monoisotopic signal of a single peptide into SIX
%   overlapping "isoform" components (e.g., six different PTM states).
%   We keep the TOTAL quantity from MS1 (area) and use MS2 diagnostic ions
%   to estimate the RELATIVE MIX (fractions) along the chromatographic peak.
%
% Inputs (minimal mental model):
%   - MS1_index [N x ?]: Per-scan metadata for MS1.
%       MS1_index(:,1) : MS1 scan id (integer)
%       MS1_index(:,2) : retention time (RT, minutes)
%       MS1_index(:,3) : row offset into MS1_peaks for this scan (used by GetProfiles)
%   - MS1_peaks  [M x 2]: All MS1 peaks concatenated by scans (mz, intensity).
%   - MS2_index [K x ?]: Per-MS2 spectrum metadata.
%       MS2_index(:,1) : parent MS1 scan id
%       MS2_index(:,2) : RT (minutes)
%       MS2_index(:,4) : precursor m/z (or DIA channel center m/z)
%       MS2_index(:,6) : instrument type code (odd=ion trap (wide tol), even=FT (narrow tol))
%       MS2_index(:,7) : row offset into MS2_peaks for this MS2 scan
%   - MS2_peaks  [P x 2]: All MS2 peaks concatenated (mz, intensity).
%   - ptol (ppm): m/z tolerance in ppm (special rule: if 100 → use 10).
%   - unitdiff (Da): isotopic spacing (for 13C ~ 1.003355 Da).
%   - Mods: structure defining diagnostic ions per isoform pair.
%   - His: structure with peptide/isoform metadata:
%       .pep_mz(hno:hno+5, :)  → m/z per charge state for each isoform
%       .pep_ch(hno:hno+5, :)  → charge for each isoform
%       .rt_ref(hno:hno+5)     → reference RTs (minutes)
%       .mod_short{hno:hno+5}  → short labels for isoforms (strings)
%       .outpath, .outfile     → output folder/name for QC PDF
%   - hno: index of the FIRST isoform in the block of SIX (hno..hno+5).
%   - special: flags/modes
%       .nDAmode  : 1 = DDA, 2 = DIA
%       .nsource, .nsubtype : enable hydrogen-mass variants for key ions
%
% Outputs:
%   - cur_rts         [6 x nCharge]: apex RT (min) for each isoform (replicated across charges).
%   - cur_intens      [6 x nCharge]: area per isoform/charge (sums ~ total MS1 area).
%   - cur_mono_isointens [nMS1Scans x 6]: time-resolved monoisotopic sub-traces per isoform.
%
% Dependencies (provided elsewhere in the codebase):
%   GetProfiles, GetTopBottom11 (or GetTopBottom), smooth, spline,
%   get_key_ions / get_key_ionsH, get_key_ions2 / get_key_ions2H.
%
% Design notes:
%   1) Obtain TOTAL area from MS1 via get_histone1 (robust, charge-aware).
%   2) In a small RT window covering the six isoforms, extract the MS1
%      monoisotopic trace.
%   3) Use MS2 diagnostic ions to compute five pairwise ratios that inform
%      the six-way composition (per scan). Solve a constrained system to
%      get six non-negative fractions summing to 1.
%   4) Multiply the experimental monoisotopic trace by those fractions to
%      obtain six sub-traces; integrate each to obtain areas and apex RTs.
%   5) Save a QC PDF overlaying the experimental trace and the six sub-traces.
%
% Behaviour preserved: identical logic/thresholds/fallbacks as the original.
% Only documentation and structure were improved for readability.
% -------------------------------------------------------------------------

    %% ---- Initialization: fixed-size outputs and a safety check ----------
    npart = 6;  % number of isoforms in this block
    [~, ncharge] = size(His.pep_mz);
    cur_rts    = zeros([npart, ncharge]);           % apex RT per isoform
    cur_intens = zeros([npart, ncharge]);           % area per isoform
    num_MS1    = size(MS1_index, 1);
    cur_mono_isointens = zeros([num_MS1, npart]);   % sub-traces on MS1 grid

    %% ---- Step 1: Total MS1 area via get_histone1 ------------------------
    % get_histone1 returns a single-component RT and area distribution per charge.
    [h_rts, h_intens] = get_histone1(MS1_index, MS1_peaks, ptol, unitdiff, His, hno);
    if h_intens(1) == 0
        % Nothing to deconvolve (no measurable signal) → return zeros.
        return;
    end

    %% ---- Step 2–8: Estimate six-way split and build sub-traces ----------
    [s_rts, ratio, cur_mono_isointens] = get_ratio_6iso( ...
        MS1_index, MS1_peaks, MS2_index, MS2_peaks, ...
        ptol, unitdiff, Mods, His, hno, h_rts, special);

    % Replicate the six apex RTs across charge states
    cur_rts(1,1:ncharge) = repmat(s_rts(1), [1, ncharge]);
    cur_rts(2,1:ncharge) = repmat(s_rts(2), [1, ncharge]);
    cur_rts(3,1:ncharge) = repmat(s_rts(3), [1, ncharge]);
    cur_rts(4,1:ncharge) = repmat(s_rts(4), [1, ncharge]);
    cur_rts(5,1:ncharge) = repmat(s_rts(5), [1, ncharge]);
    cur_rts(6,1:ncharge) = repmat(s_rts(6), [1, ncharge]);

    % Scale the total MS1 area per charge by each isoform fraction
    cur_intens(1,1:ncharge) = h_intens * ratio(1);
    cur_intens(2,1:ncharge) = h_intens * ratio(2);
    cur_intens(3,1:ncharge) = h_intens * ratio(3);
    cur_intens(4,1:ncharge) = h_intens * ratio(4);
    cur_intens(5,1:ncharge) = h_intens * ratio(5);
    cur_intens(6,1:ncharge) = h_intens * ratio(6);

    % Done. Plotting/QC is handled inside get_ratio_6iso (to keep inputs/outputs
    % of this wrapper clean and not pollute the namespace).

end


% =========================================================================
% INTERNAL: Compute six-way ratios and sub-traces using MS2 diagnostics
% =========================================================================
function [s_rts, ratio, cur_mono_isointens] = get_ratio_6iso( ...
    MS1_index, MS1_peaks, MS2_index, MS2_peaks, ...
    ptol, unitdiff, Mods, His, hno, h_rts, special)
%GET_RATIO_6ISO
%   Builds the monoisotopic MS1 trace in a small RT window that covers the
%   six isoforms; then estimates, for each MS1 scan, the six-way fraction
%   mix using co-temporal MS2 diagnostic ions. Finally, it integrates each
%   of the six sub-traces to produce isoform-specific areas and apex RTs.

    npart = 6;
    s_rts  = zeros([1, npart]);                % apex RTs (one per isoform)
    ratio  = zeros([1, npart]);                % global isoform fractions (sum=1)
    num_MS1 = size(MS1_index, 1);
    cur_mono_isointens = zeros([num_MS1, npart]);

    % --- RT window: small, centered to cover hno..hno+5
    delta = 0.5;
    p  = find(MS1_index(:,2) >= His.rt_ref(hno)   - delta);
    rt_i1 = p(1);
    pp = find(MS1_index(:,2) <= His.rt_ref(hno+5) + delta);
    rt_i2 = pp(end);

    % --- Extract monoisotopic MS1 trace (use reference isotopic pattern)
    c_mz  = His.pep_mz(hno,1);  % anchor at first isoform, first charge column
    c_ch  = His.pep_ch(hno,1);
    c_ref_isomzs = [c_mz - unitdiff/c_ch, c_mz, c_mz + unitdiff/c_ch, c_mz + 2*unitdiff/c_ch];

    if ptol == 100
        ptol = 10;  % legacy behaviour
    end
    nC13 = (ptol > 100 && c_ch >= 3);  % use +1 C13 model if very wide tol and high charge

    [c_isorts, c_ref_isointens] = GetProfiles( ...
        MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2);

    j = 2;  % index of the monoisotopic column within c_ref_isointens
    c_mono_isointens = c_ref_isointens(:, j);

    % --- MS2 selection in this RT window (different for DDA vs DIA)
    IX  = rt_i1:rt_i2;
    rt1 = c_isorts(IX(1));
    rt2 = c_isorts(IX(end));

    if special.nDAmode == 2
        % DIA: look for the MS2 "channel" closest to c_mz and lock to it
        premzs = unique(MS2_index(:,4));
        [~, ii] = min(abs(premzs - c_mz));
        target  = premzs(ii);

        num_MS2 = size(MS2_index, 1);
        flag = zeros([num_MS2, 1]);
        p  = find(MS2_index(:,2) >= rt1);
        if isempty(p), return; end
        i1 = p(1);
        pp = find(MS2_index(:,2) <= rt2);
        if isempty(pp), return; end
        i2 = pp(end);

        for i = i1:i2
            cen_mz = MS2_index(i,4);
            if cen_mz - target == 0
                flag(i) = 1;
            end
        end

        ms2pos = find(flag == 1);
        if isempty(ms2pos), return; end

        IX = zeros([length(ms2pos), 1]);
        for i = 1:length(ms2pos)
            IX(i) = find(MS1_index(:,1) == MS2_index(ms2pos(i),1));
        end
        rt1 = c_isorts(IX(1));
        rt2 = c_isorts(IX(end));
    end

    % --- Match MS2 and compute per-MS2 six-way compositions
    [ms2pos, ms2ratios] = MatchMS2( ...
        MS2_index, MS2_peaks, c_mz, c_ch, ptol, unitdiff, Mods, His, hno, rt1, rt2, special);

    % --- Fallback: too few MS2? Use equal split (visual QA only).
    nlen = length(IX);
    if length(ms2pos) <= 3
        if nlen <= 3, return; end
        s_rts(1:npart)  = [His.rt_ref(hno), His.rt_ref(hno+1), His.rt_ref(hno+2), ...
                           His.rt_ref(hno+3), His.rt_ref(hno+4), His.rt_ref(hno+5)];
        ratio(1:npart)  = repmat(1/npart, [1, npart]);
        n_inten1        = c_mono_isointens(IX) / npart;

        if special.nDAmode ~= 2
            % DDA: place equally split trace directly on the MS1 grid
            for k = 1:npart
                cur_mono_isointens(:,k) = [zeros([IX(1)-1,1]); n_inten1; zeros([num_MS1-IX(end),1])];
            end
        else
            % DIA: spline-interpolate split trace back onto MS1 grid
            xx = c_isorts(IX(1):IX(end));
            yy = spline(c_isorts(IX), n_inten1, xx);
            for k = 1:npart
                cur_mono_isointens(:,k) = [zeros([IX(1)-1,1]); yy; zeros([num_MS1-IX(end),1])];
            end
        end
        return;
    end

    % --- Build per-MS1-scan ratio matrix (six columns), via MS2-to-MS1 mapping
    ratio1  = zeros([nlen, npart]);      % rows = MS1 scans within window
    ms1scans = MS2_index(ms2pos, 1);     % parent MS1 scan ids hosting the MS2 events

    for i = 1:nlen
        c_scan = MS1_index(IX(i), 1);
        [tf, loc] = ismember(c_scan, ms1scans);
        if tf
            % Copy six-way composition measured at that co-temporal MS2
            ratio1(i, 1:npart) = ms2ratios(1:npart, loc)';
        end
    end

    % --- Smooth columns to reduce per-scan noise; then row-normalize
    for k = 1:npart
        ratio1(:,k) = smooth(ratio1(:,k), 3);
    end
    sumrow = sum(ratio1, 2);
    for k = 1:npart
        ratio1(:,k) = ratio1(:,k) ./ (eps + sumrow);
    end

    % --- Split the experimental mono trace into six sub-traces
    n_rt     = c_isorts(IX);
    n_inten1 = c_mono_isointens(IX);
    n_inten2 = n_inten1 .* ratio1(:,1);
    n_inten3 = n_inten1 .* ratio1(:,2);
    n_inten4 = n_inten1 .* ratio1(:,3);
    n_inten5 = n_inten1 .* ratio1(:,4);
    n_inten6 = n_inten1 .* ratio1(:,5);
    n_inten7 = n_inten1 .* ratio1(:,6);   % note: variable name keeps original style

    % Apex indices (per isoform)
    [~, x1] = max(n_inten2);
    [~, x2] = max(n_inten3);
    [~, x3] = max(n_inten4);
    [~, x4] = max(n_inten5);
    [~, x5] = max(n_inten6);
    [~, x6] = max(n_inten7);

    % --- Areas via spline integration on a fine RT mesh
    xx = n_rt(1):0.005:n_rt(end);
    yy = spline(n_rt, n_inten2, xx); area1 = sum(abs(yy));
    yy = spline(n_rt, n_inten3, xx); area2 = sum(abs(yy));
    yy = spline(n_rt, n_inten4, xx); area3 = sum(abs(yy));
    yy = spline(n_rt, n_inten5, xx); area4 = sum(abs(yy));
    yy = spline(n_rt, n_inten6, xx); area5 = sum(abs(yy));
    yy = spline(n_rt, n_inten7, xx); area6 = sum(abs(yy));
    sumarea = eps + area1 + area2 + area3 + area4 + area5 + area6;

    % --- Outputs: apex RTs and global fractions (sum=1)
    s_rts(1:npart) = [n_rt(x1), n_rt(x2), n_rt(x3), n_rt(x4), n_rt(x5), n_rt(x6)];
    ratio(1:npart) = [area1, area2, area3, area4, area5, area6] ./ sumarea;

    % --- Write sub-traces back onto the full MS1 scan grid
    if special.nDAmode ~= 2
        % DDA: direct placement
        cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]); n_inten2; zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]); n_inten3; zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,3) = [zeros([IX(1)-1,1]); n_inten4; zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,4) = [zeros([IX(1)-1,1]); n_inten5; zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,5) = [zeros([IX(1)-1,1]); n_inten6; zeros([num_MS1-IX(end),1])];
        cur_mono_isointens(:,6) = [zeros([IX(1)-1,1]); n_inten7; zeros([num_MS1-IX(end),1])];
    else
        % DIA: spline back to the exact MS1 RT locations within [IX(1), IX(end)]
        xx = c_isorts(IX(1):IX(end));
        yy = spline(n_rt, n_inten2, xx);
        cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]); yy; zeros([num_MS1-IX(end),1])];
        yy = spline(n_rt, n_inten3, xx);
        cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]); yy; zeros([num_MS1-IX(end),1])];
        yy = spline(n_rt, n_inten4, xx);
        cur_mono_isointens(:,3) = [zeros([IX(1)-1,1]); yy; zeros([num_MS1-IX(end),1])];
        yy = spline(n_rt, n_inten5, xx);
        cur_mono_isointens(:,4) = [zeros([IX(1)-1,1]); yy; zeros([num_MS1-IX(end),1])];
        yy = spline(n_rt, n_inten6, xx);
        cur_mono_isointens(:,5) = [zeros([IX(1)-1,1]); yy; zeros([num_MS1-IX(end),1])];
        yy = spline(n_rt, n_inten7, xx);
        cur_mono_isointens(:,6) = [zeros([IX(1)-1,1]); yy; zeros([num_MS1-IX(end),1])];
    end

    % --- QC plot (non-intrusive: figures stay hidden, auto-saved as PDF)
    set(gcf, 'visible', 'off');
    xx_ac = strfind(His.mod_short{hno}, 'ac');  %#ok<NASGU> % preserved behaviour
    out_file1 = fullfile(His.outpath, ['Iso_', His.outfile, '_', num2str(length(xx_ac)), 'ac', '_', His.mod_short{hno}, '.pdf']);

    plot(n_rt, n_inten1, 'linestyle','-','linewidth',2,'color','k'); hold on;
    plot(n_rt, n_inten2, 'linestyle','--','linewidth',2,'color','r');
    plot(n_rt, n_inten3, 'linestyle','-.','linewidth',2,'color','b');
    plot(n_rt, n_inten4, 'linestyle','--','linewidth',2,'color','g');
    plot(n_rt, n_inten5, 'linestyle','-.','linewidth',2,'color','m');
    plot(n_rt, n_inten6, 'linestyle','--','linewidth',2,'color','c');
    plot(n_rt, n_inten7, 'linestyle','-.','linewidth',2,'color','y');
    xlabel('time (min)'); ylabel('intensity');
    legend('experiment', His.mod_short{hno}, His.mod_short{hno+1}, His.mod_short{hno+2}, ...
           His.mod_short{hno+3}, His.mod_short{hno+4}, His.mod_short{hno+5});

    r = floor(100 * ratio);                       % integer percentages
    r(6) = 100 - sum(r(1:5));                     % keep sum close to 100
    title([His.mod_short{hno},'/',His.mod_short{hno+1},'/',His.mod_short{hno+2},'/', ...
           His.mod_short{hno+3},'/',His.mod_short{hno+4},'/',His.mod_short{hno+5}, ...
           ':', num2str(r(1)),'%:', num2str(r(2)),'%:', num2str(r(3)),'%:', ...
           num2str(r(4)),'%:', num2str(r(5)),'%:', num2str(r(6)),'%']);
    print('-dpdf', out_file1);
    close();

end


% =========================================================================
% INTERNAL: Find MS2 in the RT window; compute six-way compositions per MS2
% =========================================================================
function [ms2pos, ms2ratios] = MatchMS2(MS2_index, MS2_peaks, c_mz, c_ch, ...
    ptol, unitdiff, Mods, His, hno, rt1, rt2, special)
%MATCHMS2
%   Select MS2 spectra within the given RT window, match the expected
%   precursor(s) (DDA) or DIA channel, and compute per-MS2 pairwise ratios
%   using diagnostic ions. Then solve for a six-way composition vector X.

    % --- Hydrogen-mass variant selector (preserved)
    if 4 == special.nsource && (special.nsubtype ~= 0 && special.nsubtype ~= 2)
        nhmass = 1;
    else
        nhmass = 0;
    end

    npart   = 6;
    num_MS2 = size(MS2_index, 1);

    % --- Select MS2 in RT window, matching precursor/channel
    if special.nDAmode == 1
        % DDA: allow M, M+1, M+2 isotopologues (same charge)
        sets = [0, 1, 2];
        mzs  = c_mz + sets * unitdiff / c_ch;
        flag = zeros([num_MS2, 1]);

        p = find(MS2_index(:,2) >= rt1);  if isempty(p), ms2pos=[]; ms2ratios=[]; return; end
        i1 = p(1);
        pp = find(MS2_index(:,2) <= rt2); if isempty(pp), ms2pos=[]; ms2ratios=[]; return; end
        i2 = pp(end);

        for i = i1:i2
            cen_mz = MS2_index(i,4);
            ix = find(abs(mzs - cen_mz) < ptol * cen_mz * 1e-6);
            if ~isempty(ix), flag(i) = 1; end
        end

        ms2pos = find(flag == 1);
        if isempty(ms2pos), ms2ratios=[]; return; end

    elseif special.nDAmode == 2
        % DIA: exact match on channel m/z (fixed window)
        premzs = unique(MS2_index(:,4));
        [~, ii] = min(abs(premzs - c_mz));
        target  = premzs(ii);

        flag = zeros([num_MS2, 1]);
        p = find(MS2_index(:,2) >= rt1);  if isempty(p), ms2pos=[]; ms2ratios=[]; return; end
        i1 = p(1);
        pp = find(MS2_index(:,2) <= rt2); if isempty(pp), ms2pos=[]; ms2ratios=[]; return; end
        i2 = pp(end);

        for i = i1:i2
            if MS2_index(i,4) - target == 0
                flag(i) = 1;
            end
        end

        ms2pos = find(flag == 1);
        if isempty(ms2pos), ms2ratios=[]; return; end

    else
        ms2pos   = [];
        ms2ratios = [];
        return;
    end

    % --- Instrument model: define activation and mass tolerance per scan
    instruments = MS2_index(ms2pos, 6);
    if length(unique(instruments)) == 1
        c_instrument = instruments(1);
        if c_instrument == 3 || c_instrument == 4
            ActiveType = 'ETD';
        else
            ActiveType = 'CID';
        end
        if mod(c_instrument, 2) == 1
            tol = 0.4;    % ion trap-like
        else
            tol = 0.02;   % FT-like
        end

        % --- Diagnostic ion sets for five pairwise comparisons
        if nhmass
            [K11,K12] = get_key_ionsH (His, hno,   hno+1, Mods, ActiveType);
            [K21,K22] = get_key_ionsH (His, hno+1, hno+2, Mods, ActiveType);
            [K31,K32] = get_key_ionsH (His, hno+2, hno+4, Mods, ActiveType);
            [K41,K42] = get_key_ionsH (His, hno+4, hno+5, Mods, ActiveType);
            [K51,K52] = get_key_ions2H(His, hno+5, hno+2, Mods, ActiveType);
        else
            [K11,K12] = get_key_ions  (His, hno,   hno+1, Mods, ActiveType);
            [K21,K22] = get_key_ions  (His, hno+1, hno+2, Mods, ActiveType);
            [K31,K32] = get_key_ions  (His, hno+2, hno+4, Mods, ActiveType);
            [K41,K42] = get_key_ions  (His, hno+4, hno+5, Mods, ActiveType);
            [K51,K52] = get_key_ions2 (His, hno+5, hno+2, Mods, ActiveType);
        end
    end

    % --- Iterate MS2 scans and compute six-way composition vectors
    index = [1; MS2_index(1:num_MS2, 7)];
    ms2ratios = zeros([npart, length(ms2pos)]);

    for i = 1:length(ms2pos)
        cno    = ms2pos(i);
        newpos = cno;                 % preserved behaviour (single scan, not ±1 neighbourhood)

        for pno = newpos
            if pno < 1 || pno > num_MS2, continue; end

            % If instruments are mixed, adapt per-scan
            if length(unique(instruments)) > 1
                c_instrument = MS2_index(pno,6);
                if c_instrument == 3 || c_instrument == 4
                    ActiveType = 'ETD';
                else
                    ActiveType = 'CID';
                end
                if mod(c_instrument,2) == 1
                    tol = 0.4;
                else
                    tol = 0.02;
                end

                if nhmass
                    [K11,K12] = get_key_ionsH (His, hno,   hno+1, Mods, ActiveType);
                    [K21,K22] = get_key_ionsH (His, hno+1, hno+2, Mods, ActiveType);
                    [K31,K32] = get_key_ionsH (His, hno+2, hno+4, Mods, ActiveType);
                    [K41,K42] = get_key_ionsH (His, hno+4, hno+5, Mods, ActiveType);
                    [K51,K52] = get_key_ions2H(His, hno+5, hno+2, Mods, ActiveType);
                else
                    [K11,K12] = get_key_ions  (His, hno,   hno+1, Mods, ActiveType);
                    [K21,K22] = get_key_ions  (His, hno+1, hno+2, Mods, ActiveType);
                    [K31,K32] = get_key_ions  (His, hno+2, hno+4, Mods, ActiveType);
                    [K41,K42] = get_key_ions  (His, hno+4, hno+5, Mods, ActiveType);
                    [K51,K52] = get_key_ions2 (His, hno+5, hno+2, Mods, ActiveType);
                end
            end

            % Extract current MS2 spectrum peaks
            IX = index(pno):index(pno+1)-1;
            mz    = MS2_peaks(IX, 1);
            inten = MS2_peaks(IX, 2);

            % Five pairwise ratios (sum(A)/(sum(A)+sum(B)) for each pair)
            r1 = get_ratio(mz, inten, tol, K11, K12);
            r2 = get_ratio(mz, inten, tol, K21, K22);
            r3 = get_ratio(mz, inten, tol, K31, K32);
            r4 = get_ratio(mz, inten, tol, K41, K42);
            r5 = 1.1 * get_ratio(mz, inten, tol, K51, K52);  % intentional scaling

            % Solve for the six-way composition X = [a b c d e f]'
            X = resove_linear(r1, r2, r3, r4, r5);

            % Keep the "most informative" composition (fewest zeros replaced)
            if length(find(ms2ratios(1:npart, i) == 0)) > length(find(X == 0))
                ms2ratios(1:npart, i) = X;
            end
        end
    end

    % --- Guard: if no fragments resolved, set a minimal non-degenerate vector
    x = find(ms2ratios(npart, :) < 1);  %#ok<NASGU> % preserved check
    if isempty(x)
        ms2pos    = [];
        ms2ratios = [];
    else
        xx = find(ms2ratios(npart, :) == 1);
        if ~isempty(xx)
            for ino = 1:length(xx)
                ms2ratios(1:npart, xx(ino)) = [0, 0, 0, 0, 0.01, 0.99]';
            end
        end
    end
end


% =========================================================================
% INTERNAL: Pairwise ratio (sum(A) / (sum(A)+sum(B))) using diagnostic ions
% =========================================================================
function ratio1 = get_ratio(mz, inten, tol, K1, K2)
%GET_RATIO
%   Given an MS2 spectrum (mz, inten) and two diagnostic ion lists (K1,K2),
%   return sum(intensity@K1) / (sum(K1)+sum(K2)) with per-ion nearest-peak
%   matching under tolerance "tol".

    intens1 = zeros([1, length(K1)]);
    intens2 = zeros([1, length(K1)]);
    for j = 1:length(K1)
        ix1 = find(abs(mz - K1(j)) <= tol);
        ix2 = find(abs(mz - K2(j)) <= tol);
        [~, x1] = min(abs(mz(ix1) - K1(j)));  % nearest peak to K1(j)
        [~, x2] = min(abs(mz(ix2) - K2(j)));  % nearest peak to K2(j)

        if ~isempty(ix1), intens1(j) = inten(ix1(x1)); end
        if ~isempty(ix2), intens2(j) = inten(ix2(x2)); end
    end
    ratio1 = sum(intens1) / (eps + sum(intens1) + sum(intens2));
end


% =========================================================================
% INTERNAL: Piece-wise constrained solver for six-way composition
% =========================================================================
function X = resove_linear(r1, r2, r3, r4, r5)
%RESOVE_LINEAR
%   Recover a six-component non-negative vector X = [a b c d e f]' that
%   sums to 1 and is consistent with five pairwise ratios r1..r5.
%   This is intentionally a piece-wise solver mirroring the original logic
%   (kept verbatim for behaviour equivalence).

    if r2 == 0
        a=0; b=0; d=0; c=r3;
        if r3 >= r4, e=0;        f=1-r3;
        else         e=r4-r3;    f=1-r4; end
        X=[a;b;c;d;e;f];

    elseif r3 == 0
        a=0; b=0; c=0; d=r2;
        if r2 >= r4, e=0;        f=1-r2;
        else         e=r4-r2;    f=1-r4; end
        X=[a;b;c;d;e;f];

    elseif r4 == 0
        a=r2; b=0; c=0; d=0; e=0; f=1-r2;
        X=[a;b;c;d;e;f];

    elseif r1 == 0
        a=0;
        if r3 >= r4
            d=0; e=0; b=r2; c = max(r3-r2,0); f=1-b-c;
        elseif r2 >= r4
            c=0; e=0; b=r3; d = max(r2-r3,0); f=1-b-d;
        else
            if r2 > r3
                c=0; b=r3; d=r2-r3; e=r4-r2;
            else
                d=0; b=r2; c=r3-r2; e=r4-r3;
            end
            f=1-r4;
        end
        X=[a;b;c;d;e;f];

    elseif r5 == 0
        a=r1; f=0;
        if r1 >= r2
            b=0; d=0;
            if r1 >= r3
                c=0; e=1-r1;
            else
                c=r3-r1; e=1-r3;
            end
        elseif r1 >= r3
            b=0; c=0; d=r2-r1; e=1-r2;
        else
            % reduce similar to r1==0 case
            r2=r2-r1; r3=r3-r1; r4=1-r1;
            if r3 >= r4
                d=0; e=0; b=r2; c = (r1+b>=1) * 0 + (r1+b<1) * (1-r1-b);
            elseif r2 >= r4
                c=0; e=0; b=r3; d = (r1+b>=1) * 0 + (r1+b<1) * (1-r1-b);
            else
                if r2 > r3
                    c=0; b=r3; d=r2-r3; e=r4-r2;
                else
                    d=0; b=r2; c=r3-r2; e=r4-r3;
                end
                a=1-r4;
            end
        end
        X=[a;b;c;d;e;f];

    else
        % general (all > 0)
        a = r4/(r4 + 1/r1 - 1);
        f = 1 - r4/(r1*r4 + 1 - r1);

        if a + f >= 1
            b=0; c=0; d=0; e=0; f=1-a;

        elseif a >= r2
            b=0; d=0;
            if a >= r3
                c=0; e=1-a-f;
            else
                c = r3 - a;
                if a + f + c >= 1
                    c = 1 - a - f; e=0;
                else
                    e = 1 - a - f - c;
                end
            end

        elseif a >= r3
            b=0; c=0; d = r2 - a;
            if a + f + d >= 1
                d = 1 - a - f; e=0;
            else
                e = 1 - a - f - d;
            end

        else
            % reduce similar to r1==0 case
            r2=r2-a; r3=r3-a; r4=1-a-f;
            if r3 >= r4
                d=0; e=0; b=r2; c = (a+f+b>=1) * 0 + (a+f+b<1) * (1-a-f-b);
            elseif r2 >= r4
                c=0; e=0; b=r3; d = (a+f+b>=1) * 0 + (a+f+b<1) * (1-a-f-b);
            else
                if r2 > r3
                    c=0; b=r3; d=r2-r3; e=r4-r2;
                else
                    d=0; b=r2; c=r3-r2; e=r4-r3;
                end
            end
        end
        X=[a;b;c;d;e;f];
    end
end
