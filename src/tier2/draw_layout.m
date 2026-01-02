function draw_layout(cur_outpath, out_filename, His, pep_rts, pep_intens, isorts, mono_isointens, MS2_index, MS2_peaks, special)
%%
% ========================================================================
% draw_layout  —  Produce a multi-panel PDF showing MS1 XICs per peptide
%                 (histone modified forms) and optionally overlays MS2
%                 fragment evidence (DIA-like). Marks local peaks and
%                 integration boundaries for each row/variant.
%
% Audience note (plain English, non-proteomics):
%   - MS1 is a time series per mass/charge (m/z) tracking the abundance of
%     a precursor ion while it elutes in chromatography (XIC = intensity vs time).
%   - MS2 are fragmentation spectra acquired around certain times; their
%     fragment ions (b/y for CID, c/z for ETD) help confirm identity.
%   - We draw one subplot per peptide “row” (modified form), showing its MS1
%     trace, a labeled apex, vertical lines delimiting the extraction window,
%     and (optionally) small overlay traces indicating MS2 fragments.
%
% This version includes additional safeguards compared to a prior variant:
%   - Defensive checks when indexing the time window (p1/p2 possibly empty).
%   - Seed position for MS1 can fall back to the first time index.
%   - In the nested MatchMS2(), we scan pno = cpos-1:cpos+1 (±1 neighbor)
%     instead of only pno=cpos, and re-evaluate activation/tolerance if
%     instruments vary across scans. See README “Differences” section.
% ========================================================================

% -----------------------------
% 1) Filter peptide rows to plot:
%    Condition: (a) RT > 4 min (to avoid early noise) AND
%               (b) His.display == 1 (pre-approved rows)
npep = size(His.pep_mz, 1);
ix = find(pep_rts(1:npep,1) > 4 & His.display(1:npep) == 1);
if isempty(ix)
    % No valid rows to plot → nothing to draw for this window.
    return;
end

% -----------------------------
% 2) For each selected row, compute a seed position on the time axis,
%    try to pick a local apex and an integration window ("terminus").
nplot = length(ix);
localmax_rt    = zeros(nplot, 1);  % apex time to display (locked to pep_rts)
localmax_inten = zeros(nplot, 1);  % apex intensity (from local search or 0)
terminus       = zeros(nplot, 2);  % [startIdx, endIdx] in isorts (time index)

for ino = 1:nplot
    cno = ix(ino);

    % 2a) Find the last time index not exceeding the row RT (seed for search).
    p = find(isorts <= pep_rts(cno,1));
    if isempty(p)
        % If RT lies below the minimum time in isorts, defensively seed at 1.
        warning('draw_layout: pep_rts(%d,1)=%.4g < min(isorts). Forcing c_ms1pos = 1.', cno, pep_rts(cno,1));
        c_ms1pos = 1;
    else
        c_ms1pos = p(end);
    end

    % 2b) Extract the mono-isotopic MS1 trace for this row.
    c_mono_isointens = mono_isointens(:, cno);

    if pep_intens(cno,1) > 0
        % 2c) Try to define a plausible local peak (and its boundaries).
        [nt, nb] = GetTopBottom11(c_mono_isointens); %#ok (used internally by GetLocal)
        [localmax_rt(ino), localmax_inten(ino), IX] = GetLocal(c_ms1pos, isorts, c_mono_isointens, nb);

        if isempty(IX)
            % No robust subrange identified → degenerate window at the seed.
            warning('draw_layout: empty IX for cno=%d. Forcing terminus = [seed seed].', cno);
            term_start = c_ms1pos;
            term_end   = c_ms1pos;
        else
            term_start = IX(1);
            term_end   = IX(end);
        end
        terminus(ino,1:2) = [term_start, term_end];

        % IMPORTANT: The displayed apex time is locked to the nominal RT.
        localmax_rt(ino) = pep_rts(cno,1);
    else
        % 2d) If the integrated scalar intensity is zero, skip local search.
        term_start           = c_ms1pos;
        term_end             = c_ms1pos;
        terminus(ino,1:2)   = [term_start, term_end];
        localmax_rt(ino)    = pep_rts(cno,1);
        localmax_inten(ino) = 0;
    end
end  % for ino

% -----------------------------
% 3) Compute global X-limits for plotting:
%    Take the earliest/largest boundary across rows and extend for context.
st = isorts(min(terminus(:,1)));  % leftmost boundary time
st = max(st - 5, 1);              % pad to the left (at least 1)
tm = isorts(max(terminus(:,2))) + 25; % pad to the right

% -----------------------------
% 4) Prepare an offscreen figure and output PDF path.
out_file1 = fullfile(fileparts(cur_outpath), [out_filename, '.pdf']);
warning off all;
set(gcf, 'visible', 'off');  % run headless (no window popping)

% Normalize filename if it starts with 'HH' (historical prefix).
if strcmp(out_filename(1:2), 'HH')
    out_filename = out_filename(2:end);
end

% Title: parse pieces from out_filename + peptide sequence + charge.
p = strfind(out_filename, '_');
cur_title = [ ...
    out_filename(1:p(1)-1), ' ', ...
    out_filename(p(2)+1:p(3)-1), '-', ...
    out_filename(p(3)+1:end), ' ', ...
    His.pep_seq, ' +', num2str(His.pep_ch(1,1)), ' ions' ...
];

% Mode flags for chemistry/instrument logic downstream.
nhmass = special.nhmass;

% Retrieve modification catalog (used in MS2 matching).
Mods   = GetMods();

% Color palette for fragment overlays (wraps every 6 entries).
colors = {'k','r','g','b','c','m'};

% -----------------------------
% 5) Draw one subplot per selected row:
for ino = 1:nplot
    cno = ix(ino);

    % Subplot wrapper:
    subplot(nplot, 1, ino);

    % 5b) Draw the MS1 mono-isotopic XIC (blue line). Hide ticks for compactness.
    plot(isorts, mono_isointens(:, cno), 'color', 'b', 'linewidth', 1);
    set(gca, 'xtick', [], 'ytick', []);
    hold on;

    % 5c) Set X-limits for this panel:
    xlim([st tm]);

    % 5d) Compute safe in-window indices [p1(end) .. p2(end)]:
    p1 = find(isorts <= st);
    p2 = find(isorts <= tm);

    if isempty(p1)
        warning('draw_layout: p1 empty (st=%.4g < min time). Forcing p1 = 1.', st);
        p1 = 1;
    end
    if isempty(p2)
        warning('draw_layout: p2 empty (tm=%.4g > max time). Forcing p2 = numel(isorts).', tm);
        p2 = numel(isorts);
    end
    IX = p1(end) : p2(end);

    % 5h) Choose Y-limits from the local maximum inside [st, tm].
    tmp_maxinten = max(mono_isointens(IX, cno));
    if tmp_maxinten > 0
        ylim([0 1.05 * tmp_maxinten]);
    end

    % 5i) Mark the apex (magenta dot/line) and annotate with a red label.
    plot(localmax_rt(ino), localmax_inten(ino), 'color', 'm', 'linestyle', '-', 'linewidth', 1);
    cur_txt = [ ...
        His.mod_short{cno}, '(', ...
        num2str(His.pep_mz(cno,1), '%.4f'), ',+', ...
        num2str(His.pep_ch(cno,1)), '), ', ...
        num2str(localmax_rt(ino), '%.2f'), ', ', ...
        num2str(localmax_inten(ino), '%.2e') ...
    ];
    if pep_intens(cno,1) > 0
        text(localmax_rt(ino), 1.05 * tmp_maxinten, cur_txt, 'color', 'r', 'fontsize', 7);
    else
        text(localmax_rt(ino), localmax_inten(ino), cur_txt, 'color', 'r', 'fontsize', 7);
    end

    % 5j) Draw vertical boundary lines at the integration limits.
    plot([isorts(terminus(ino,1)), isorts(terminus(ino,1))], [0 1], 'color', 'm', 'linestyle', '-', 'linewidth', 1);
    plot([isorts(terminus(ino,2)), isorts(terminus(ino,2))], [0 1], 'color', 'm', 'linestyle', '-', 'linewidth', 1);

    % Title on the first subplot only (to avoid repetition).
    if ino == 1
        title(cur_title);
    end

    % 5k) Optional MS2 overlay (DIA mode): draw small fragment traces/labels.
    if special.nDAmode == 2
        % Define a narrow MS2 time window aligned to the integration limits:
        rt1 = isorts(terminus(ino,1)) - 0.001;
        rt2 = isorts(terminus(ino,2)) + 0.001;

        % Run MS2 matching in that window for this row/precursor:
        [ms2pos, ms2rts, ms2intens, posn, posc, ActiveType] = MatchMS2( ...
            MS2_index, MS2_peaks, Mods, His, cno, rt1, rt2, nhmass);

        if isempty(ms2pos)
            % Nothing to overlay in this time window.
            continue;
        end

        % Label sets depending on activation type (CID→b/y, ETD→c/z).
        if strcmp(ActiveType, 'CID')
            strn = 'b';
            strc = 'y';
        else
            strn = 'c';
            strc = 'z';
        end

        % Scale MS2 intensities to fit under the MS1 panel height:
        new_maxinten = max(ms2intens(:));
        if new_maxinten > 0
            fold = (tmp_maxinten / new_maxinten) / (1 + length(posn));
        else
            fold = 1 / (1 + length(posn));
        end

        % Draw N-series (b or c) slightly left-shifted (−4 minutes).
        for kno = 1:length(posn)
            plot( ...
                ms2rts(ms2pos) - 4, ...
                fold * ms2intens(ms2pos, kno) + (kno - 1) * fold * new_maxinten, ...
                'color', colors{mod(kno, 6) + 1}, 'linestyle', '-', 'linewidth', 0.5 ...
            );
            text( ...
                ms2rts(ms2pos(end)) - 4, ...
                fold * ms2intens(ms2pos(end), kno) + (kno - 1) * fold * new_maxinten, ...
                [strn, num2str(posn(kno))], ...
                'color', colors{mod(kno, 6) + 1}, 'fontsize', 7 ...
            );
        end

        % Draw C-series (y or z) slightly right-shifted (+3 minutes).
        for kno = 1:length(posc)
            qno = kno + length(posn);
            plot( ...
                ms2rts(ms2pos) + 3, ...
                fold * ms2intens(ms2pos, qno) + (kno - 1) * fold * new_maxinten, ...
                'color', colors{mod(kno, 6) + 1}, 'linestyle', '-', 'linewidth', 0.5 ...
            );
            text( ...
                ms2rts(ms2pos(end)) + 3, ...
                fold * ms2intens(ms2pos(end), qno) + (kno - 1) * fold * new_maxinten, ...
                [strc, num2str(posc(kno))], ...
                'color', colors{mod(kno, 6) + 1}, 'fontsize', 7 ...
            );
        end
    end  % if DIA
end  % for ino

% -----------------------------
% 6) Finalize axes, label, write PDF, and close figure.
set(gca, 'xtickMode', 'auto');
xlabel('Time (min)');
ylabel('Abundance');
print('-dpdf', out_file1);
close();

end  % function draw_layout


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NESTED HELPER  —  MatchMS2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ms2pos, ms2rts, ms2intens, posn, posc, ActiveType] = MatchMS2( ...
    MS2_index, MS2_peaks, Mods, His, hno, rt1, rt2, nhmass)
%%
% PURPOSE (English, non-proteomics):
%   Given a target precursor m/z (from His.pep_mz(hno,1)) and a time window
%   [rt1, rt2], select relevant MS2 scans and extract intensities for key
%   fragment ions that help confirm the peptide (b/y for CID, c/z for ETD).
%
% Inputs
%   MS2_index : table with per-scan metadata:
%       col 2: scan retention time (minutes)
%       col 4: precursor m/z
%       col 6: instrument code (1..6 → CID/ETD/HCD, IT/FT)
%       col 7: starting row index into MS2_peaks for this scan
%   MS2_peaks : [all_scans_concatenated_rows x 2] matrix of [m/z, intensity]
%   Mods, His : domain metadata (peptide sequence, charge, mods, etc.)
%   hno       : row index in this peptide window
%   rt1, rt2  : lower/upper time bounds
%   nhmass    : choose H-mass variant for fragment lists if ==1
%
% Outputs
%   ms2pos    : scan indices deemed relevant in [rt1,rt2] (exact m/z match
%               if found; otherwise all scans in range as a fallback)
%   ms2rts    : vector of RTs for all scans (used for plotting)
%   ms2intens : [num_scans x num_key_ions] matched intensities
%   posn,posc : positions for N- and C-terminal fragment series (for labels)
%   ActiveType: 'CID' or 'ETD' (drives b/y vs c/z series and tolerances)
%
% Notes on THIS version (differs from the other variant):
%   - If instruments differ across scans, we *recompute* ActiveType/tol and
%     key ions per-scan (and we iterate pno=cpos-1:cpos+1).
%   - We handle empty “instruments” defensively (set default tol/type).
% -------------------------------------------------------------------------

% 1) Choose the closest observed precursor m/z to the theoretical one.
num_MS2 = size(MS2_index, 1);
c_mz    = His.pep_mz(hno, 1);
premzs  = unique(MS2_index(:, 4));
[~, ii] = min(abs(premzs - c_mz));      %#ok (ii used to pick target)
target  = premzs(ii);

% 2) Time filter in [rt1, rt2], with early exit if no overlap.
p  = find(MS2_index(:, 2) >= rt1);
pp = find(MS2_index(:,  2) <= rt2);
if isempty(p) || isempty(pp) || p(1) > pp(end)
    ms2pos     = [];
    ms2rts     = [];
    ms2intens  = [];
    posn       = [];
    posc       = [];
    ActiveType = [];
    return;
end
i1 = p(1);
i2 = pp(end);

% 3) Prefer exact matches to the chosen target precursor m/z within [i1:i2].
flag = zeros(num_MS2, 1);
for i = i1:i2
    if MS2_index(i, 4) == target
        flag(i) = 1;
    end
end
ms2pos = find(flag == 1);
if isempty(ms2pos)
    % No exact precursor → take all scans in the window (lenient fallback).
    ms2pos = i1:i2;
end

% 4) Vector of all scan RTs (global time base for overlays).
ms2rts = MS2_index(:, 2);

% 5) Activation type (CID vs ETD) and tolerances (absolute m/z):
%    Odd codes → IT (wider tol=0.4), even codes → FT (tighter tol=0.02).
instruments = MS2_index(ms2pos, 6);
if isempty(instruments)
    ActiveType = '';
    tol        = 0.02;
else
    c_instrument = instruments(1);
    if any(c_instrument == [3, 4])
        ActiveType = 'ETD';
    else
        ActiveType = 'CID';
    end
    if mod(c_instrument, 2) == 1
        tol = 0.4;
    else
        tol = 0.02;
    end
end

% 6) Build the key fragment ion list (b/y or c/z) for the current peptide.
if nhmass == 1
    [K1, posn, posc] = get_key_ions1H(His, hno, Mods, ActiveType);
else
    [K1, posn, posc] = get_key_ions1(His, hno, Mods, ActiveType);
end

% 7) “index” points to where each scan’s peaks start inside MS2_peaks.
index = [1; MS2_index(1:num_MS2, 7)];

% 8) For each selected scan, extract matched intensities for all K1 ions.
ms2intens = zeros(num_MS2, length(K1));

for idx_i = 1:length(ms2pos)
    cpos = ms2pos(idx_i);

    % Compared to the other variant, consider neighbor scans (±1) as well:
    for pno = cpos-1 : cpos+1
        if pno < 1 || pno > num_MS2
            continue;
        end

        % If instrument types vary across scans, adapt per scan:
        if length(unique(instruments)) > 1
            c_instrument = MS2_index(pno, 6);
            if any(c_instrument == [3, 4])
                ActiveType = 'ETD';
            else
                ActiveType = 'CID';
            end
            if mod(c_instrument, 2) == 1
                tol = 0.4;
            else
                tol = 0.02;
            end
            % Rebuild fragment list for the new ActiveType.
            if nhmass == 1
                [K1, posn, posc] = get_key_ions1H(His, hno, Mods, ActiveType);
            else
                [K1, posn, posc] = get_key_ions1(His, hno, Mods, ActiveType);
            end
        end

        % Slice peaks belonging to scan pno and match against K1 with tolerance.
        IX   = index(pno) : (index(pno + 1) - 1);
        mz   = MS2_peaks(IX, 1);
        inten= MS2_peaks(IX, 2);

        for j = 1:length(K1)
            ix1 = find(abs(mz - K1(j)) <= tol);
            if ~isempty(ix1)
                [~, x1] = min(abs(mz(ix1) - K1(j)));  %#ok (nearest in m/z)
                if ms2intens(cpos, j) < inten(ix1(x1))
                    ms2intens(cpos, j) = inten(ix1(x1)); % keep the strongest match
                end
            end
        end
    end
end

end  % function MatchMS2
