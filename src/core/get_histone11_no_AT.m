function [cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
%% get_histone11 — Extracts MS1 monoisotopic profiles around a reference RT for a given histone peptide,
%                  selects the apex near the reference, and integrates the peak area per charge state.
% -----------------------------------------------------------------------------------------------------------------
% INPUTS
%   MS1_index   : [nScans x 2] matrix. Column 1 = scan number (or index), Column 2 = retention time (minutes).
%   MS1_peaks   : Project-specific structure/array with centroided MS1 peaks (see GetProfiles for expected layout).
%   ptol        : Mass tolerance setting (coarse/fine threshold used by GetProfiles; special handling for ITMS below).
%   unitdiff    : Mass difference between isotopologues per charge unit (typically ~1.00335 u for 13C; divided by z).
%   His         : Structure with histone peptide metadata:
%                   - His.pep_mz(hno, j)  : theoretical m/z for peptide hno at charge state j
%                   - His.pep_ch(hno, j)  : charge state (z) for peptide hno at column j
%                   - His.rt_ref(hno)     : reference RT (minutes) for peptide hno
%   hno         : Peptide index in His.* matrices (row selector).
%
% OUTPUTS
%   cur_rts            : [1 x ncharge] array of selected apex RTs per charge (minutes).
%   cur_intens         : [1 x ncharge] array of integrated areas per charge (computed via get_area).
%   cur_mono_isointens : [nTimepoints x 1] monoisotopic intensity trace (time-series) for the FIRST charge only.
%                        (Returned as-is from GetProfiles' monoisotopic column; useful for diagnostics/plotting.)
%
% NOTES
%   - The function scans a narrow RT window around His.rt_ref(hno) (±0.5 min), extracts isotopic traces centered at
%     [M−1, M, M+1, M+2] m/z positions (see c_ref_isomzs), and uses GetTopBottom11 to locate candidate apex indices.
%   - If multiple apex candidates lie within ±0.5 min of the reference RT, the one with the largest integrated signal
%     (inten_sum) is chosen for robust quantification.
%   - For ITMS-like coarse tolerance (ptol > 100) and high charge (z >= 3), a flag (nC13 = 1) is passed to GetProfiles
%     to adjust isotopic extraction behavior (implementation detail inside GetProfiles).
%   - If no valid apex is found near the reference RT, the function returns RT = His.rt_ref(hno) and area = 0 for that
%     charge state, to keep downstream matrices well-formed.
% -----------------------------------------------------------------------------------------------------------------

% Determine peptide/charge matrix dimensions from the histone library:
[npep,ncharge] = size(His.pep_mz); %#ok  % npep is unused here; ncharge defines the number of charge states (columns).

% Pre-allocate outputs for all charge states of peptide hno:
cur_rts    = zeros([1,ncharge]);  % Selected apex RT per charge.
cur_intens = zeros([1,ncharge]);  % Integrated area per charge.

% Define a half-window (in minutes) around the reference RT for local extraction:
delta = 0.5;

% Find the first MS1 scan index at or after (rt_ref - delta):
p    = find( MS1_index(:,2) >= His.rt_ref(hno) - delta );
rt_i1 = p(1);  % Start scan index of the local RT window.

% Find the last MS1 scan index at or before (rt_ref + delta):
pp    = find( MS1_index(:,2) <= His.rt_ref(hno) + delta );
rt_i2 = pp(end);  % End scan index of the local RT window.

% Special-case: if ptol equals 100, internally narrow it to 10.
% [Inference] This maps a legacy "100" setting to a more stringent extraction tolerance.
if ptol == 100
    ptol = 10;
end

% Iterate over all possible charge states (columns) for peptide hno:
for jno = 1:ncharge

    % ----- Build theoretical m/z anchors for isotopic profile extraction -----
    % Get theoretical precursor m/z and charge for this peptide/charge column:
    c_mz = His.pep_mz(hno, jno);  % Theoretical m/z at charge state jno.
    c_ch = His.pep_ch(hno, jno);  % Charge state (z) at column jno.

    % Construct reference isotopologue m/z array:
    % [M−1, M, M+1, M+2] spaced by (unitdiff / z).
    % Note: The "M−1" slot is included here as a symmetric reference (useful for baseline/noise checks).
    c_ref_isomzs = [ c_mz - unitdiff/c_ch, c_mz, c_mz + unitdiff/c_ch, c_mz + 2*unitdiff/c_ch ];

    % Toggle isotopic behavior for coarse tolerance + high charge:
    % If extraction tolerance is coarse (ptol > 100) and z >= 3, enable 'nC13' flag.
    % [Inference] Inside GetProfiles, nC13 may alter how 13C envelopes are matched/clustered.
    if ptol > 100 && c_ch >= 3
        nC13 = 1;
    else
        nC13 = 0;
    end

    % ----- Extract time-aligned isotopic intensity traces within the local RT window -----
    % c_isorts           : vector of RTs (minutes) for scans rt_i1:rt_i2.
    % c_ref_isointens    : [nTimepoints x 4] matrix of intensities for [M−1, M, M+1, M+2], in this order.
    [c_isorts, c_ref_isointens] = GetProfiles(...
        MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2);

    % By convention in this codebase, column j=2 corresponds to monoisotopic "M".
    j = 2;
    c_mono_isointens = c_ref_isointens(:, j);

    % For diagnostics, store the monoisotopic trace for the first charge only:
    if jno == 1
        cur_mono_isointens = c_mono_isointens;
    end

    % ----- Locate candidate apex positions and integration bounds -----
    % GetTopBottom11 returns:
    %   nt        : indices (into c_isorts) of local apex candidates
    %   nb        : boundary indices (peak start/end per segment)
    %   top1_idx  : (not used directly here) primary apex candidate
    %   inten_sum : integrated intensity per candidate region (area-like score)
    [nt, nb, top1_idx, inten_sum] = GetTopBottom11(c_mono_isointens); %#ok

    % Build a mask for apex candidates that lie within ±0.5 min of the reference RT:
    flag = zeros([1, length(nt)]);
    for i = 1:length(nt)
        if abs( c_isorts(nt(i)) - His.rt_ref(hno) ) <= 0.5
            flag(i) = 1;
        end
    end
    x = find(flag == 1);  % Indices into the candidate list 'nt' that satisfy the RT proximity criterion.

    if ~isempty(x)
        % If multiple valid candidates exist near the reference RT, pick the one with largest area proxy (inten_sum).
        % (Alternative, commented in original, would be the closest in RT; here we prefer intensity-driven robustness.)
        % [~, id] = min(abs(c_isorts(nt(x)) - His.rt_ref(hno))); %#ok  % (previous strategy by RT proximity)
        [~, id] = max(inten_sum(x)); %#ok  % Select the most abundant candidate among those near the reference.

        top1_idx = x(id);           % Winner index in the candidate list.
        cur_pos  = nt(top1_idx);    % Position in c_isorts (scan index within the local window).

        % Report apex RT for this charge:
        cur_rts(jno) = c_isorts(cur_pos);

        % If there are ≥2 candidates in the window, coalesce integration bounds to cover from first to last:
        if length(x) >= 2
            % nb holds boundary indices per candidate; here we form a combined [left, right] span.
            nb = [ nb(x(1)), nb(x(end) + 1) ];
        end

        % Integrate area for this charge using the (possibly coalesced) bounds and apex position:
        cur_intens(jno) = get_area(...
            c_isorts, c_ref_isointens, nb, cur_pos, ...
            c_mz, c_ch, MS1_index, MS1_peaks, unitdiff, ptol);

    else
        % No valid apex near the reference RT — fall back to reference RT and area 0 for this charge.
        cur_rts(jno)    = His.rt_ref(hno);
        cur_intens(jno) = 0;
    end
end
