function [cur_rts,cur_intens,cur_mono_isointens] = get_histone12(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
%% get_histone12 — Extracts MS1 monoisotopic profiles in a *tighter* RT window (±0.2 min)
%                  around a peptide's reference RT, selects an apex near rt_ref, and integrates area.
% ----------------------------------------------------------------------------------------------------
% INPUTS
%   MS1_index   : [nScans x 2] matrix; column 2 holds retention times (minutes).
%   MS1_peaks   : structure/array of centroided MS1 peaks (layout expected by GetProfiles).
%   ptol        : mass tolerance parameter; if exactly 100, it is internally reset to 10.
%   unitdiff    : isotopic mass spacing per charge (≈1.00335 u / z, passed in as unitdiff).
%   His         : histone library with fields:
%                   - pep_mz(hno, j) : theoretical m/z for peptide hno at charge state j
%                   - pep_ch(hno, j) : charge state for peptide hno at column j
%                   - rt_ref(hno)    : reference retention time (minutes)
%   hno         : row index selecting a peptide in His.pep_mz/pep_ch and rt_ref.
%
% OUTPUTS
%   cur_rts            : [1 x ncharge] selected apex RTs per charge state.
%   cur_intens         : [1 x ncharge] integrated areas per charge state.
%   cur_mono_isointens : monoisotopic intensity time series (for the FIRST charge only).
%
% NOTE
%   Compared to get_histone11, this function uses a *narrower* extraction window (delta=0.2)
%   but still filters candidates within ±0.5 min of rt_ref (the latter becomes non-restrictive
%   because the extracted RT grid already lies within ±0.2). Peak detection helper differs:
%   here it calls GetTopBottom (not GetTopBottom11).
% ----------------------------------------------------------------------------------------------------

[npep,ncharge] = size(His.pep_mz); %#ok  % Determine number of peptides (rows) and charge states (columns).
cur_rts = zeros([1,ncharge]);           % Preallocate output for apex RT per charge.
cur_intens = zeros([1,ncharge]);        % Preallocate output for integrated area per charge.

delta = 0.2;                            % Half-window (minutes) around rt_ref for MS1 extraction (tighter than 0.5).

% Find the first scan with RT >= (rt_ref - delta):
p = find( MS1_index(:,2) >= His.rt_ref(hno) - delta );
rt_i1 = p(1);                           % Start index of the local RT window.

% Find the last scan with RT <= (rt_ref + delta):
pp = find( MS1_index(:,2) <= His.rt_ref(hno) + delta );
rt_i2 = pp(end);                        % End index of the local RT window.

% Legacy adjustment: if ptol equals 100, internally reduce it to 10 (stricter tolerance).
if ptol==100
    ptol = 10;
end

% Iterate over all possible charge states for this peptide:
for jno=1:ncharge
    % ----- Build theoretical m/z anchors for [M−1, M, M+1, M+2] -----
    c_mz = His.pep_mz(hno,jno);        % Theoretical monoisotopic m/z at this charge column.
    c_ch = His.pep_ch(hno,jno);        % Charge state (z) for this column.

    % Isotopologue positions spaced by unitdiff/z:
    c_ref_isomzs = [c_mz-unitdiff/c_ch, c_mz, c_mz+unitdiff/c_ch, c_mz+2*unitdiff/c_ch];

    % Heuristic for coarse tolerance + high charge: signal envelopes treated differently by GetProfiles.
    if ptol>100 && c_ch>=3
        nC13 = 1;
    else
        nC13 = 0;
    end

    % ----- Extract time-aligned isotopic intensity traces within [rt_i1:rt_i2] -----
    [c_isorts,c_ref_isointens] = GetProfiles( ...
        MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2);

    % Select the monoisotopic trace (column j=2 corresponds to M):
    j = 2;
    c_mono_isointens = c_ref_isointens(:,j);

    % Keep the monoisotopic time series for the first charge only (for diagnostics/plotting):
    if 1==jno
        cur_mono_isointens = c_mono_isointens;
    end

    % ----- Locate candidate apex positions and boundaries using GetTopBottom (variant) -----
    [nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens); %#ok

    % Build a mask of candidate tops that lie within ±0.5 min of the reference RT:
    flag = zeros([1,length(nt)]);
    for i=1:length(nt)
        if abs( c_isorts(nt(i)) - His.rt_ref(hno) ) <= 0.5
            flag(i) = 1;
        end
    end
    x = find(flag==1);                  % Candidate indices that pass the RT proximity filter.

    if 0==isempty(x)
        % Among candidates near the reference RT, pick the one with maximum area proxy:
        % (Alternative "closest in RT" is commented in the original code.)
        [tmp,id] = max(inten_sum(x)); %#ok
        top1_idx = x(id);
        cur_pos = nt(top1_idx);

        % Report the apex RT for this charge:
        cur_rts(jno) = c_isorts(cur_pos);

        % If multiple candidates are present, coalesce integration bounds from first to last:
        if length(x)>=2
            nb = [nb(x(1)) nb(x(end)+1)];
        end

        % Integrate area using get_area with chosen bounds and apex position:
        cur_intens(jno) = get_area( ...
            c_isorts, c_ref_isointens, nb, cur_pos, ...
            c_mz, c_ch, MS1_index, MS1_peaks, unitdiff, ptol);

    else
        % No valid candidate near the reference RT — fall back to rt_ref and area = 0.
        cur_rts(jno) = His.rt_ref(hno);
        cur_intens(jno) = 0;
    end
end
