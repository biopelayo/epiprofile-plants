function [rts,top1_rt,inten_sum,top1_inten_sum] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2)
%% get_rts — Extract candidate apex RTs and integrated-intensity proxies over a user RT window.
% ----------------------------------------------------------------------------------------------
% INPUTS
%   MS1_index   : [nScans x 2] matrix; column 2 holds MS1 retention times (minutes).
%   MS1_peaks   : project-specific container with centroided MS1 peaks (format used by GetProfiles).
%   ptol        : mass tolerance setting; if exactly 100, it is internally reset to 10 (legacy rule).
%   unitdiff    : isotopic spacing per charge (≈ 1.00335 u / z), passed as a scalar.
%   His         : histone library struct with fields:
%                   - pep_mz(hno, j): theoretical m/z for peptide hno at charge j (j=1 used here)
%                   - pep_ch(hno, j): charge state for peptide hno at column j (j=1 used here)
%   hno         : peptide row index into His.pep_mz / His.pep_ch.
%   nsplit      : selector for peak finder: 1 → GetTopBottom, otherwise → GetTopBottom11.
%   t1, t2      : user-provided RT bounds (minutes) for extraction window [t1, t2].
%
% OUTPUTS
%   rts            : vector of candidate apex RTs (minutes) returned by the peak finder (subset of c_isorts).
%   top1_rt        : single RT (minutes) of the selected primary apex (by the peak finder).
%   inten_sum      : vector of area-like scores per candidate apex (aligned with rts).
%   top1_inten_sum : area-like score for the selected primary apex.
%
% NOTES
%   - This function extracts only the FIRST charge column (j=1) for the peptide library.
%   - Peak detection branch is controlled by nsplit:
%         nsplit == 1  → GetTopBottom  (typically a single-peak picker)
%         else         → GetTopBottom11 (variant allowing multi-peak segmentation)
%   - Early return if the user RT window is invalid or if the chosen top apex elutes too early (< 4 min).
% ----------------------------------------------------------------------------------------------

% Number of MS1 scans (rows in MS1_index)
num_MS1 = size(MS1_index,1);

% ---- Guard: sanity check for the requested RT window [t1, t2] ----
% If t1 is beyond the last RT, or t2 is negative, or t2 < t1, there is nothing to extract.
if t1 > MS1_index(num_MS1,2) || t2 < 0 || t2 < t1
    rts = [];
    top1_rt = [];
    inten_sum = [];
    top1_inten_sum = [];
    return;
end

% check error rt_ref fuera de [minRT,maxRT] ⇒ Causa 1/5.
% maxRT ~6000 con t1,t2 ~5–10 ⇒ unidades (Causa 2).
% n(in)=0 pero los otros paneles tienen n(in)>0 con mismo delta ⇒ Causa 3/4.
% sorted?=0 o NaNs?=1 ⇒ orden/limpieza de RT (menos probable aquí).

minRT = MS1_index(1,2); maxRT = MS1_index(end,2);
fprintf('RUN RT=[%.3f, %.3f] min | t1=%.3f t2=%.3f | rt_ref(hno)=%d: %.3f\n', ...
        minRT, maxRT, t1, t2, hno, His.rt_ref(hno));
fprintf('sorted? %d | NaNs? %d | n(>=t1)=%d | n(<=t2)=%d | n(in)= %d\n', ...
        all(diff(MS1_index(:,2))>=0), any(isnan(MS1_index(:,2))), ...
        sum(MS1_index(:,2)>=t1), sum(MS1_index(:,2)<=t2), ...
        sum(MS1_index(:,2)>=t1 & MS1_index(:,2)<=t2));
fprintf('HNO=%d | pep seq panel=%s? | rt_ref=%.3f\n', hno, 'H3_18_26 o H3_27_40', His.rt_ref(hno));

% ---- Pick peptide m/z and charge (FIRST charge column only) ----
c_mz = His.pep_mz(hno,1);  % theoretical monoisotopic m/z for peptide hno, charge column #1
c_ch = His.pep_ch(hno,1);  % corresponding charge state

% ---- Map the user RT window [t1, t2] to scan indices ----
%  p  : indices of scans with RT >= t1; choose the first as window start
%  pp : indices of scans with RT <= t2; choose the last as window end
p  = find( MS1_index(:,2) >= t1 );
rt_i1 = p(1);
pp = find( MS1_index(:,2) <= t2 );
rt_i2 = pp(end);

% ---- Build isotopologue anchors: [M−1, M, M+1, M+2] ----
c_ref_isomzs = [ c_mz - unitdiff/c_ch, c_mz, c_mz + unitdiff/c_ch, c_mz + 2*unitdiff/c_ch ];

% ---- Legacy tolerance rule ----
if ptol == 100
    ptol = 10;  % internally narrow if a "100" placeholder was passed
end

% ---- Coarse tolerance + high charge heuristic for GetProfiles ----
if ptol > 100 && c_ch >= 3
    nC13 = 1;
else
    nC13 = 0;
end

% ---- Extract time-aligned isotopic profiles within [rt_i1 : rt_i2] ----
% c_isorts        : time grid (minutes) for the selected scans
% c_ref_isointens : [nTime x 4] intensities for M−1, M, M+1, M+2 (column 2 is monoisotopic M)
[c_isorts, c_ref_isointens] = GetProfiles( ...
    MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2);

% ---- Select the monoisotopic trace (column 2) for peak detection ----
c_mono_isointens = c_ref_isointens(:,2);

% ---- Peak detection branch: choose helper according to nsplit ----
if 1 == nsplit
    % Simpler peak finder / segmentation
    [nt, nb, top1_idx, inten_sum] = GetTopBottom(c_mono_isointens); %#ok
else
    % Alternative variant (often more permissive with multiple maxima)
    [nt, nb, top1_idx, inten_sum] = GetTopBottom11(c_mono_isointens); %#ok
end

% Map candidate indices nt to RT values on the time grid:
rts = c_isorts(nt);

% Primary apex RT (selected by the helper) and its area-like score:
top1_rt         = c_isorts( nt(top1_idx) );
top1_inten_sum  = inten_sum( top1_idx );

% ---- Early-elution guard: discard peaks eluting before 4 min ----
% [Inference] This threshold likely removes solvent front / void-volume artifacts.
if top1_rt < 4
    rts = [];
    top1_rt = [];
    inten_sum = [];
    top1_inten_sum = [];
    return;
end
