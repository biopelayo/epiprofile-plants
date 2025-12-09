function His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His)
% =========================================================================
% relocateD
% -------------------------------------------------------------------------
% PURPOSE (plain English, non-proteomics audience):
%   "Re-center" (update/refine) the reference retention time (RT) of each
%   histone peptide entry using MS1 signal around its current RT estimate.
%
%   For every entry hno in His.rt_ref (starting from the 2nd element), the
%   function opens a small retention-time window [rt_ref±delta] in MS1 data
%   and calls a helper finder (get_rts) that returns candidate RTs and the
%   top (best) RT peak (top1_rt2). If no peak is found, the reference RT is
%   invalidated (set to 0). Otherwise, it is updated to top1_rt2.
%
% WHAT THIS FUNCTION CHANGES:
%   - It updates the field His.rt_ref in place (side effect) and returns
%     the modified struct.
%
% DEPENDENCIES:
%   - A helper function named GET_RTS (not included here) with signature:
%       [rts2, top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff, ...
%                                  His, hno, nsplit, t1, t2)
%     * This helper must scan MS1 within [t1,t2] and identify candidate RTs.
%     * 'top1_rt2' is expected to be the single best RT for entry 'hno'.
%
% INPUTS:
%   MS1_index : index/meta for MS1 scans (e.g., scan→RT map). Exact schema
%               depends on the upstream pipeline; it is passed to get_rts.
%   MS1_peaks : MS1 signal container (profile/centroided peaks, intensities
%               and m/z). Exact structure is pipeline-specific; used by get_rts.
%   ptol      : mass tolerance (ppm or similar) used when matching m/z in
%               MS1. IMPORTANT: if ptol == 100, it is coerced to 10 (legacy
%               convention in this codebase).
%   unitdiff  : isotopic unit mass difference per charge (or similar), also
%               passed to get_rts to build m/z search windows.
%   His       : struct holding histone-peptide metadata. Must contain:
%                 - rt_ref : vector of reference RTs (minutes). This
%                            function reads and overwrites entries 2..end.
%
% OUTPUT:
%   His       : same struct, but with rt_ref(hno) refined or set to 0 if
%               nothing was found around the expected location.
%
% KEY PARAMETERS INSIDE:
%   delta     : ±window size in minutes around the current rt_ref to search
%               for a better MS1 apex. Fixed to 1 here.
%
% SAFETY / EDGE CASES:
%   - If get_rts returns empty (no peaks), the corresponding rt_ref is set
%     to 0, explicitly marking it as "invalid / not found".
%   - Loop starts at hno = 2 (entry #1 is intentionally NOT relocated).
%
% COMPLEXITY:
%   O(N * cost(get_rts)), where N = length(His.rt_ref) - 1.
% =========================================================================

% ----------------------------
% Normalize tolerance: a legacy guardrail used elsewhere in the codebase.
% If ptol is exactly 100, force it down to 10 (narrower search).
if ptol==100
    ptol = 10;
end;

% ----------------------------
% Half-window (in minutes) around the current reference RT to search for a
% better MS1 apex. The search interval is [rt_ref - delta, rt_ref + delta].
delta = 1;

% ----------------------------
% Iterate all peptide entries EXCEPT the first one (starts at 2 by design).
for hno = 2:length(His.rt_ref)

    % Define the local RT window centered on the current reference RT.
    % t1 = left bound; t2 = right bound (in minutes).
    t1 = His.rt_ref(hno)-delta;
    t2 = His.rt_ref(hno)+delta;

    % Query MS1 in that window. nsplit=1 means "no splitting"/single slice.
    % get_rts is responsible for extracting candidate RTs and picking the
    % top candidate (top1_rt2) according to its own scoring logic.
    [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);

    % If there are no candidate RTs, invalidate this reference.
    if 1==isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        % Otherwise, accept the best candidate as the new reference RT.
        His.rt_ref(hno) = top1_rt2;
    end;
end;
