function [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
%%
% GET_HISTONE1  Targeted extraction of RT and integrated MS1 area for one histone peptide across all charge states.
%
% Syntax:
%   [cur_rts, cur_intens, cur_mono_isointens] = get_histone1( ...
%       MS1_index, MS1_peaks, ptol, unitdiff, His, hno)
%
% Purpose (overview):
%   For the peptide indicated by index 'hno' in the 'His' structure, this function:
%     • Builds MS1 extracted-ion chromatograms (XICs) around four isotopic positions
%       {M−1, M, M+1, M+2} for each available charge state of the peptide.
%     • Restricts the scan range to a narrow ±2 min window around the peptide’s reference RT.
%     • Detects candidate apex scans using GetTopBottom11 on the **monoisotopic** XIC.
%     • Selects the candidate nearest to the reference RT (±0.5 min gate, then max intensity).
%     • Integrates the chromatographic area at the chosen apex via get_area (with potential overlap handling).
%     • Returns the apex RTs and areas for all charge states; also outputs the full mono XIC for jno==1.
%
% Inputs:
%   MS1_index : [nScans × ≥3] indexer. Convention used here:
%                 - MS1_index(k,2) = retention time (RT) of scan k
%                 - MS1_index(k,3) = starting row in MS1_peaks for scan k
%   MS1_peaks : [nRows × 2] double matrix with columns [m/z, intensity] for all MS1 centroids.
%   ptol      : scalar (ppm tolerance). Historical convention: if ptol==100, it is coerced to 10.
%   unitdiff  : scalar (Da). Isotopic mass difference per 13C step (e.g., ~1.003355 Da).
%   His       : struct with at least:
%                 - pep_mz  [nPeptides × nCharge]: target m/z per charge state
%                 - pep_ch  [nPeptides × nCharge]: charge (z) per state
%                 - rt_ref  [nPeptides × 1]     : reference RT per peptide
%   hno       : scalar index of the peptide within His.
%
% Outputs:
%   cur_rts            : [1 × nCharge] apex RT per charge state (or ref RT if none found).
%   cur_intens         : [1 × nCharge] integrated peak area per charge (0 if not found).
%   cur_mono_isointens : (column vector) monoisotopic XIC used for charge state jno==1.
%
% Dependencies (external to this file):
%   - GetProfiles(...)     → builds XICs for specified m/z seeds within a scan range.
%   - GetTopBottom11(xic)  → detects candidate apex scans, returns (nt, nb, top1_idx, inten_sum).
%   - get_area(...)        → performs overlap-aware area integration at a chosen apex (via judgeOverlap1).
%
% Notes:
%   - The scan window is **fixed** at ±2 minutes around His.rt_ref(hno) for every charge state.
%   - Mono XIC is taken from column 2 of the multi-isotope XIC matrix (consistent with the codebase).
%   - If there are ≥2 qualifying candidates near the reference RT, 'nb' is tightened to bracket them.
%   - No changes to inputs/outputs or control flow have been made; only documentation added.

[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([1,ncharge]);
cur_intens = zeros([1,ncharge]);

% --------------------------------------------------------------------------------
% Define a narrow scan range around the peptide reference RT: [ref−2, ref+2].
% This focuses the search and avoids spurious distant peaks.
% --------------------------------------------------------------------------------
delta = 2;
p = find( MS1_index(:,2)>=His.rt_ref(hno)-delta );
rt_i1 = p(1);
pp = find( MS1_index(:,2)<=His.rt_ref(hno)+delta );
rt_i2 = pp(end);

% Historical convention: treat ptol==100 as 10 ppm (tighter extraction).
if ptol==100
    ptol = 10;
end;

for jno=1:ncharge
    % --------------------------------------------------------------------
    % Build MS1 profiles for current charge state:
    %   Seeds are four isotopic positions centered on the target m/z:
    %     {M−1, M, M+1, M+2} with spacing unitdiff / charge.
    %   Optionally set nC13=1 when tolerance is coarse (>100 ppm) and z≥3.
    % --------------------------------------------------------------------
    c_mz = His.pep_mz(hno,jno);
    c_ch = His.pep_ch(hno,jno);
    c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
    if ptol>100 && c_ch>=3
        nC13 = 1;
    else
        nC13 = 0;
    end;
    [c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
    j = 2;                                % Column 2 is the monoisotopic XIC by convention.
    c_mono_isointens = c_ref_isointens(:,j);
    if 1==jno
        % Expose the full monoisotopic XIC only for the first charge state.
        cur_mono_isointens = c_mono_isointens;
    end;

    % --------------------------------------------------------------------
    % Candidate peak picking on the mono XIC:
    %   - GetTopBottom11 returns candidate indices 'nt', segment bounds 'nb',
    %     a primary index 'top1_idx' and a per-candidate integrated measure 'inten_sum'.
    %   - We accept candidates within ±0.5 min of the reference RT, and among them,
   %     pick the one with the highest 'inten_sum'.
    % --------------------------------------------------------------------
    [nt,nb,top1_idx,inten_sum] = GetTopBottom11(c_mono_isointens);%#ok
    flag = zeros([1,length(nt)]);
    for i=1:length(nt)
        if abs( c_isorts(nt(i))-His.rt_ref(hno) )<=0.5
            flag(i) = 1;
        end;
    end;
    x = find(flag==1);
    if 0==isempty(x)
        %[tmp,id] = min(abs(c_isorts(nt(x))-His.rt_ref(hno)));%#ok
        [tmp,id] = max(inten_sum(x));%#ok
        top1_idx = x(id);
        cur_pos = nt(top1_idx);
        cur_rts(jno) = c_isorts(cur_pos);

        % If multiple candidates near ref RT, tighten 'nb' to bracket them:
        if length(x)>=2
            nb = [nb(x(1)) nb(x(end)+1)];
        end;

        % Integrate peak area at the chosen apex using overlap-aware get_area:
        cur_intens(jno) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
    else
        % No acceptable candidate near the reference → fallback RT and zero area.
        cur_rts(jno) = His.rt_ref(hno);
        cur_intens(jno) = 0;
    end;
end;
