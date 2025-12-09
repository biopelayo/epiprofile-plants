function [cur_rts,cur_intens,cur_mono_isointens] = get_histone13(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
%%
% GET_HISTONE13  Ultra-narrow RT-anchored extraction and integration for a histone peptide across charge states.
%
% Synopsis:
%   [cur_rts, cur_intens, cur_mono_isointens] = get_histone13(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
%
% High-level purpose:
%   Quantify the MS1 chromatographic peak area for a specific histone peptide (row 'hno' in 'His')
%   across ALL its charge states, using a *very tight* scan window around its reference RT:
%     • Extraction window for scans: ±0.1 min around His.rt_ref(hno) (ultra-focused search).
%     • Candidate acceptance gate: ±0.5 min around His.rt_ref(hno) (slightly looser for safety).
%   Peak area is integrated via GET_AREA (with isotopic overlap handling when applicable).
%
% Inputs:
%   MS1_index : [nScans × ≥3] numeric. Conventions used here:
%                 - MS1_index(:,2) = retention time (minutes) for each MS1 scan
%                 - MS1_index(:,3) = starting row index in MS1_peaks for that scan
%   MS1_peaks : [nRows × 2] numeric. Concatenated MS1 centroids as [m/z, intensity].
%   ptol      : scalar (ppm). Historical rule: if ptol==100, it is coerced to 10 (tighter).
%   unitdiff  : scalar (Da). Isotopic spacing per 13C step (≈1.003355 Da), used as unitdiff/charge.
%   His       : struct with at least:
%                 • pep_mz  [nPeptides × nCharge] : target m/z per charge
%                 • pep_ch  [nPeptides × nCharge] : charge state per column
%                 • rt_ref  [nPeptides × 1]       : reference RT per peptide
%   hno       : scalar. Peptide row index in 'His'.
%
% Outputs:
%   cur_rts            : [1 × nCharge] apex RT per charge (fallback = His.rt_ref(hno) if none).
%   cur_intens         : [1 × nCharge] integrated MS1 area per charge (0 if no accepted candidate).
%   cur_mono_isointens : column vector. The full monoisotopic XIC used for the *first* charge state.
%
% Key traits / differences vs. get_histone1/get_histone10:
%   • Uses an *ultra-narrow* scan extraction window (±0.1 min), ideal when RT reference is highly reliable.
%   • Candidate gate remains ±0.5 min to allow slight deviations but extraction itself is constrained to ±0.1.
%   • Candidate selection uses GetTopBottom (not GetTopBottom11), then chooses the **max inten_sum** within gate.
%
% External dependencies:
%   GetProfiles(...)  → builds XICs around seeds {M−1, M, M+1, M+2}
%   GetTopBottom(xic) → detects candidate apex indices (nt), segment boundaries (nb), and intensity sums
%   get_area(...)     → overlap-aware peak area integration (may call judgeOverlap1 internally)
%

[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([1,ncharge]);
cur_intens = zeros([1,ncharge]);

delta = 0.1;                                   % Ultra-tight extraction window around rt_ref (±0.1 min).
p = find( MS1_index(:,2)>=His.rt_ref(hno)-delta );
rt_i1 = p(1);
pp = find( MS1_index(:,2)<=His.rt_ref(hno)+delta );
rt_i2 = pp(end);

if ptol==100
    ptol = 10;                                  % Historical convention: coerce 100 → 10 ppm for extraction.
end;

for jno=1:ncharge
    % ------------------------------------------------------------
    % Build MS1 profiles (XICs) for current charge state 'jno'.
    % Seeds: four isotope lines {M−1, M, M+1, M+2} spaced by unitdiff/charge.
    % nC13 hint is enabled for coarse ppm (>100) and z≥3, as elsewhere in codebase.
    % ------------------------------------------------------------
    c_mz = His.pep_mz(hno,jno);
    c_ch = His.pep_ch(hno,jno);
    c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
    if ptol>100 && c_ch>=3
        nC13 = 1;
    else
        nC13 = 0;
    end;
    [c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
    j = 2;                                      % By convention, column 2 is the monoisotopic XIC.
    c_mono_isointens = c_ref_isointens(:,j);
    if 1==jno
        cur_mono_isointens = c_mono_isointens;  % Expose mono XIC only for the first charge state.
    end;

    % ------------------------------------------------------------
    % Candidate picking on the monoisotopic XIC using GetTopBottom.
    % Acceptance gate: |RT − His.rt_ref(hno)| ≤ 0.5 min (wider than extraction window).
    % Selection: choose the candidate with maximum inten_sum among accepted ones.
    % If ≥2 candidates are accepted, tighten 'nb' to bracket them before integration.
    % ------------------------------------------------------------
    [nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens);%#ok
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
        if length(x)>=2
            nb = [nb(x(1)) nb(x(end)+1)];       % Tighten integration segment to bracket accepted candidates.
        end;
        cur_intens(jno) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
    else
        cur_rts(jno) = His.rt_ref(hno);         % Fallback: set RT to reference and area to zero if none accepted.
        cur_intens(jno) = 0;
    end;
end;
