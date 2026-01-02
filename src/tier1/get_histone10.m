function [cur_rts,cur_intens,cur_mono_isointens] = get_histone10(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
%%
% GET_HISTONE10  Targeted extraction of apex RT and integrated MS1 area across all charge states
%                using a *nearest-to-reference RT* selection rule (as opposed to max-intensity).
%
% Synopsis
%   [cur_rts, cur_intens, cur_mono_isointens] = get_histone10( ...
%       MS1_index, MS1_peaks, ptol, unitdiff, His, hno)
%
% High-level purpose
%   Quantify a histone peptide (row 'hno' in struct 'His') for each available charge state by:
%     1) Building MS1 XICs around four isotope positions {M−1, M, M+1, M+2}.
%     2) Restricting the scan range to a tight ±2-min window around the peptide's reference RT.
%     3) Detecting candidate apex scans on the monoisotopic XIC (col 2) with GetTopBottom11.
%     4) Selecting the candidate whose RT is *closest* to His.rt_ref(hno) (within ±0.5 min).
%     5) Integrating area at the chosen apex via get_area (with overlap handling if applicable).
%
% Key difference vs. get_histone1
%   • Candidate selection: here we choose the *closest RT* to the reference (min |ΔRT|),
%     whereas get_histone1 chooses the *largest inten_sum* within the ±0.5 window.
%
% Inputs
%   MS1_index  : [nScans × ≥3], with MS1_index(:,2) = scan RT (min), MS1_index(:,3) = start row in MS1_peaks.
%   MS1_peaks  : [nRows × 2], concatenated MS1 centroids as [m/z, intensity] across scans.
%   ptol       : ppm tolerance. Historical convention: if ptol==100, it is coerced to 10.
%   unitdiff   : isotopic spacing per 13C step in Da (≈1.003355 Da), used as unitdiff/charge.
%   His        : struct with at least:
%                  • pep_mz  [nPeptides × nCharge]  — target m/z per charge state
%                  • pep_ch  [nPeptides × nCharge]  — charge (z) per state
%                  • rt_ref  [nPeptides × 1]        — reference RT per peptide (minutes)
%   hno        : scalar index of the peptide within His.
%
% Outputs
%   cur_rts            : [1 × nCharge] apex RT per charge (or reference RT if no candidate).
%   cur_intens         : [1 × nCharge] integrated area per charge (0 if no acceptable candidate).
%   cur_mono_isointens : column vector – monoisotopic XIC used for jno==1 (first charge state only).
%
% Assumptions / conventions
%   - Column 2 of the multi-trace XIC matrix corresponds to the monoisotopic trace.
%   - Scan window is symmetrically clamped around His.rt_ref(hno) by ±2 min for *all* charges.
%   - Overlap-aware integration is delegated to get_area (which itself may call judgeOverlap1).
%
% External dependencies
%   GetProfiles, GetTopBottom11, get_area
%

[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([1,ncharge]);
cur_intens = zeros([1,ncharge]);

% Define a tight ±2-min scan window around the peptide's reference RT.
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
    % -----------------------------
    % Build XICs for current charge
    % -----------------------------
    c_mz = His.pep_mz(hno,jno);
    c_ch = His.pep_ch(hno,jno);
    % Four isotope seeds: {M−1, M, M+1, M+2}, spaced by unitdiff/charge.
    c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
    % For coarse tolerances (>100 ppm) and higher charges (z>=3), enable a 13C-related mode in GetProfiles.
    if ptol>100 && c_ch>=3
        nC13 = 1;
    else
        nC13 = 0;
    end;
    [c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
    j = 2; % Column 2 holds the monoisotopic XIC by convention.
    c_mono_isointens = c_ref_isointens(:,j);
    if 1==jno
        cur_mono_isointens = c_mono_isointens; % Expose mono XIC for the first charge state only.
    end;

    % ----------------------------------------------
    % Candidate picking on the monoisotopic XIC
    % ----------------------------------------------
    [nt,nb,top1_idx,inten_sum] = GetTopBottom11(c_mono_isointens);%#ok
    % Filter candidates to those within ±0.5 min of the reference RT:
    flag = zeros([1,length(nt)]);
    for i=1:length(nt)
        if abs( c_isorts(nt(i))-His.rt_ref(hno) )<=0.5
            flag(i) = 1;
        end;
    end;
    x = find(flag==1);
    if 0==isempty(x)
        % Selection rule (distinctive of get_histone10):
        %   prefer the candidate with *minimum absolute deviation* from the reference RT.
        [tmp,id] = min(abs(c_isorts(nt(x))-His.rt_ref(hno)));%#ok
        %[tmp,id] = max(inten_sum(x));%#ok  % (alternative used in get_histone1)
        top1_idx = x(id);
        cur_pos = nt(top1_idx);
        cur_rts(jno) = c_isorts(cur_pos);

        % If there are at least two qualifying candidates, tighten the nb window to bracket them.
        if length(x)>=2
            nb = [nb(x(1)) nb(x(end)+1)];
        end;

        % Overlap-aware area integration at the apex (delegated to get_area).
        cur_intens(jno) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
    else
        % No acceptable candidate near the reference RT → fallback RT and zero area.
        cur_rts(jno) = His.rt_ref(hno);
        cur_intens(jno) = 0;
    end;
end;
