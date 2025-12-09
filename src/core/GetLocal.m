function [localmax_rt,localmax_inten,IX] = GetLocal(cur_ms1pos,isorts,c_ref_isointens,nb)
%% GetLocal — Per-isotope local-maximum finder around a given MS1 apex index
%
% PLAIN-ENGLISH PURPOSE (for non-experts)
%   Given:
%     • cur_ms1pos        — index of a current MS1 apex (e.g., one candidate peak)
%     • isorts            — retention-time vector aligned to MS1 scans
%     • c_ref_isointens   — MS1 intensity matrix [nScans × nIsotopes] for the target
%                           m/z channels (e.g., M−1, M, M+1, M+2 ... columns)
%     • nb                — boundary indices of segmented regions (from a peak picker)
%
%   This function locates the *local region* (segment) that contains cur_ms1pos,
%   then, within that region, finds the **maximum intensity** and its **RT** for
%   every isotope channel. It returns:
%     localmax_rt(j)     — RT at which channel j reaches its local max (within region)
%     localmax_inten(j)  — intensity of that local max (channel j, within region)
%     IX                 — index range [i1:i2] of the region used for the search
%
% CONTEXT
%   In isotope-aware MS1 analysis, you often need to align apexes across isotopic
%   channels (M, M+1, M+2...) within a *single chromatographic feature*. This
%   helper extracts, for the region containing a known apex index, the per-channel
%   maxima so downstream code can compare M vs M+1 timing, co-elution, areas, etc.
%
% ASSUMPTIONS
%   - 'isorts' and rows of 'c_ref_isointens' share the same indexing and length.
%   - 'nb' contains region boundaries such that each feature is nb(k):nb(k+1).
%   - Column 1 of 'c_ref_isointens' is used as the guiding trace for denoising
%     and to robustly determine the active region (historical convention).
%
% BEHAVIOR NOTES
%   - Very small intensities in the guiding trace are zeroed (<1.5% of its max),
%     to avoid tiny tails splitting regions incorrectly. If that zeroing kills
%     the intensity at 'cur_ms1pos', the code *falls back* to the raw (unzeroed)
%     guiding trace to preserve the intended region.
%   - Inside each region, trailing/leading zeros are trimmed to tighten [i1,i2].
%   - Commented code includes an optional "split at deep valley" heuristic that
%     was disabled; we keep it commented for reproducibility.
%
% INPUTS
%   cur_ms1pos        scalar index (1-based) — candidate apex position in MS1 grid
%   isorts            [nScans × 1] double    — RTs corresponding to MS1 scans
%   c_ref_isointens   [nScans × nIsotopes]   — intensity per scan for each isotope channel
%   nb                [K × 1] int            — region boundary indices (length K >= 2)
%
% OUTPUTS
%   localmax_rt       [nIsotopes × 1] double — RT of local maximum per isotope channel
%   localmax_inten    [nIsotopes × 1] double — intensity of that maximum per channel
%   IX                [1 × L] int            — index vector of the region used (i1:i2);
%                                             empty if no region contains cur_ms1pos.

    % ---- 0) Initialize outputs ----
    nmz = size(c_ref_isointens,2);                 % number of isotope channels (columns)
    localmax_rt    = repmat(0,[nmz,1]);
    localmax_inten = repmat(0,[nmz,1]);
    IX = [];

    % ---- 1) Choose a guiding trace to define/clean the region ----
    % We use column 1 as "navigation" trace (historical convention). This may be M−1 or M,
    % depending on upstream GetProfiles; the goal is only to robustly localize the region.
    nw = c_ref_isointens(:,1);
    % Optional smoothing was considered historically; kept disabled for reproducibility:
    % nw = smooth(nw,3);

    % Zero out very small values (<1.5% of the max) to remove faint tails and noise
    nw(find(nw<0.015*max(nw))) = 0; %#ok<FNDSB>
    % If zeroing accidentally wipes out the apex position, revert to the raw signal
    if 0 == nw(cur_ms1pos)
        nw = c_ref_isointens(:,1);
    end

    % ---- 2) Find the region [i1,i2] in 'nb' that contains 'cur_ms1pos' ----
    bflag = 0;  % becomes 1 when we find the segment that contains cur_ms1pos
    for ino = 1:length(nb)-1
        i1 = nb(ino);
        i2 = nb(ino+1);

        % trim leading zeros within [i1,i2] in the guiding trace
        while i1<=length(nw) && 0==nw(i1)
            i1 = i1+1;
        end
        if i1>nb(ino)
            i1 = i1-1;
        end

        % trim trailing zeros within [i1,i2]
        while i2>=1 && 0==nw(i2)
            i2 = i2-1;
        end
        if i2<nb(ino+1)
            i2 = i2+1;
        end

        IX = i1:i2;

        % Does this trimmed region contain the current apex index?
        c_pos = find( IX==cur_ms1pos ); %#ok<FNDSB>
        if ~isempty(c_pos)
            bflag = 1;

            %{
            % --- Optional sub-splitting at deep valleys (disabled) ---
            % If the guiding trace has a deep valley within IX, we could split the
            % region at that valley to keep only the sub-interval containing cur_ms1pos.
            % Left here commented for traceability and potential future re-activation.
            c_max = max(w(IX));
            [t,b] = JudgeLocalmaxmin(w(IX));
            if length(b)>2
                for j=2:length(b)-1
                    mid = IX(b(j));
                    if nw(mid)<0.2*c_max
                        s1 = sum( w(i1:mid) );
                        s2 = sum( w(mid:i2) );
                        m1 = min(s1,s2);
                        m2 = max(s1,s2);
                        if m1>0.2*m2
                            if cur_ms1pos<mid
                                IX = i1:mid;
                            else
                                IX = mid:i2;
                            end
                            break;
                        end
                    end
                end
            end
            %}

            break;  % region found
        end
    end

    % If no region contains 'cur_ms1pos', leave outputs as zeros/empty IX and return
    if 0==bflag
        return;
    end

    % ---- 3) For each isotope channel, find its local max within IX ----
    for j = 1:nmz
        c_isointens = c_ref_isointens(:,j);
        [localmax_inten(j), x] = max(c_isointens(IX));   % local maximum and relative index in IX
        localmax_rt(j) = isorts( IX(x) );                % map to RT via isorts
    end
end
