function [nt,nb,top1_idx,inten_sum] = GetTopBottom(w)
%% GetTopBottom — Split a 1D intensity trace into peak windows and compute their apex & area
%
% PLAIN-ENGLISH PURPOSE (for non-experts)
%   Given a 1D signal 'w' (e.g., an MS1 monoisotopic extracted ion chromatogram),
%   this function:
%     1) suppresses near-zero noise,
%     2) finds continuous non-zero stretches (candidate peaks),
%     3) optionally splits stretches at "valleys" (deep local minima),
%     4) for each final stretch, computes:
%          - 'nt(i)': an apex index (here, intensity-weighted centroid)
%          - 'inten_sum(i)': the integrated intensity (area under the curve) in the stretch,
%     5) returns 'top1_idx' as the index of the most intense stretch (largest area).
%
% OUTPUTS
%   nt        : [nPeak×1] integer indices (per-peak), intensity-weighted centroids inside each window
%   nb        : [nPeak+1×1] boundary indices; each peak i spans nb(i) ... nb(i+1)
%   top1_idx  : scalar, the peak-number (1..nPeak) with the largest 'inten_sum'
%   inten_sum : [nPeak×1] area (sum of w) per peak window
%
% IMPORTANT
%   • This routine relies on 'JudgeLocalmaxmin' to propose interior split points by
%     locating local minima/maxima. If you swap that helper, tune thresholds below.
%   • All indices are 1-based (MATLAB). Boundaries are inclusive.
%
% HEURISTICS / THRESHOLDS
%   • Noise gate:     w(w < 0.015*max(w)) = 0    -> keep only ≥1.5% of global max
%   • Valley depth:   split at a local minimum if w(min) < 0.2 * (local segment max)
%   • Balance check:  accept the split if the areas on both sides are not too skewed:
%                     min(areaLeft, areaRight) > 0.05 * max(areaLeft, areaRight)
%   • Apex choice:    intensity-weighted centroid, not the hard argmax (more stable for flat tops)
%
% EDGE CASES
%   • If all values vanish after noise gating -> return empty outputs.
%   • If a window contains interior zeros near edges, endpoints are nudged inward to
%     exclude trailing/leading zeros.
%
% COMPLEXITY
%   • Linear in the number of points plus the cost of local-peak detection in the helper.

% --- 1) Noise gate: keep only ≥1.5% of the global maximum; treat the rest as zeros (separators).
% w = smooth(w,3);  % (optional smoothing disabled in this version)
w(find(w<0.015*max(w))) = 0; %#ok<FNDSB>

% --- 2) If everything is zero, nothing to do.
ix = find(w==0);
if length(w)==length(ix)
    nt = [];
    nb = [];
    top1_idx = [];
    inten_sum = [];
    return;
end

% --- 3) Build initial boundaries 'nb' from runs of zeros vs non-zeros.
% Locate gaps in zero positions (dx>1 indicates separation between zero blocks).
x1  = ix(1:end-1);
x2  = ix(2:end);
dx  = x2 - x1;
idx = find(dx>1);

if isempty(idx)
    % Entire trace is one non-zero block: peak spans from 1 to end.
    nb = [1; length(w)];
else
    % Start each non-zero block right after a zero run; always close with last index.
    nb = [ix(idx); length(w)];
end

% --- 4) First pass: split non-zero blocks at deep, well-balanced valleys.
for i = 1:length(nb)-1
    i1 = nb(i);
    i2 = nb(i+1);
    IX = i1:i2;

    c_max = max(w(IX));
    [t,b] = JudgeLocalmaxmin(w(IX)); %#ok<ASGLU>
    if length(b) > 2
        % Inspect interior candidates b(2...end-1) as potential valleys.
        for j = 2:length(b)-1
            nmid = IX(b(j));               % candidate split index in global coordinates
            if w(nmid) < 0.2*c_max         % deep enough valley (≤20% of local max)
                s1 = sum( w(i1:nmid) );    % area to the left
                s2 = sum( w(nmid:i2) );    % area to the right
                m1 = min(s1,s2);
                m2 = max(s1,s2);
                if m1 > 0.05*m2           % not too imbalanced (≥5% of the larger side)
                    nb = [nb; nmid]; %#ok<AGROW>  % accept split
                    break;
                end
            end
        end
    end
end
nb = sort(nb);

% --- 5) Second pass: re-check after adding boundaries (catch a missed interior valley once).
for i = 1:length(nb)-1
    i1 = nb(i);
    i2 = nb(i+1);
    IX = i1:i2;

    c_max = max(w(IX));
    [t,b] = JudgeLocalmaxmin(w(IX)); %#ok<ASGLU>
    if length(b) > 2
        for j = 2:length(b)-1
            nmid = IX(b(j));
            if w(nmid) < 0.2*c_max
                s1 = sum( w(i1:nmid) );
                s2 = sum( w(nmid:i2) );
                m1 = min(s1,s2);
                m2 = max(s1,s2);
                if m1 > 0.05*m2
                    nb = [nb; nmid]; %#ok<AGROW>
                    break;
                end
            end
        end
    end
end
nb = sort(nb);

% --- 6) For each final window, trim leading/trailing zeros, compute apex and area.
nt        = zeros([length(nb)-1,1]);
inten_sum = zeros([length(nb)-1,1]);

for i = 1:length(nb)-1
    i1 = nb(i);
    i2 = nb(i+1);

    % Nudge left boundary rightwards if it sits on zero; keep at least original i1.
    while i1<=length(w) && w(i1)==0
        i1 = i1 + 1;
    end
    if i1 > nb(i)
        i1 = i1 - 1;
    end

    % Nudge right boundary leftwards if it sits on zero; keep at least original i2.
    while i2>=1 && w(i2)==0
        i2 = i2 - 1;
    end
    if i2 < nb(i+1)
        i2 = i2 + 1;
    end

    IX = i1:i2;

    % Apex index as intensity-weighted centroid (robust for flat/shouldered peaks).
    % Equivalent to: argmax of convolution with index, normalized by sum(w(IX)).
    % Using 'floor' to return an integer index within IX.
    nt(i) = floor(IX * w(IX) / sum(w(IX)));

    % Update the tightened boundaries in 'nb' (so that nb(i):nb(i+1) matches the window).
    nb(i)   = i1;
    nb(i+1) = i2;

    % Integrated intensity (area) of the window.
    inten_sum(i) = sum( w(nb(i):nb(i+1)) );
end

% --- 7) Pick the most intense window by area.
[~,top1_idx] = max(inten_sum);
