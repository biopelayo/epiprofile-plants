function [nt,nb,top1_idx,inten_sum] = GetTopBottom11(w)
%% GetTopBottom11 — Minimal peak windowing without valley-based splitting
%
% PLAIN-ENGLISH PURPOSE (for non-experts)
%   Given a 1D intensity trace 'w' (e.g., an MS1 monoisotopic XIC), this
%   function segments the signal into contiguous non-zero windows and, for
%   each window, computes:
%     • 'nt(i)'       : an apex index (here, intensity-weighted centroid)
%     • 'inten_sum(i)': the window area (sum of intensities)
%   It also returns 'top1_idx' — the index of the most intense window.
%
% HOW THIS DIFFERS FROM GetTopBottom
%   • This “11” variant is a *simplified splitter*: it does NOT further split
%     a non-zero window at interior valleys (the code to do that is commented
%     out below). In other words, windows are defined purely by runs of
%     non-zeros created by a noise gate.
%   • Use this when you prefer conservative segmentation (fewer windows) and
%     maximum speed/reproducibility over aggressive deconvolution of shoulders.
%
% OUTPUTS
%   nt        : [nWin×1] apex indices (intensity-weighted centroids) per window
%   nb        : [nWin+1×1] boundaries; window i spans nb(i) ... nb(i+1) (inclusive)
%   top1_idx  : scalar, window index (1..nWin) with the largest 'inten_sum'
%   inten_sum : [nWin×1] area (sum of w) per window
%
% HEURISTICS
%   • Noise gate: w(w < 0.015*max(w)) = 0  (1.5% global max as a floor to
%     create zero gaps and suppress baseline noise).
%   • Apex: centroid rather than argmax -> more stable on flat/plateau peaks.
%
% EDGE CASE
%   • If everything becomes zero after the noise gate, all outputs are empty.

% -- 1) Noise gate: discard points below 1.5% of the global max (become zeros).
% w = smooth(w,3);  % optional smoothing is disabled here by design
w(find(w<0.015*max(w))) = 0; %#ok<FNDSB>

% -- 2) Trivial all-zero case
ix = find(w==0);
if length(w)==length(ix)
    nt = [];
    nb = [];
    top1_idx = [];
    inten_sum = [];
    return;
end

% -- 3) Build initial boundaries 'nb' from non-zero runs created by the gate
x1  = ix(1:end-1);
x2  = ix(2:end);
dx  = x2-x1;
idx = find(dx>1);

if isempty(idx)
    % one single non-zero run from 1 to end
    nb = [1; length(w)];
else
    % each gap in zeros indicates the start of a non-zero run
    nb = [ix(idx); length(w)];
end

% -------------------------------------------------------------------------
% (Disabled) Optional valley-based splitting inside non-zero windows:
% The following block—duplicated twice in GetTopBottom for robustness—scans
% each window for interior deep local minima to split overlapped peaks.
% In this “11” variant it is intentionally commented out to keep segmentation
% conservative and fast. If you need more aggressive peak deconvolution,
% consider using GetTopBottom instead, or re-enable this logic and tune
% thresholds (valley depth & side-area balance).
% -------------------------------------------------------------------------
%{
% First pass
% for i=1:length(nb)-1
%     i1 = nb(i); i2 = nb(i+1); IX = i1:i2;
%     c_max = max(w(IX));
%     [t,b] = JudgeLocalmaxmin(w(IX));
%     if length(b)>2
%         for j=2:length(b)-1
%             nmid = IX(b(j));
%             if w(nmid)<0.2*c_max
%                 s1 = sum( w(i1:nmid) ); s2 = sum( w(nmid:i2) );
%                 m1 = min(s1,s2); m2 = max(s1,s2);
%                 if m1>0.05*m2
%                     nb = [nb; nmid]; %#ok<AGROW>
%                     break;
%                 end
%             end
%         end
%     end
% end
% nb = sort(nb);
%
% % Second pass (after inserting boundaries)
% for i=1:length(nb)-1
%     i1 = nb(i); i2 = nb(i+1); IX = i1:i2;
%     c_max = max(w(IX));
%     [t,b] = JudgeLocalmaxmin(w(IX));
%     if length(b)>2
%         for j=2:length(b)-1
%             nmid = IX(b(j));
%             if w(nmid)<0.2*c_max
%                 s1 = sum( w(i1:nmid) ); s2 = sum( w(nmid:i2) );
%                 m1 = min(s1,s2); m2 = max(s1,s2);
%                 if m1>0.05*m2
%                     nb = [nb; nmid]; %#ok<AGROW>
%                     break;
%                 end
%             end
%         end
%     end
% end
% nb = sort(nb);
%}

% -- 4) For each window: tighten edges (remove leading/trailing zeros),
%       compute centroid apex and area.
nt        = zeros([length(nb)-1,1]);
inten_sum = zeros([length(nb)-1,1]);
for i=1:length(nb)-1
    i1 = nb(i);
    i2 = nb(i+1);

    % Trim zeros at the left edge (keep at least original i1)
    while i1<=length(w) && w(i1)==0
        i1 = i1+1;
    end
    if i1>nb(i)
        i1 = i1-1;
    end

    % Trim zeros at the right edge (keep at least original i2)
    while i2>=1 && w(i2)==0
        i2 = i2-1;
    end
    if i2<nb(i+1)
        i2 = i2+1;
    end

    IX = i1:i2;

    % Intensity-weighted centroid as apex (robust on flat/shouldered peaks)
    % Note: returns an integer index via floor()
    nt(i) = floor(IX * w(IX) / sum(w(IX)));

    % Update tightened window back into nb
    nb(i)   = i1;
    nb(i+1) = i2;

    % Area under the window
    inten_sum(i) = sum( w(nb(i):nb(i+1)) );
end

% -- 5) Pick the most intense (largest area) window
[~,top1_idx] = max(inten_sum);
