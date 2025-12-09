function xt = find_triple(rts1,top1_rt1,rts2,rts3,inten_sum1,inten_sum2,inten_sum3)
%%
% FIND_TRIPLE  Build a 3×3 RT grid "xt" anchoring on list #1 and pairing lists #2/#3.
%
% Purpose:
%   Given three candidate RT lists (with their per-candidate intensities), assemble a small
%   3×3 map "xt" that places:
%     - Column 1: the primary anchor from rts1 (top1_rt1) and a near-time partner from rts2.
%     - Column 2: right-of-anchor picks from rts1/rts2/rts3 with specific windowing rules.
%     - Column 3: later candidates from rts2/rts3 seeded by xt(1,3) (derived in col 2).
%   This function **does not** infer anything beyond the original heuristics; it preserves
%   the exact thresholds/windows and decision flow of the source you provided.
%
% Inputs:
%   rts1, rts2, rts3     : vectors of candidate retention times for 3 lists.
%   top1_rt1             : scalar anchor RT for list #1 (must be non-empty to proceed).
%   inten_sum1..3        : intensity vectors aligned with rts1..3 (same lengths respectively).
%
% Output:
%   xt (3×3 double)      : RT grid. Zeros mean "not assigned".
%                          Row 1 ← rts1, Row 2 ← rts2, Row 3 ← rts3.
%                          Columns are left/middle/right temporal positions built by the rules below.
%
% Invariants kept from your original code:
%   - Same windows and thresholds:
%       * rts2 in (top1_rt1-1.2 , top1_rt1+0.7) for col #1 pairing
%       * rts1 ≥ top1_rt1+0.7 to seek the second anchor on the right
%       * dual-peak logic with >= 1/10 intensity ratio and < 1.3 RT apart in rts2 (col #2)
%       * xt(1,3) = max(two rts2 peaks) or first rts2 + 0.02 seed
%       * rts3 windows mirrored to rts2 around top1_rt1 and xt(1,2)
%       * optional midpoint snap for xt(1,2) toward mean(xt(2,2), xt(3,2))
%       * early return if column #2 is inconsistent
%       * col #3 windows and ±0.5 preference for rts3 around xt(2,3)
%   - No refactors affecting behavior; only comments added.
%

% Abort early if the primary anchor is missing.
if 1==isempty(top1_rt1)
    xt = zeros([3,3]);
    return;
end

% --- Original (commented) alternative anchor tweak preserved verbatim for provenance ---
%{
[tmp_sum,ix] = sort(inten_sum1,'descend');
tmp_rts = rts1(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<12
    top1_rt1 = min([tmp_rts(2),tmp_rts(1)]);
end
%}

% Delegate to the worker that follows exactly the original selection logic.
xt = find_triple_once(rts1,top1_rt1,rts2,rts3,inten_sum1,inten_sum2,inten_sum3);

end % <-- end main function


function xt = find_triple_once(rts1,top1_rt1,rts2,rts3,inten_sum1,inten_sum2,inten_sum3)
%%
% FIND_TRIPLE_ONCE  Internal worker implementing the original heuristic grid assembly.
%
% Column 1:
%   - xt(1,1) := top1_rt1 (seed from list #1)
%   - xt(2,1) := best rts2 within (top1_rt1-1.2 , top1_rt1+0.7)
%
% Column 2:
%   - xt(1,2) := rts1 candidate to the RIGHT (≥ top1_rt1+0.7)
%                If the second-strongest is reasonably intense (>= 1/10 of strongest),
%                prefer the earlier of the two provided rts2 has a neighbor within 0.7.
%   - xt(2,2) := rts2 best in [top1_rt1+0.7 , xt(1,2)+1.2] if xt(1,2)>0,
%                otherwise in [top1_rt1+0.7 , top1_rt1+1.9].
%                If two strong/close rts2 peaks exist (>=1/10 and <1.3 apart),
%                xt(2,2) = earlier; xt(1,3) = later. Else xt(2,2)=best; xt(1,3)=best+0.02.
%   - xt(3,2) := rts3 best in the symmetric window (mirroring the rts2 window).
%   - Optional: snap xt(1,2) toward the midpoint of xt(2,2) and xt(3,2) when a close rts1 exists.
%   - Early exit if col #2 is incoherent: (xt(1,2)=0 AND |xt(2,2)-xt(3,2)|>0.7) OR missing xt(2,2) or xt(3,2).
%
% Column 3:
%   - xt(2,3) := best rts2 with rts2 ≥ xt(1,3)+1
%   - xt(3,3) := rts3 in [xt(1,3)+1, xt(2,3)+1.2]; prefer within ±0.5 of xt(2,3);
%                otherwise coerce to xt(2,3).

xt = zeros([3,3]);

% --------------------------
% 1) First column (col #1)
% --------------------------
xt(1,1) = top1_rt1;

id = find( rts2>top1_rt1-1.2 & rts2<top1_rt1+0.7 );
if 0==isempty(id)
    [~,ix] = max(inten_sum2(id)); % pick the strongest rts2 near the anchor
    xt(2,1) = rts2(id(ix));
    % Note: the original left 'top1_rt1' unchanged; we keep that invariant.
end

% --------------------------
% 2) Second column (col #2)
% --------------------------
% Row 1 (rts1): look to the right of the anchor
id = find(rts1>=top1_rt1+0.7);
if 0==isempty(id)
    % Consider top-2 by intensity to optionally prefer the earlier if rts2 supports it.
    [tmp_sum,ix] = sort(inten_sum1(id),'descend');
    tmp_rts = rts1(id(ix));
    if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/10
        new_tmp_rt = min([tmp_rts(2),tmp_rts(1)]);
        if 0==isempty(find( abs(rts2-new_tmp_rt)<0.7 ))
            xt(1,2) = new_tmp_rt;
        else
            xt(1,2) = tmp_rts(1);
        end
    else
        xt(1,2) = tmp_rts(1);
    end
end

% Row 2 (rts2): window depends on whether xt(1,2) was set
if xt(1,2)>0
    id = find(rts2>=top1_rt1+0.7 & rts2<=xt(1,2)+1.2);
else
    id = find(rts2>=top1_rt1+0.7 & rts2<=top1_rt1+1.9);
end
if 0==isempty(id)
    [tmp_sum,ix] = sort(inten_sum2(id),'descend');
    tmp_rts = rts2(id(ix));
    if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/10 && abs(tmp_rts(2)-tmp_rts(1))<1.3
        % Two strong & close peaks: earlier goes to xt(2,2), later seeds xt(1,3)
        xt(2,2) = min([tmp_rts(2),tmp_rts(1)]);
        xt(1,3) = max([tmp_rts(2),tmp_rts(1)]);
    else
        % Single dominant peak: use it, and seed xt(1,3) just to the right (+0.02)
        xt(2,2) = tmp_rts(1);
        xt(1,3) = tmp_rts(1)+0.02;
    end
end

% Row 3 (rts3): symmetric windowing relative to top1_rt1 and xt(1,2)
if xt(1,2)>0
    id = find(rts3>=top1_rt1+0.7 & rts3<=xt(1,2)+1.2);
else
    id = find(rts3>=top1_rt1+0.7 & rts3<=top1_rt1+1.9);
end
if 0==isempty(id)
    [~,ix] = max(inten_sum3(id));
    xt(3,2) = rts3(id(ix));
end

% Optional midpoint snap for xt(1,2): try to place it near mean(xt(2,2), xt(3,2))
id = find(abs(rts1-(xt(2,2)+xt(3,2))/2)<0.7);
if 0==isempty(id)
    [~,ix] = min(abs(rts1(id)-(xt(2,2)+xt(3,2))/2));
    xt(1,2) = rts1(id(ix));
end

% Early exit: require a coherent 2nd column
if (0==xt(1,2) && abs(xt(2,2)-xt(3,2))>0.7) || 0==xt(2,2) || 0==xt(3,2)
    return;
end

% --------------------------
% 3) Third column (col #3)
% --------------------------
% Row 2 (rts2): later candidate
id = find(rts2>=xt(1,3)+1);
if 0==isempty(id)
    [~,ix] = max(inten_sum2(id));
    xt(2,3) = rts2(id(ix));
end

% Row 3 (rts3): prefer within ±0.5 of xt(2,3); else coerce to xt(2,3)
id = find(rts3>=xt(1,3)+1 & rts3<=xt(2,3)+1.2);
if 0==isempty(id)
    j = find(abs(rts3(id)-xt(2,3))<=0.5);
    if 1==isempty(j)
        xt(3,3) = xt(2,3);
    else
        [~,ix] = max(inten_sum3(id(j)));
        xt(3,3) = rts3(id(j(ix)));
    end
end

end % <-- end nested worker
