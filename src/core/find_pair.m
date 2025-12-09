function [t3,t4] = find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3,rts4,top1_rt4,inten_sum4,top1_inten_sum4)
%%
% FIND_PAIR  Resolve a matched RT pair (t3, t4) between two RT lists using asymmetric windows & intensity.
%
% Syntax:
%   [t3, t4] = find_pair(rts3, top1_rt3, inten_sum3, top1_inten_sum3, ...
%                        rts4, top1_rt4, inten_sum4, top1_inten_sum4)
%
% Purpose:
%   Given two candidate retention-time (RT) landscapes (channel/group "3" vs "4"),
%   choose a **paired** RT (t3 for group 3, t4 for group 4). The decision rule:
%     1) If both lists are non-empty:
%          - If the top RTs are already "close" (with an **asymmetric** window),
%            pair them directly:  t3 = top1_rt3; t4 = top1_rt4.
%          - Else, the group whose top peak is **more intense** anchors the pair:
%               • If group 3 is stronger: fix t3 = top1_rt3 and **search in rts4**
%                 within a window centered at top1_rt3 to pick the **most intense** candidate.
%                 If none found, fallback to t4 = top1_rt3 − 0.5 (synthetic offset).
%               • If group 4 is stronger: symmetric idea with roles swapped (search in rts3
%                 around top1_rt4; else fallback t3 = top1_rt4 + 0.5; and set t4 = top1_rt4).
%     2) If only one list is non-empty:
%          - If only rts3 exists:  t3 = top1_rt3;  t4 = top1_rt3 − 0.5 (synthetic mate).
%          - If only rts4 exists:  t3 = top1_rt4 + 0.5;  t4 = top1_rt4.
%     3) If both are empty: return zeros (t3 = 0; t4 = 0).
%
% Inputs (all required):
%   - rts3 : vector of candidate RTs for group 3.
%   - top1_rt3 : scalar RT of the top (best) candidate in group 3.
%   - inten_sum3 : vector of intensities aligned with rts3 (same length).
%   - top1_inten_sum3 : scalar intensity of the top candidate in group 3.
%   - rts4 : vector of candidate RTs for group 4.
%   - top1_rt4 : scalar RT of the top (best) candidate in group 4.
%   - inten_sum4 : vector of intensities aligned with rts4 (same length).
%   - top1_inten_sum4 : scalar intensity of the top candidate in group 4.
%
% Output:
%   - t3, t4 : selected paired RTs (double). May be a real candidate or a synthetic
%              offset (±0.5) when a mate cannot be found in the opposite list.
%
% Matching windows (asymmetric by design):
%   - When anchoring on top1_rt3 (searching in rts4):    [top1_rt3 − 2,  top1_rt3 + 0.5]
%   - When anchoring on top1_rt4 (searching in rts3):    [top1_rt4 − 0.5, top1_rt4 + 2]
%   Rationale [Inference]: this asymmetry likely models systematic elution shifts between
%   the two channels/groups (e.g., labeling or chromatography-dependent bias).
%
% Tie-breaking:
%   - Within a window, choose the **maximum intensity** (max(inten_sum?)) among candidates.
%   - If no candidate exists inside the window, use a **fixed offset** (±0.5) relative to the anchor.
%
% Dependencies:
%   - MATLAB built-ins: isempty, find, max.
%
% Caveats:
%   - Assumes inten_sum? and rts? vectors are aligned element-wise.
%   - No unit is enforced for RTs; minutes are typical in LC–MS (context-dependent).
%   - If multiple maxima exist, MATLAB max returns the first index.
%

if 0==isempty(rts3) && 0==isempty(rts4)
    % Case A: both lists have candidates → try direct proximity first, then intensity-anchored search.
    if top1_rt4>top1_rt3-2 && top1_rt4<top1_rt3+0.5
        % Direct proximity: top RTs are "close enough" within the asymmetric window around top1_rt3
        % → accept the pair as-is.
        t3 = top1_rt3;
        t4 = top1_rt4;
    elseif top1_inten_sum3>top1_inten_sum4
        % Group 3 dominates by top intensity → anchor on top1_rt3 and search in rts4.
        t3 = top1_rt3;
        id = find( rts4>top1_rt3-2 & rts4<top1_rt3+0.5 );
        if 0==isempty(id)
            [tmp,ix] = max(inten_sum4(id));%#ok
            t4 = rts4(id(ix));
        else
            % No mate found in the window → synthetic fallback offset (-0.5) from anchor.
            t4 = top1_rt3-0.5;
        end;
    else
        % Group 4 dominates by top intensity → anchor on top1_rt4 and search in rts3.
        id = find( rts3>top1_rt4-0.5 & rts3<top1_rt4+2 );
        if 0==isempty(id)
            [tmp,ix] = max(inten_sum3(id));%#ok
            t3 = rts3(id(ix));
        else
            % No mate found in the window → synthetic fallback offset (+0.5) from anchor.
            t3 = top1_rt4+0.5;
        end;
        t4 = top1_rt4;
    end;
elseif 0==isempty(rts3)
    % Case B: only rts3 has candidates → pair with a synthetic mate at -0.5.
    t3 = top1_rt3;
    t4 = top1_rt3-0.5;
elseif 0==isempty(rts4)
    % Case C: only rts4 has candidates → pair with a synthetic mate at +0.5.
    t3 = top1_rt4+0.5;
    t4 = top1_rt4;
else
    % Case D: both lists empty → default zeros.
    t3 = 0;
    t4 = 0;
end;
