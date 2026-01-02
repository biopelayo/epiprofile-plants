function [t3,t4] = find_pair_new(top1_rt3,rts4,top1_rt4,inten_sum4,rt3_up)
%%
% FIND_PAIR_NEW  Pair a reference RT from group 3 with a mate in group 4 using directional windows.
%
% Syntax:
%   [t3, t4] = find_pair_new(top1_rt3, rts4, top1_rt4, inten_sum4, rt3_up)
%
% Purpose:
%   Given the **top RT** of group 3 (top1_rt3) and a **list of RT candidates** for group 4 (rts4)
%   with their intensities (inten_sum4), decide a **paired** RT tuple (t3, t4). The pairing uses:
%     - An **asymmetric search window** around top1_rt3 whose orientation depends on rt3_up.
%     - A **max-intensity** selection within that window for t4.
%     - If no candidate exists in the window, a **synthetic fallback** offset (±0.5) from the anchor.
%   Edge cases where one of the top RTs is missing are handled with sensible defaults.
%
% Inputs:
%   - top1_rt3   : scalar double. RT of the top candidate in group 3 (anchor).
%   - rts4       : vector double. Candidate RTs for group 4.
%   - top1_rt4   : scalar double. RT of the top candidate in group 4 (used in edge cases).
%   - inten_sum4 : vector double. Intensities aligned index-wise with rts4.
%   - rt3_up     : scalar logical/flag (0 or 1). **Directionality of expected elution order**:
%                    1 → group 3 expected **later** than group 4 (t3 > t4).
%                        Search window: [top1_rt3 − 2, top1_rt3 + 0.5]; fallback t4 = top1_rt3 − 0.5.
%                    0 → group 3 expected **earlier** than group 4 (t3 < t4).
%                        Search window: [top1_rt3 − 0.5, top1_rt3 + 2]; fallback t4 = top1_rt3 + 0.5.
%
% Outputs:
%   - t3, t4 : scalar doubles. Paired RTs. If no mate is found, one member is a synthetic ±0.5 offset.
%
% Notes:
%   - This routine **always anchors on top1_rt3** when both top RTs are present; it does not compare
%     intensities across groups (unlike other variants). The selection by intensity happens **within rts4** only.
%   - Uses the idiom `0==isempty(x)` to mean “x is present/non-empty”.
%   - Assumes `numel(rts4)==numel(inten_sum4)` and aligned.
%

if 1==rt3_up
    % Expected order: t3 (group 3) elutes LATER than t4 (group 4).
    if 0==isempty(top1_rt3) && 0==isempty(top1_rt4)
        % Both top RTs available → anchor on top1_rt3 and try to find t4 near it (earlier-biased window).
        t3 = top1_rt3;
        id = find( rts4>top1_rt3-2 & rts4<top1_rt3+0.5 );
        if 0==isempty(id)
            [tmp,ix] = max(inten_sum4(id));%#ok
            t4 = rts4(id(ix));
        else
            % No candidate in window → synthetic mate earlier by 0.5.
            t4 = top1_rt3-0.5;
        end;
    elseif 0==isempty(top1_rt3)
        % Only group 3 present → pair with synthetic earlier mate.
        t3 = top1_rt3;
        t4 = top1_rt3-0.5;
    elseif 0==isempty(top1_rt4)
        % Only group 4 present → synthesize t3 later by 0.5 around top1_rt4.
        t3 = top1_rt4+0.5;
        t4 = top1_rt4;
    else
        % Neither present.
        t3 = 0;
        t4 = 0;
    end;
else
    % Expected order: t3 (group 3) elutes EARLIER than t4 (group 4).
    if 0==isempty(top1_rt3) && 0==isempty(top1_rt4)
        % Both top RTs available → anchor on top1_rt3 and try to find t4 near it (later-biased window).
        t3 = top1_rt3;
        id = find( rts4>top1_rt3-0.5 & rts4<top1_rt3+2 );
        if 0==isempty(id)
            [tmp,ix] = max(inten_sum4(id));%#ok
            t4 = rts4(id(ix));
        else
            % No candidate in window → synthetic mate later by 0.5.
            t4 = top1_rt3+0.5;
        end;
    elseif 0==isempty(top1_rt3)
        % Only group 3 present → pair with synthetic later mate.
        t3 = top1_rt3;
        t4 = top1_rt3+0.5;
    elseif 0==isempty(top1_rt4)
        % Only group 4 present → synthesize t3 earlier by 0.5 around top1_rt4.
        t3 = top1_rt4-0.5;
        t4 = top1_rt4;
    else
        % Neither present.
        t3 = 0;
        t4 = 0;
    end;
end;
