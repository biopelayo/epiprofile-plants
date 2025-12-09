function [K1,K2] = get_key_ions2(His, pno, qno, Mods, ActiveType)
%GET_KEY_IONS2
% ------------------------------------------------------------------------------
% Purpose (non-proteomics wording):
%   Build two *aligned* lists of theoretical fragment ion m/z values (K1, K2)
%   for the SAME peptide sequence under TWO different modification states
%   (isoforms) indexed by 'pno' and 'qno', **using ONLY C-terminal series**
%   fragments. Each index in K1 corresponds to the SAME backbone cut as the
%   same index in K2. These pairs serve as "diagnostic" ion matches to
%   discriminate between the two isoforms in MS/MS data.
%
% IMPORTANT: Unlike get_key_ions.m, this function intentionally *ignores*
%            N-terminal fragments (the original N-series loop is omitted,
%            exactly as in your source). We preserve that behavior.
%
% How it works (mental model):
%   1) Decode modifications for isoforms p and q on the same sequence:
%        positions (modpos) and types (modtype) via get_mod_postype.
%      Encode each site as POS*10000 + TYPE to compare sets easily.
%   2) Identify positions where p and q differ (set difference over codes),
%      and derive a contiguous position window [t1 .. t2-1] around them.
%   3) Compute theoretical fragment ions for both isoforms with get_theo_mz,
%      then split them into N-series (codes <2000) and C-series (codes >2000).
%      Positions are recovered from "type codes" as:
%         N-series: pos = floor((tp - 1000)/10)
%         C-series: pos = floor((tp - 2000)/10)
%   4) **C-terminal only**: For each C-terminal backbone position mirrored
%      in the peptide (posc), look up the matching theoretical ion in both
%      isoforms and store the paired m/z into K1/K2 at the same index.
%   5) Trim out entries with no match (zeros) to keep aligned K1/K2 vectors.
%
% Inputs:
%   His         : struct with peptide/isoform info
%                   .pep_seq      → peptide sequence (char array)
%                   .mod_type{}   → cell array of mod strings per isoform
%   pno, qno    : 1-based indices of the two isoforms to compare
%   Mods        : modification definitions (library/labels used downstream)
%   ActiveType  : fragmentation type string ('CID','ETD',...) passed to get_theo_mz
%
% Outputs:
%   K1, K2      : row vectors (1×L) of aligned theoretical m/z values
%                 from C-terminal ion series only. K1(i) pairs with K2(i).
%
% Notes:
%   - No numerical behavior is altered vs your original function.
%   - Charge fed to get_theo_mz is 1 (as in original code path).
%   - If K1/K2 end up empty, the C-terminal positions in the selected window
%     may not yield shared ion types between the two isoforms.
% ------------------------------------------------------------------------------

    % ----- Inputs unpacking
    c_seq  = His.pep_seq;
    c_mod1 = His.mod_type{pno};
    c_mod2 = His.mod_type{qno};

    % --------------------------------------------------------------------------
    % 1) Map modifications on sequence for both isoforms (position & type)
    %    Encode POS+TYPE → POS*10000 + TYPE to perform set algebra
    % --------------------------------------------------------------------------
    [modpos1, modtype1] = get_mod_postype(c_seq, c_mod1, Mods);
    [modpos2, modtype2] = get_mod_postype(c_seq, c_mod2, Mods);

    newtype1  = modpos1 * 10000 + modtype1;    % isoform p codes
    newtype2  = modpos2 * 10000 + modtype2;    % isoform q codes
    comtype   = intersect(newtype1, newtype2); % shared sites (pos & type)
    uniontype = union(newtype1, newtype2);     % any site in either isoform
    lefttype  = setdiff(uniontype, comtype);   % differences only

    % Derive inclusive window bounds [t1, t2) in backbone positions
    t1 = floor(min(lefttype) / 10000);
    if 0 == t1
        t1 = 1;  % guard against 0 (positions are 1-based)
    end
    t2 = floor(max(lefttype) / 10000);

    % Preallocate output length: twice the span (N- and C-series lanes).
    % We will fill ONLY the C-terminal portion (by design).
    nlen = 2 * (t2 - t1);
    K1   = repmat(0, [1, nlen]);
    K2   = repmat(0, [1, nlen]);

    % Backbone position ranges:
    %   N-terminal (unused in this function): t1 .. t2-1 → posn
    %   C-terminal (used): mirrored indices at peptide C-end → posc
    peplen = length(c_seq);
    posn   = t1:(t2 - 1);                      %#ok<NASGU>  % kept for symmetry
    posc   = (peplen - t2 + 1) : (peplen - t1);            % used below

    % --------------------------------------------------------------------------
    % 2) Theoretical ions for isoforms p and q (activation-dependent)
    % --------------------------------------------------------------------------
    [theo_mz1, theo_tp1] = get_theo_mz(c_seq, c_mod1, 1, Mods, ActiveType);
    [theo_mz2, theo_tp2] = get_theo_mz(c_seq, c_mod2, 1, Mods, ActiveType);

    % Split by series and decode positions
    %   For isoform p:
    n1   = find(theo_tp1 < 2000);
    tpn1 = floor((theo_tp1(n1) - 1000) / 10);   %#ok<NASGU>  % N-series (unused)
    c1   = find(theo_tp1 > 2000);
    tpc1 = floor((theo_tp1(c1) - 2000) / 10);   % C-series (used)

    %   For isoform q:
    n2   = find(theo_tp2 < 2000);
    tpn2 = floor((theo_tp2(n2) - 1000) / 10);   %#ok<NASGU>  % N-series (unused)
    c2   = find(theo_tp2 > 2000);
    tpc2 = floor((theo_tp2(c2) - 2000) / 10);   % C-series (used)

    % --------------------------------------------------------------------------
    % 3) C-terminal only: collect aligned m/z pairs across mirrored positions
    %    (We intentionally keep the original N-series loop removed/commented.)
    % --------------------------------------------------------------------------

    % ORIGINAL N-series loop intentionally omitted to match source behavior:
    %{
    for ino = 1:length(posn)
        [tf1, loc1] = ismember(posn(ino), tpn1);
        [tf2, loc2] = ismember(posn(ino), tpn2);
        if tf1 == 1 && tf2 == 1
            K1(ino) = theo_mz1(n1(loc1));
            K2(ino) = theo_mz2(n2(loc2));
        end
    end
    %}

    % C-terminal (mirrored) positions
    for ino = 1:length(posc)
        [tf1, loc1] = ismember(posc(ino), tpc1);
        [tf2, loc2] = ismember(posc(ino), tpc2);
        if tf1 == 1 && tf2 == 1
            % shift by length(posn) to place C-series after N-series lane
            K1(ino + length(posn)) = theo_mz1(c1(loc1));  % isoform p
            K2(ino + length(posn)) = theo_mz2(c2(loc2));  % isoform q
        end
    end

    % --------------------------------------------------------------------------
    % 4) Keep only valid aligned entries (drop zeros)
    % --------------------------------------------------------------------------
    x  = find(K1 > 0);
    K1 = K1(x);
    K2 = K2(x);
end
