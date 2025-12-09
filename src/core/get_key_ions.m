function [K1,K2] = get_key_ions(His, pno, qno, Mods, ActiveType)
%GET_KEY_IONS
% -------------------------------------------------------------------------
% Purpose (non-proteomics wording):
%   Build two *aligned* lists of theoretical fragment ion m/z values (K1, K2)
%   for the SAME peptide sequence but TWO different modification states
%   (isoforms) indexed by 'pno' and 'qno'. Each index in K1 corresponds to
%   the SAME backbone cut site as the same index in K2. These pairs are
%   used downstream to compute "diagnostic ion" ratios between the two states.
%
% How it works (mental model):
%   1) Identify where and what modifications are present for each isoform:
%        [modpos, modtype] for isoform p and q (via get_mod_postype).
%      We encode each site as POS*10000 + TYPE to compactly compare.
%   2) Compute the set of *different* modification "codes" between p and q.
%      From the min/max "codes" we derive a contiguous window of positions
%      [t1 .. t2-1] (N-terminal side) and the mirrored C-terminal positions.
%   3) Compute theoretical fragment ions for p and q (via get_theo_mz),
%      conditioned on the activation type (CID/ETD). The function returns:
%         - theo_mz: m/z values for all candidate ions of that isoform
%         - theo_tp: encoded ion "type" and position, using a simple scheme:
%            * N-terminal series  → codes in [1000..2000), position encoded
%            * C-terminal series  → codes > 2000,   position encoded
%      We map theo_tp back to integer backbone positions (tpn*, tpc*).
%   4) For every backbone position within the chosen window, we look up
%      the matching theoretical ion in isoform p and isoform q and store
%      the two m/z into K1 and K2 at the same index.
%   5) Finally, we remove any positions where K1==0 (no ion found), keeping
%      K1 and K2 aligned.
%
% Inputs:
%   His         : structure that includes
%                   .pep_seq      → peptide sequence (char array)
%                   .mod_type{}   → cell array of mod strings, one per isoform
%   pno, qno    : indices (1-based) of the two isoforms to compare (same peptide)
%   Mods        : structure describing modification definitions/labels
%   ActiveType  : string 'CID' or 'ETD' (or similar) passed to get_theo_mz
%
% Outputs:
%   K1, K2      : row vectors (1 × L) of m/z values; same length L, aligned
%                 by backbone position. K1(i) pairs with K2(i).
%
% Behaviour:
%   This function preserves the original EpiProfile logic/thresholds and
%   array shapes. Only comments and readability were improved.
% -------------------------------------------------------------------------

    % ----- Basic fields from inputs
    c_seq  = His.pep_seq;             % peptide sequence (same for p and q)
    c_mod1 = His.mod_type{pno};       % modification string for isoform p
    c_mod2 = His.mod_type{qno};       % modification string for isoform q

    % ---------------------------------------------------------------------
    % 1) Map modifications to (position, type) on the sequence
    %    newtype encodes both in a single integer: POS*10000 + TYPE
    %    This allows easy set algebra (intersect/union/diff).
    % ---------------------------------------------------------------------
    [modpos1, modtype1] = get_mod_postype(c_seq, c_mod1, Mods);
    [modpos2, modtype2] = get_mod_postype(c_seq, c_mod2, Mods);

    newtype1  = modpos1 * 10000 + modtype1;     % isoform p codes
    newtype2  = modpos2 * 10000 + modtype2;     % isoform q codes
    comtype   = intersect(newtype1, newtype2);  % shared (same pos & type)
    uniontype = union(newtype1, newtype2);      % any present in either
    lefttype  = setdiff(uniontype, comtype);    % the "differences" only

    % If there are differences, derive the inclusive position window [t1, t2)
    % by projecting code→position (floor(code/10000)).
    t1 = floor(min(lefttype) / 10000);
    if 0 == t1
        t1 = 1;  % guard for N-terminal edge (positions are 1-based)
    end
    t2 = floor(max(lefttype) / 10000);

    % nlen is twice the span (we will collect both N- and C-terminal series)
    nlen = 2 * (t2 - t1);
    K1   = repmat(0, [1, nlen]);  % isoform p theoretical m/z per site
    K2   = repmat(0, [1, nlen]);  % isoform q theoretical m/z per site

    % Backbone position ranges (1-based):
    %  - N-terminal cuts: positions t1 .. t2-1
    %  - C-terminal cuts: mirrored at the C-terminus
    peplen = length(c_seq);
    posn   = t1:(t2 - 1);
    posc   = (peplen - t2 + 1) : (peplen - t1);

    % ---------------------------------------------------------------------
    % 2) Compute theoretical ions for both isoforms (activation-dependent)
    %    get_theo_mz returns:
    %      theo_mzX: m/z values of candidate fragments
    %      theo_tpX: encoded ion "type" + position:
    %         [1000..2000): N-series  → position = floor((tp-1000)/10)
    %         > 2000      : C-series  → position = floor((tp-2000)/10)
    %    Charge state is fixed to 1 here (per original code path).
    % ---------------------------------------------------------------------
    [theo_mz1, theo_tp1] = get_theo_mz(c_seq, c_mod1, 1, Mods, ActiveType);
    [theo_mz2, theo_tp2] = get_theo_mz(c_seq, c_mod2, 1, Mods, ActiveType);

    % Split indexes by series (N- vs C-terminal) for isoform p
    n1   = find(theo_tp1 < 2000);
    tpn1 = floor((theo_tp1(n1) - 1000) / 10);    % N-series positions for p
    c1   = find(theo_tp1 > 2000);
    tpc1 = floor((theo_tp1(c1) - 2000) / 10);    % C-series positions for p

    % Same for isoform q
    n2   = find(theo_tp2 < 2000);
    tpn2 = floor((theo_tp2(n2) - 1000) / 10);    % N-series positions for q
    c2   = find(theo_tp2 > 2000);
    tpc2 = floor((theo_tp2(c2) - 2000) / 10);    % C-series positions for q

    % ---------------------------------------------------------------------
    % 3) For each backbone cut within the window, pick the matching
    %    theoretical ion in p and q and store paired m/z into K1/K2.
    %    First loop: N-terminal series (posn), second loop: C-terminal (posc).
    % ---------------------------------------------------------------------

    % N-terminal (positions t1 .. t2-1)
    for ino = 1:length(posn)
        % Is this position available in both isoforms' N-series?
        [tf1, loc1] = ismember(posn(ino), tpn1);
        [tf2, loc2] = ismember(posn(ino), tpn2);
        if tf1 == 1 && tf2 == 1
            % map local positions back to indices in theo_mz arrays
            K1(ino) = theo_mz1(n1(loc1));  % isoform p
            K2(ino) = theo_mz2(n2(loc2));  % isoform q
        end
    end

    % C-terminal (mirrored positions at the peptide C-terminus)
    for ino = 1:length(posc)
        [tf1, loc1] = ismember(posc(ino), tpc1);
        [tf2, loc2] = ismember(posc(ino), tpc2);
        if tf1 == 1 && tf2 == 1
            K1(ino + length(posn)) = theo_mz1(c1(loc1));  % isoform p
            K2(ino + length(posn)) = theo_mz2(c2(loc2));  % isoform q
        end
    end

    % ---------------------------------------------------------------------
    % 4) Drop any empty slots (K1==0) to keep only valid aligned pairs
    % ---------------------------------------------------------------------
    x  = find(K1 > 0);
    K1 = K1(x);
    K2 = K2(x);
end
