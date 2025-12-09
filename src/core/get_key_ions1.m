function [K1,posn,posc] = get_key_ions1(His,pno,Mods,ActiveType)
%GET_KEY_IONS1
% ------------------------------------------------------------------------------
% Purpose (non-proteomics wording):
%   Compute a *single* aligned list of theoretical fragment ion m/z values (K1)
%   for ONE peptide isoform (index 'pno'), covering both N-terminal and
%   C-terminal backbone cuts that lie around the region where this isoform
%   carries modifications.
%
%   In addition to K1, the function returns the *backbone positions* used for:
%     - N-terminal series (posn): forward cuts (1..L-1, subset)
%     - C-terminal series (posc): mirrored cuts from the peptide C-end
%
%   Each entry in K1 corresponds to a *specific backbone cut*:
%     K1(1:length(posn))       -> N-series ions at positions 'posn'
%     K1(length(posn)+1:end)   -> C-series ions at positions 'posc'
%
% What this is for:
%   You can later search experimental MS/MS spectra near these m/z to confirm
%   the presence of the isoform and/or to compare with other isoforms that have
%   different modification patterns.
%
% Key ideas (mental model):
%   1) Resolve the peptide sequence for isoform p (including special case
%      "unmod" → extract from His.mod_short{pno} after the final dot).
%   2) Map the *positions* of modifications on the sequence (no types needed).
%   3) Define a compact backbone window [t1..t2) that spans the modified area
%      (or the whole peptide if mods are absent/very close).
%   4) Generate theoretical fragments for this isoform (get_theo_mz), and decode
%      which entries are N-series (code<2000) or C-series (code>2000), and which
%      backbone position each one refers to.
%   5) For each backbone position in the chosen window, if a theoretical ion
%      exists, record its m/z in K1 (first N-lane, then C-lane).
%
% Inputs
%   His         : struct with isoform info
%                   .pep_seq        → peptide sequence (char)
%                   .mod_short{}    → human-readable labels per isoform
%                   .mod_type{}     → machine labels per isoform (for engines)
%   pno         : index (1-based) of the isoform to compute
%   Mods        : modification definitions library
%   ActiveType  : fragmentation type ('CID','ETD',...) forwarded to get_theo_mz
%
% Outputs
%   K1          : row vector of m/z values; first N-series then C-series
%   posn        : integer backbone positions used for N-series lane
%   posc        : integer backbone positions used for C-series lane (mirrored)
%
% Notes
%   - Absolutely no numerical behavior changed vs your original source.
%   - Charge fed to get_theo_mz is fixed to 1 (as in your code).
%   - get_mod_postype is called with a single output: only positions are needed.
%   - If some positions have no matching theoretical ion, zeros are left in K1.
% ------------------------------------------------------------------------------

    % ----- 1) Resolve the peptide sequence for isoform p -----------------------
    c_seq = His.pep_seq;

    % Special case: some pipelines encode the raw 'sequence' as 'unmod' and
    % keep the actual unmodified sequence within mod_short{pno} after a dot.
    if 1==strcmp(c_seq,'unmod')
        c_seq = His.mod_short{pno};
        p = strfind(c_seq,'.');
        if 0==isempty(p)
            c_seq = c_seq(p(end)+1:end);
        end;
    end;

    % Isoform p modification label (engine-oriented string)
    c_mod1 = His.mod_type{pno};

    % ----- 2) Map modification POSITIONS on the sequence -----------------------
    % Only positions are needed here (no mod types for this function).
    % get_mod_postype returns the 1-based positions where modifications apply.
    modpos1 = get_mod_postype(c_seq,c_mod1,Mods);

    peplen = length(c_seq);
    ix = find(modpos1>0);

    % Define the window [t1..t2) of backbone cuts to focus on:
    % - If no modifications or they are within <=3 residues apart:
    %     use the whole peptide span.
    % - Otherwise:
    %     use min..max modified positions to focus on the modified region.
    if 1==isempty(ix) || modpos1(ix(end))-modpos1(ix(1))<=3
        t1 = 1;
        t2 = peplen;
    else
        t1 = min( modpos1(ix) );
        t2 = max( modpos1(ix) );
    end;

    % Preallocate K1 for both N and C lanes over the chosen span.
    % The number of cuts in [t1..t2) is (t2-t1), and we reserve space for
    % two lanes (N and C), hence 2*(t2-t1).
    nlen = 2*(t2-t1);
    K1 = zeros([1,nlen]);

    % Backbone positions to probe:
    %   posn → N-terminal forward cuts within the window (t1..t2-1)
    %   posc → mirrored C-terminal cuts w.r.t. peptide length
    posn = t1:t2-1;
    posc = (peplen-t2+1):(peplen-t1);

    % ----- 3) Get theoretical fragments for this isoform -----------------------
    % theo_mz1 : m/z values
    % theo_tp1 : integer codes storing both the ion series and the backbone pos
    %            Convention used here:
    %               N-series → code < 2000,   pos = floor((code-1000)/10)
    %               C-series → code > 2000,   pos = floor((code-2000)/10)
    [theo_mz1,theo_tp1] = get_theo_mz(c_seq,c_mod1,1,Mods,ActiveType);

    % Split into N and C series and decode their backbone positions
    n1   = find(theo_tp1<2000);
    tpn1 = floor( (theo_tp1(n1)-1000)/10 );   % positions for N-series ions
    c1   = find(theo_tp1>2000);
    tpc1 = floor( (theo_tp1(c1)-2000)/10 );   % positions for C-series ions

    % ----- 4) Fill the N-series lane in K1 (if ion exists at that position) ----
    for ino=1:length(posn)
        [tf1,loc1] = ismember(posn(ino),tpn1);
        if 1==tf1
            K1(ino) = theo_mz1(n1(loc1));
        end;
    end;

    % ----- 5) Fill the C-series lane in K1 (shifted after the N-lane) ----------
    for ino=1:length(posc)
        [tf1,loc1] = ismember(posc(ino),tpc1);
        if 1==tf1
            K1(ino+length(posn)) = theo_mz1(c1(loc1));
        end;
    end;
end
