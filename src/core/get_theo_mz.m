function [theo_mz,theo_tp] = get_theo_mz(c_sequence,c_modification,chg,Mods,ActiveType)
%% get_theo_mz — Theoretical fragment ion m/z calculator (b/y for CID, c/z for ETD)
%
% PURPOSE (plain English, non-expert friendly)
%   Given a peptide sequence and its modifications, this function computes the
%   *theoretical fragment ion masses* for the main fragmentation series used in
%   proteomics:
%     - CID/HCD-like:   b-ions (N-terminus) and y-ions (C-terminus)
%     - ETD-like:       c-ions (N-terminus) and z-ions (C-terminus)
%
%   It then converts neutral fragment masses to m/z values for charge states
%   z = 1..chg (inclusive), assuming protonation, and returns both:
%     - theo_mz : all computed fragment m/z values (concatenated over charges)
%     - theo_tp : an integer code that *labels each m/z* with its ion series,
%                 ordinal (b1, b2, … / y1, y2, …), and charge.
%
% INPUTS
%   c_sequence      char(1×L)   Peptide sequence in uppercase 1-letter code (A..Z).
%   c_modification  (type used by your codebase; typically a per-residue spec)
%                   Description of variable/fixed mods for this peptide. This is
%                   interpreted by 'get_mod_mass' together with 'Mods'.
%   chg             scalar int  Maximum charge for fragments to compute (z = 1..chg).
%   Mods            struct      Modification catalog used by get_mod_mass().
%   ActiveType      char        'CID'  -> compute b & y series
%                                'ETD' -> compute c & z series
%
% OUTPUTS
%   theo_mz  1×K double   Theoretical fragment m/z values across all series and charges.
%   theo_tp  1×K double   Encoded type for each theo_mz.
%
%   Type encoding (TP) convention:
%     For 'CID'  (b,y):
%       base code for b_i  : 1000 + 10*i
%       base code for y_i  : 2000 + 10*i
%     For 'ETD' (c,z):
%       base code for c_i  : 1000 + 10*i
%       base code for z_i  : 2000 + 10*i
%     Then, for charge z (1..chg), we add the charge as +z, so:
%       b3 at z=2   -> TP = 1000 + 10*3 + 2 = 1032
%       y7 at z=1   -> TP = 2000 + 10*7 + 1 = 2071
%     You can decode series = floor(TP/1000) (1 = N-series, 2 = C-series),
%     ordinal = mod(floor(TP/10),100), and charge = mod(TP,10).
%
% ASSUMPTIONS
%   - Amino-acid masses returned by Getaamass() are *residue* masses (i.e., the
%     mass the residue contributes when inside a peptide chain; typically AA mono
%     mass minus H2O for polymerization). This matches the b/y, c/z formulas used.
%   - get_mod_mass(c_sequence, c_modification, Mods) returns a per-residue delta-mass
%     vector aligned to c_sequence, which is added to the residue masses.
%   - Proton mass used is 1.007276 Da (common convention for m/z conversion).
%
% LIMITATIONS
%   - No immonium ions, a-ions, internal fragments, or neutral losses are computed
%     here—only the backbone series (b/y for CID; c/z for ETD).
%   - All m/z are monoisotopic and computed for z = 1..chg with +z protons added.

%% ---- 1) Define elemental constants used by the series formulas ----
% element masses (monoisotopic): [C, H, N, O, S]
element = [12, 1.0078246, 14.0030732, 15.9949141, 31.972070];

% Shorthand building blocks frequently seen in fragment formulas
mH   = element(2);
mCO  = element(1) + element(4);         %#ok<NASGU>  % C + O (not used below, kept for clarity)
mHO  = element(2) + element(4);                    % H + O
mH2O = element(2)*2 + element(4);                 % 2H + O  (water)
mNH3 = element(2)*3 + element(3);                 % 3H + N  (ammonia)
mNH2 = element(2)*2 + element(3);                 % 2H + N

%% ---- 2) Base peptide masses (residue-wise) + modifications ----
% Get per-letter residue masses (assumed to be "residue masses" inside a peptide)
aamass = Getaamass();              % expects a 26×1 (or ≥26) table for 'A'..'Z'
idx = c_sequence - 'A' + 1;        % map 'A'..'Z' → 1..26
residuemass = aamass(idx,1)';      % 1×L vector

% Modification deltas for each position (same length as c_sequence)
deltam = get_mod_mass(c_sequence, c_modification, Mods);

% Residue masses with modifications applied
residuemass = residuemass + deltam;

%% ---- 3) Fragment series neutral masses (before adding charge) ----
% Proton mass used to convert neutral mass → m/z: add z*pmass and divide by z
pmass = 1.007276;

if 1 == strcmp(ActiveType,'CID')
    % CID/HCD backbone series:
    %   b_i : cumulative sum from N-terminus, residues 1..i
    %   y_i : cumulative sum from C-terminus, residues L-i+1..L, PLUS H2O
    %         (y-ions carry the C-terminus, which includes -COOH → net +H2O vs residues)
    Mb = cumsum(residuemass(1:end-1));                                 % b1..b(L-1)
    My = cumsum(rot90(rot90(residuemass(2:end)))) + mH2O;              % y1..y(L-1)

    % Concatenate and build type codes (series+ordinal; charge is added later)
    M  = [Mb, My];
    TP = [ones(1,length(Mb))*1000 + (1:length(Mb))*10, ...
          2*ones(1,length(My))*1000 + (1:length(My))*10];

    % Remove duplicates if any (e.g., very short sequences / degeneracies)
    [M,ix] = unique(M);
    TP = TP(ix);

    % For each charge 1..chg, convert to m/z and encode charge in theo_tp
    theo_mz = repmat(0,[1,length(M)*chg]);
    theo_tp = repmat(0,[1,length(M)*chg]);
    for ino = 1:chg
        % m/z = (neutral_mass + z*proton_mass) / z
        theo_mz((ino-1)*length(M)+1 : ino*length(M)) = (M + ino*pmass) / ino;
        % add the charge as +ino on top of the base code
        theo_tp((ino-1)*length(M)+1 : ino*length(M)) = TP + ino;
    end
    % Final deduplication by m/z (keeps first occurrence and aligned TP)
    [theo_mz,ix] = unique(theo_mz);
    theo_tp = theo_tp(ix);

elseif 1 == strcmp(ActiveType,'ETD')
    % ETD/ECD backbone series:
    %   c_i : cumulative sum from N-terminus, residues 1..i, PLUS NH3
    %   z_i : cumulative sum from C-terminus, residues L-i+1..L, PLUS HO - NH2 + H
    %         (standard c/z definitions; the adjustment accounts for radical sites)
    Mc = cumsum(residuemass(1:end-1)) + mNH3;                          % c1..c(L-1)
    Mz = cumsum(rot90(rot90(residuemass(2:end)))) + mHO - mNH2 + mH;   % z1..z(L-1)

    M  = [Mc, Mz];
    TP = [ones(1,length(Mc))*1000 + (1:length(Mc))*10, ...
          2*ones(1,length(Mz))*1000 + (1:length(Mz))*10];
    [M,ix] = unique(M);
    TP = TP(ix);

    theo_mz = repmat(0,[1,length(M)*chg]);
    theo_tp = repmat(0,[1,length(M)*chg]);
    for ino = 1:chg
        theo_mz((ino-1)*length(M)+1 : ino*length(M)) = (M + ino*pmass) / ino;
        theo_tp((ino-1)*length(M)+1 : ino*length(M)) = TP + ino;
    end
    [theo_mz,ix] = unique(theo_mz);
    theo_tp = theo_tp(ix);

else
    % Unknown fragmentation mode: return empties to signal caller
    theo_mz = [];
    theo_tp =  [];
end
end
