function aamass = Getaamass()
%% Getaamass — Residue mass look-up table for peptide calculations (mono & average)
%
% PURPOSE (plain English, non-expert friendly)
%   Returns a 26×2 numeric table with *residue* masses (not free amino acid masses),
%   aligned to uppercase ASCII letters 'A'..'Z'. This is used to convert a peptide
%   sequence string into per-residue masses by indexing with:
%       idx    = c_sequence - 'A' + 1;      % map 'A'..'Z' -> 1..26
%       masses = aamass(idx, 1);            % column 1 = monoisotopic residue mass
%   Column 1 (mono) is intended for high-accuracy MS calculations (fragment ions),
%   Column 2 (avge) offers average masses if ever needed (rare in MS/MS scoring).
%
% IMPORTANT DEFINITIONS
%   *Residue mass* means “mass contribution inside a peptide chain”. In standard
%   proteomics, that is the amino-acid monoisotopic mass MINUS H2O (18.01056… Da)
%   to account for peptide bond formation. This matches b/y (CID) and c/z (ETD)
%   fragment ion formulas used elsewhere in the codebase.
%
% TABLE LAYOUT (by letter):
%   Row i corresponds to letter char('A'+i-1). Non-standard letters are included
%   for convenience/legacy; some map to surrogate values or 0.0 (unsupported).
%   The exact meaning of ambiguous letters (B, J, O, U, X, Z) is project-specific
%   and inherited from upstream code. See notes below.
%
%   Letter :  mono (Da)   avge (Da)    Notes
%   ---------------------------------------------------------------
%   A :  71.03711   71.0788     Alanine (Ala), residue mass
%   B : 114.53494  114.5962     Ambiguous/legacy (project-specific placeholder)
%   C : 103.00919  103.1388     Cysteine (Cys), *unmodified* residue
%   D : 115.02694  115.0886     Aspartic acid (Asp)
%   E : 129.04259  129.1155     Glutamic acid (Glu)
%   F : 147.06841  147.1766     Phenylalanine (Phe)
%   G :  57.02146   57.0519     Glycine (Gly)
%   H : 137.05891  137.1411     Histidine (His)
%   I : 113.08406  113.1594     Isoleucine (Ile)
%   J : 114.04293  114.1038     Legacy mapping (matches Asparagine mono mass)
%   K : 128.09496  128.1741     Lysine (Lys)
%   L : 113.08406  113.1594     Leucine (Leu) — same as I (isobaric I/L)
%   M : 131.04049  131.1926     Methionine (Met)
%   N : 114.04293  114.1038     Asparagine (Asn)
%   O : 114.07931  114.1472     Non-canonical (e.g., Ornithine; project-specific)
%   P :  97.05276   97.1167     Proline (Pro)
%   Q : 128.05858  128.1307     Glutamine (Gln)
%   R : 156.10111  156.1875     Arginine (Arg)
%   S :  87.03203   87.0782     Serine (Ser)
%   T : 101.04768  101.1051     Threonine (Thr)
%   U :   0.00000    0.0000     Reserved/unsupported (Selenocysteine not encoded here)
%   V :  99.06841   99.1326     Valine (Val)
%   W : 186.07931  186.2132     Tryptophan (Trp)
%   X : 113.08406  113.1594     Unknown → defaulted to I/L residue mass (legacy)
%   Y : 163.06333  163.1760     Tyrosine (Tyr)
%   Z : 128.55059  128.6231     Ambiguous/legacy (often Glx placeholder)
%
% PRACTICAL NOTES
%   - Ambiguous letters (B, J, O, X, Z) are *not* universal standards. They are
%     kept for backward compatibility with upstream data. If your project needs
%     strict Unimod/PSI conventions, redefine these rows accordingly.
%   - U (Selenocysteine) is set to 0.0 here (unsupported). If you need it,
%     replace row 21 with the Selenocysteine *residue* mass.
%   - Cysteine is *unmodified* here. Fixed/variable modifications (e.g., carbamidomethyl)
%     must be added via the modification layer (e.g., get_mod_mass).
%
% USAGE EXAMPLE
%   seq  = 'PEPTIDE';                     % uppercase peptide sequence
%   idx  = seq - 'A' + 1;                 % 1..26 indices
%   mono = Getaamass();                   % fetch table
%   mvec = mono(idx,1)';                  % monoisotopic residue masses for the seq
%   pep_mass_residues = sum(mvec);        % sum of residues (neutral backbone mass)
%   % For precursors/fragments you’ll add the appropriate end-group terms elsewhere.
%
% RETURNS
%   aamass : 26×2 double
%       col 1 = monoisotopic residue masses
%       col 2 = average residue masses

% aamass table
% two columns: (1)mono, (2)avge
aamass = [ ...
     71.03711   71.0788  % A
    114.53494  114.5962  % B (ambiguous/legacy)
    103.00919  103.1388  % C
    115.02694  115.0886  % D
    129.04259  129.1155  % E
    147.06841  147.1766  % F
     57.02146   57.0519  % G
    137.05891  137.1411  % H
    113.08406  113.1594  % I
    114.04293  114.1038  % J (legacy → Asn mass)
    128.09496  128.1741  % K
    113.08406  113.1594  % L (same as I)
    131.04049  131.1926  % M
    114.04293  114.1038  % N
    114.07931  114.1472  % O (non-canonical; e.g., Ornithine)
     97.05276   97.1167  % P
    128.05858  128.1307  % Q
    156.10111  156.1875  % R
     87.03203   87.0782  % S
    101.04768  101.1051  % T
      0.00000    0.0000  % U (selenocysteine not supported here)
     99.06841   99.1326  % V
    186.07931  186.2132  % W
    113.08406  113.1594  % X (unknown → I/L)
    163.06333  163.1760  % Y
    128.55059  128.6231  % Z (ambiguous/legacy; often Glx)
]; %#ok
end
