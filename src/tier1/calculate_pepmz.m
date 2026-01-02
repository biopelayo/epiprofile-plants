function pep_mz = calculate_pepmz(His)
%%
% CALCULATE_PEPMZ  Compute theoretical m/z values for a set of peptides across multiple charge states.
%
% Syntax:
%   pep_mz = calculate_pepmz(His)
%
% Description:
%   Computes the theoretical mass-to-charge (m/z) ratios for histone peptides—possibly bearing
%   post-translational modifications—over a matrix of charge states. Neutral peptide mass is
%   formed by summing per-residue monoisotopic masses (from Getaamass) plus one H2O (terminals).
%   Modification mass shifts are added residue-wise via get_mod_mass(c_seq, c_mod, Mods).
%   Final m/z uses:
%       m/z = (Mr + z * m_proton) / z
%   where Mr is the neutral peptide mass (with H2O and deltas), z is charge, and m_proton = 1.007276 Da.
%
% Input (struct 'His' expected fields):
%   - His.mod_type : cell array (N×1 or N×M) of modification specs; each row i is passed to get_mod_mass.
%   - His.pep_ch   : numeric matrix (N×K) with charge states; entry (i,j) is the charge used for row i, col j.
%   - His.pep_seq  : char. If NOT equal to 'unmod', this is the *common* base sequence for all rows.
%                    If equal to 'unmod', per-row sequences are taken from His.mod_short{i} (post-trimming, see below).
%   - His.mod_short: cell array (N×1) used only if His.pep_seq == 'unmod'. For each i, take His.mod_short{i}
%                    and, if it contains dots '.', keep ONLY the substring *after the last dot*.
%
% Output:
%   - pep_mz : double matrix (N×K). Entry (i,j) is the theoretical m/z for peptide i at charge column j.
%
% Dependencies (called inside):
%   - Mods   = GetMods();              % returns modification definitions
%   - aamass = Getaamass();            % returns monoisotopic AA masses indexed by 'A'..'Z' (≥26×1)
%   - deltam = get_mod_mass(c_seq, c_mod, Mods);  % returns per-residue delta vector (same length as c_seq)
%
% Constants:
%   - element = [12 1.0078246 14.0030732 15.9949141 31.972070]; % C, H, N, O, S monoisotopic (H=idx 2, O=idx 4)
%   - mH2O    = 2*H + O = element(2)*2 + element(4);           % one water per peptide for termini
%   - pmass   = 1.007276;                                      % proton mass (Da)
%
% Important notes & caveats:
%   (1) Sequence source:
%       - If His.pep_seq ~= 'unmod': the SAME sequence is used for every row.
%       - If His.pep_seq == 'unmod': sequence is taken row-wise from His.mod_short{i}, then trimmed after the last '.'.
%   (2) Dot-trimming:
%       - Keeps substring AFTER the last '.' (e.g., 'K.PEPTIDE.R' -> 'R'). Confirm this matches your data convention.
%   (3) Residue indexing:
%       - idx = c_seq - 'A' + 1 assumes uppercase letters A..Z only. Any other symbol will break indexing.
%   (4) Mod deltas:
%       - get_mod_mass must return one delta per residue (length == numel(c_seq)) to add element-wise.
%
% Example (pseudo):
%   His.pep_seq   = 'STELLARPEPTIDE';
%   His.mod_type  = {'unmodified'; 'K9ac'};
%   His.mod_short = {'x'; 'y'}; % ignored here because pep_seq ~= 'unmod'
%   His.pep_ch    = [2 3; 2 3];
%   mz = calculate_pepmz(His);
%

pep_mz = repmat(0,[size(His.mod_type,1),size(His.pep_ch,2)]);
% Preallocate output (N×K) where:
%   N = size(His.mod_type,1)
%   K = size(His.pep_ch,2)

Mods = GetMods();
% Retrieve modification catalog/definitions (passed into get_mod_mass).

aamass = Getaamass();
% Monoisotopic AA masses; addressed via indices derived from letters 'A'..'Z'.

element = [12 1.0078246 14.0030732 15.9949141 31.972070];% element mass
% Monoisotopic masses (C, H, N, O, S). Only H (idx 2) and O (idx 4) are used below.

mH2O = element(2)*2 + element(4);
% Mass of one H2O to account for peptide termini.

pmass = 1.007276;
% Proton mass for m/z calculation.

if 0==strcmp(His.pep_seq,'unmod')
    % Branch A: common sequence for all rows (pep_seq ~= 'unmod').

    c_seq = His.pep_seq;
    idx = c_seq-'A'+1;
    % Map uppercase letters to 1..26 indices.

    residuemass = aamass(idx,1)';
    % Row vector of base monoisotopic masses for residues in c_seq.

    for ino=1:size(His.mod_type,1)
        c_mod = His.mod_type{ino};
        deltam = get_mod_mass(c_seq,c_mod,Mods);
        % Per-residue mass shifts for the given modification spec.

        % peptide+modification
        residuemass_new = residuemass + deltam;
        % Add deltas residue-wise.

        Mr = sum(residuemass_new)+mH2O;
        % Neutral peptide mass including one H2O.

        for jno=1:size(His.pep_ch,2)
            c_ch = His.pep_ch(ino,jno);
            pep_mz(ino,jno) = (Mr+c_ch*pmass)/c_ch;
            % m/z for this peptide and this charge column.
        end;
    end;
else
    % Branch B: per-row sequence from His.mod_short{ino} (pep_seq == 'unmod').

    for ino=1:size(His.mod_type,1)
        c_seq = His.mod_short{ino};
        p = strfind(c_seq,'.');
        if 0==isempty(p)
            c_seq = c_seq(p(end)+1:end);
            % Keep only substring AFTER the last '.' (data-format specific).
        end;
        idx = c_seq-'A'+1;
        residuemass = aamass(idx,1)';

        c_mod = His.mod_type{ino};
        deltam = get_mod_mass(c_seq,c_mod,Mods);
        % Per-residue shifts for this specific sequence.

        % peptide+modification
        residuemass_new = residuemass + deltam;

        Mr = sum(residuemass_new)+mH2O;
        % Neutral mass for this row's sequence.

        for jno=1:size(His.pep_ch,2)
            c_ch = His.pep_ch(ino,jno);
            pep_mz(ino,jno) = (Mr+c_ch*pmass)/c_ch;
            % Final m/z at the given charge state.
        end;
    end;
end;
