function deltam = get_mod_mass(c_sequence,c_modification,Mods)
%GET_MOD_MASS
% ------------------------------------------------------------------------------
% Purpose (plain English, non-proteomics wording):
%   Convert a human-readable string of "modifications at positions" into a
%   numeric vector of per-position mass shifts (one value per residue
%   in the peptide sequence).
%
% What this function returns:
%   deltam  → 1 × peplen vector (double), where deltam(i) is the total mass
%             delta (in Daltons) to apply at residue i due to the listed
%             modifications.
%
% Inputs
%   c_sequence      : char array, peptide sequence (length = peplen).
%   c_modification  : char array describing mods in the format:
%                       "pos,ModName;pos,ModName;..."
%                     Example:
%                       "3,ac;5,me1;"  → acetylation at pos 3 and mono-methyl
%                                         at pos 5 (names must exist in Mods).
%                     Special positions (exactly as in the original code):
%                       - pos = 0         → applies at the N-terminus, mapped to
%                                           residue position 1.
%                       - pos = peplen+1  → applies at the C-terminus, mapped to
%                                           residue position peplen.
%                     If the string has NO comma ',', the function returns a
%                     zero vector (no mods applied). This is the original
%                     behavior kept intact.
%   Mods            : struct with the lookup tables of modifications:
%                       .mass(1, :)  → numeric row vector of masses
%                       .name{1, :}  → cell array of modification names
%                     The i-th name corresponds to the i-th mass.
%
% Output
%   deltam          : 1 × peplen double, per-residue mass offsets.
%
% Algorithm (step-by-step):
%   1) Initialize deltam = zeros(1, peplen).
%   2) If c_modification contains no comma ',', bail out early (no mods).
%   3) Parse pairs "pos,ModName" separated by semicolons ';'.
%   4) For each pair:
%        - Parse the numeric 'pos' between the previous ';' (or string start)
%          and the comma ','.
%        - Normalize boundary 'pos' (0 → 1; peplen+1 → peplen).
%        - Extract the ModName between ',' and the next ';'.
%        - Find ModName in Mods.name; if found, add Mods.mass at that
%          position into deltam(pos). If not found, reset deltam to zeros,
%          print a message, and return (original behavior).
%
% Important notes / constraints (unchanged from original):
%   - The parser is intentionally simple and assumes the exact syntax
%     "pos,Name;" blocks. Trailing semicolon is expected.
%   - Multiple modifications at the SAME position sum their masses.
%   - If ANY mod name is unknown, the function reverts to all zeros and exits.
%   - No try/catch or extra validation is added to keep the original behavior.
%
% Complexity
%   O(L + M), where L = length(c_modification) and M = number of mods.
%
% Examples (for humans, not executed here):
%   % Given c_sequence = 'PEPTIDE' (peplen = 7)
%   % c_modification = '0,ac;3,me1;8,ph;'
%   %   → N-term ac  (0 → position 1)
%   %   → pos 3 me1
%   %   → C-term ph  (peplen+1=8 → position 7)
%   % deltam will have non-zero entries at 1, 3, and 7.
% ------------------------------------------------------------------------------

    % Peptide length and output initialization
    peplen = length(c_sequence);
    deltam = repmat(0,[1,peplen]);

    % Quick exit: if there is no comma, original logic assumes "no mods"
    pos1 = strfind(c_modification,',');
    if 1==isempty(pos1);
        return;
    end;

    % Lookup tables
    masses = Mods.mass(1,:);
    names  = Mods.name;

    % Find semicolon separators; pos2 includes a virtual "0" as start anchor
    pos2 = [0 strfind(c_modification,';')];

    % Each block is "pos,Name;" → loop over commas to delimit pairs
    for jno = 1 : length(pos1)
        % Extract numeric position between previous ';' (or start) and the comma
        % Example: "...;3,me1;..." → substring "3"
        pos = str2num( c_modification(pos2(jno)+1 : pos1(jno)-1) ); %#ok<ST2NM>
        % Boundary mapping (keep the original semantics):
        %   0           → N-term mapped to residue 1
        %   peplen + 1  → C-term mapped to residue peplen
        if 0==pos
            pos = 1;
        elseif peplen+1==pos
            pos = peplen;
        end;

        % Extract modification name between comma and the next semicolon
        % Example: "...;3,me1;..." → substring "me1"
        cmod = c_modification(pos1(jno)+1 : pos2(jno+1)-1);

        % Lookup mod name in Mods.name
        [tf,loc] = ismember(cmod, names);

        if 1==tf
            % Add this mod mass to the selected residue position
            deltam(pos) = deltam(pos) + masses(loc);
        else
            % Unknown modification name → revert to all-zeros and exit
            deltam = repmat(0,[1,peplen]);
            fprintf(1,'no %s in Mods\n',cmod);
            return;
        end;
    end;
end
