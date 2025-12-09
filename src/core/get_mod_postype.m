function [modpos,modtype] = get_mod_postype(c_sequence,c_modification,Mods)
%GET_MOD_POSTYPE
% ------------------------------------------------------------------------------
% Purpose (plain English, non-proteomics wording):
%   Parse a compact text string that lists "which positions carry which
%   modification names" and convert it into two numeric arrays:
%     - modpos  : the positions along the peptide sequence
%     - modtype : the index of each modification name in Mods.name
%
% What this function returns:
%   modpos   → 1×N vector of integer positions (N = number of "pos,Name" pairs)
%   modtype  → 1×N vector of integer codes; modtype(i) is the position of
%              the i-th modification name inside the lookup table Mods.name.
%
% Inputs:
%   c_sequence      : char array with the peptide sequence. Let peplen = length(c_sequence).
%   c_modification  : char array with pairs in the exact format:
%                       "pos,Name;pos,Name;..."
%                     Examples:
%                       "3,ac;5,me1;"  → position 3 has "ac", position 5 has "me1".
%                       "0,ac;8,ph;"   → (see notes below about boundary positions)
%   Mods            : struct with lookup tables:
%                       .name  → 1×M cell array of valid modification names
%                       .mass  → 1×M numeric vector of the corresponding masses
%                     (This function ONLY uses Mods.name; Mods.mass is not
%                     accessed here but is part of the same lookup struct.)
%
% Output details:
%   - modpos(i) is the position parsed from each "pos,Name;" pair.
%   - modtype(i) is the integer index such that Mods.name{modtype(i)} equals
%     the parsed Name. If a name is not found, the function prints a message
%     and returns immediately (leaving whatever was filled so far).
%
% Boundary handling (exactly as in the original code, kept intact):
%   - If a position equals peplen + 1 (i.e., sequence length + 1),
%       it is mapped to peplen (C-terminus mapped to the last residue).
%   - Positions equal to 0 are NOT remapped here. They remain as 0 in modpos,
%       which upstream/downstream logic may interpret as N-terminus.
%       (Other functions may later normalize or treat 0 specially.)
%
% Early exit rule (unchanged from original):
%   - If c_modification contains no comma ',', the function returns modpos = []
%     and modtype = [] (empty vectors).
%
% Parsing strategy (step-by-step):
%   1) Locate all commas (',') → each comma marks one "pos,Name" pair.
%   2) Locate all semicolons (';') and also treat position 0 as the start anchor.
%   3) For each comma:
%        - Read the integer 'pos' between the previous ';' (or start) and the comma.
%        - If pos == peplen + 1, remap to pos = peplen (C-terminus rule).
%        - Read the 'Name' between the comma and the next ';'.
%        - Find 'Name' inside Mods.name via ismember; store the index in modtype.
%        - If not found, print "no <Name> in Mods" and return immediately.
%
% Important notes:
%   - This function does not sum or combine modifications at the same site; it
%     simply enumerates them. Combination/aggregation is handled elsewhere.
%   - The input string format must strictly use ";" as pair terminator and ","
%     as the separator between position and name. A trailing ';' is expected.
%   - No extra validation is added to keep original behavior.
%
% Example (for humans, not executed here):
%   c_sequence     = 'PEPTIDE';      % peplen = 7
%   c_modification = '0,ac;3,me1;8,ph;'
%   Mods.name      = {'ac','me1','ph', ...};
%   → modpos   might be [0 3 7]  (note the C-term 8 remapped to 7)
%   → modtype  might be [1 2 3]  (indexes into Mods.name)
% ------------------------------------------------------------------------------

    % Find commas: each comma indicates one "pos,Name" pair
    pos1 = strfind(c_modification,',');
    modpos  = zeros([1,length(pos1)]);
    modtype = zeros([1,length(pos1)]);

    % Early exit: if there's no comma, we consider there are no pairs to parse
    if 1==isempty(pos1);
        return;
    end;

    % Lookup table (names only; masses exist but are not used here)
    names = Mods.name;

    % Sequence length for C-terminal remap rule
    peplen = length(c_sequence);

    % Find semicolons; also add a virtual "0" to ease substring slicing
    pos2 = [0 strfind(c_modification,';')];

    % Loop over each comma → parse "pos,Name" between markers
    for jno = 1 : length(pos1)
        % Extract the numeric position between previous ';' (or start) and the comma
        % Example block: "...;3,me1;..." → substring "3"
        pos = str2num( c_modification(pos2(jno)+1:pos1(jno)-1) ); %#ok<ST2NM>

        % Map peplen+1 (C-terminus) to peplen (last residue)
        if peplen+1==pos
            pos = peplen;
        end;
        modpos(jno) = pos;

        % Extract the mod name between ',' and the next ';'
        % Example block: "...;3,me1;..." → substring "me1"
        cmod = c_modification(pos1(jno)+1:pos2(jno+1)-1);

        % Look up the mod name in Mods.name
        [tf,loc] = ismember(cmod,names);
        if 1==tf
            modtype(jno) = loc;
        else
            % Name not found → print message and return immediately
            fprintf(1,'no %s in Mods\n',cmod);
            return;
        end;
    end;
end
