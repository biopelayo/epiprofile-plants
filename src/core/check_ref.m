function [rt_ref,ndebug] = check_ref(raw_path,new_seq,rt_ref,ndebug)
% CHECK_REF  Try to resolve a reference retention time (rt_ref) from a layouts MAT-file.
%
% Syntax:
%   [rt_ref, ndebug] = check_ref(raw_path, new_seq, rt_ref, ndebug)
%
% Purpose:
%   If debugging/normal mode allows (ndebug==0 or ndebug==-1) and the reference file
%   '<raw_path>/histone_layouts/0_ref_info.mat' exists, load it and try to find the
%   reference retention time for the peptide/modification combination encoded by 'new_seq'.
%   Matching proceeds in two stages:
%     (1) Fast lookup via a Gödel-like numeric key (field 'seq_godel').
%     (2) If multiple hits share the same key, disambiguate by exact string comparison
%         of [pep_seq, mod_type] against 'new_seq'.
%   On a successful match, overwrite 'rt_ref' with the stored value and set ndebug = -1
%   (a sentinel indicating “reference resolved/used”).
%
% Inputs:
%   - raw_path : char/string. Base directory containing subfolder 'histone_layouts'.
%   - new_seq  : char/string. Encoded sequence key for the query; must match the way
%                this code constructs 'cur_seq' from the reference:
%                  cur_seq = [AllUnHis.pep_seq{i,1}, AllUnHis.mod_type{i,1}]
%                i.e., peptide sequence concatenated with its modification spec.
%   - rt_ref   : double. Incoming/default reference RT; returned unchanged if no match is found.
%   - ndebug   : double. Control flag; lookup only runs if (ndebug==0 || ndebug==-1).
%
% Outputs:
%   - rt_ref   : double. If a match is found, updated to AllUnHis.rt_ref(i); otherwise unchanged.
%   - ndebug   : double. Set to -1 when a match is applied; unchanged if no match/lookup skipped.
%
% Reference file & expected fields:
%   '<raw_path>/histone_layouts/0_ref_info.mat' should define a struct 'AllUnHis' with fields:
%     - seq_godel : numeric vector. A Gödel-like key per entry (see below).
%     - rt_ref    : numeric vector. Stored reference retention times aligned to seq_godel entries.
%     - pep_seq   : cellstr (N×1). Peptide sequences.
%     - mod_type  : cellstr (N×1). Modification specifications.
%
% Gödel-like key (seq_godel):
%   new_godel = sum((new_seq - '0' + 49) .* log(2:1+length(new_seq)))
%   This maps each character of 'new_seq' to an integer weight ((char - '0' + 49)) and multiplies it
%   by log of successive integers (natural log), summing across positions. The result is a numeric
%   fingerprint intended to be comparable (==) against 'AllUnHis.seq_godel'.
%   Note: Comparing floating values with '==' can be sensitive to precision; here it assumes the
%   same formula/precision was used to build 'seq_godel' in the MAT-file.
%
% Early exit conditions:
%   - If ~(ndebug==0 || ndebug==-1), do nothing (lookup skipped).
%   - If the MAT-file does not exist, do nothing.
%   - If no 'seq_godel' match is found, do nothing.
%
% Side effects:
%   - Loads variables from MAT-file into function workspace (via 'load').
%
% Dependencies:
%   - Built-ins: fullfile, exist, load, sum, log, find, strcmp.
%   - Data file: '<raw_path>/histone_layouts/0_ref_info.mat' with 'AllUnHis'.
%

mat_file = fullfile(raw_path,'histone_layouts','0_ref_info.mat');
% Build path to the reference metadata file.

if (0==ndebug || -1==ndebug) && 0~=exist(mat_file,'file')
    % Proceed only if debugging mode allows and the mat-file exists.

    load(mat_file);
    % Load 'AllUnHis' (and possibly other variables) into the local workspace.

    new_godel = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
    % Compute the Gödel-like numeric key for 'new_seq'.
    %   - (new_seq - '0' + 49) maps each character to an integer weight.
    %   - log(2:1+length(new_seq)) uses natural logs of 2..(1+L), L = length(new_seq).
    %   - Sum of element-wise products yields the key.

    ix = find(AllUnHis.seq_godel==new_godel);
    % Find all entries whose stored key equals the computed one.

    if 1==isempty(ix)
        % No key match: leave rt_ref and ndebug unchanged; exit.
        return;
    elseif 1==length(ix)
        % Unique key match: adopt the stored rt_ref, mark ndebug = -1.
        rt_ref = AllUnHis.rt_ref(ix);
        ndebug = -1;
    else
        % Multiple entries share the same key: disambiguate by exact string identity.

        nflag = 0;
        for i=1:length(ix)
            cur_seq = [AllUnHis.pep_seq{ix(i),1},AllUnHis.mod_type{ix(i),1}];
            % Reconstruct the sequence+modification string for this candidate.

            if 1==strcmp(cur_seq,new_seq)
                nflag = 1;
                break;
            end;
        end;

        if 1==nflag
            % Exact string match found: adopt corresponding rt_ref; mark ndebug = -1.
            rt_ref = AllUnHis.rt_ref(ix(i));
            ndebug = -1;
        end;
        % If no exact match among candidates, leave outputs unchanged.
    end;
end;
