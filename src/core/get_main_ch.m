function His = get_main_ch(His,main_ch_idx)
%GET_MAIN_CH
% ------------------------------------------------------------------------------
% Purpose (plain English, non-proteomics wording):
%   Reorder the peptide charge-state columns so that a chosen "main" charge
%   becomes the first column across all peptide entries in the structure 'His'.
%
%   In practice, 'His.pep_mz' and 'His.pep_ch' are matrices with:
%     - rows   → peptides / isoforms
%     - cols   → charge states available for that peptide
%   The function moves the column corresponding to 'main_ch_idx' to the front,
%   and keeps the remaining columns in their original relative order.
%
% Why this matters:
%   Downstream code often assumes that the FIRST column is the reference (main)
%   charge-state to compute MS1 profiles, intensities, etc. This utility makes
%   that assumption true by reordering columns once, centrally.
%
% Inputs
%   His           : struct with at least fields
%                     .pep_mz  (npep x ncharge)  → m/z matrix per peptide/charge
%                     .pep_ch  (npep x ncharge)  → integer charge per column
%   main_ch_idx   : 1-based column index selecting the desired "main" charge
%                   from the current layout (the same charge is assumed to be
%                   present for ALL rows in that column position)
%
% Output
%   His           : same struct, but columns are re-ordered so that the charge
%                   at 'main_ch_idx' becomes column 1. Both .pep_mz and .pep_ch
%                   are permuted identically to stay aligned.
%
% Contract / Assumptions (no runtime checks added to keep original logic):
%   - 1 <= main_ch_idx <= size(His.pep_ch,2)
%   - The first row of His.pep_ch lists the available charges; columns represent
%     the same charge across all rows (no per-row inconsistency).
%   - Each charge appears only once among columns (unique per column set).
%
% Behavior details (exactly as original code):
%   1) Read the "main charge" value from the FIRST row and the selected column:
%        main_ch = His.pep_ch(1,main_ch_idx)
%   2) If that value is already in column 1, do nothing.
%   3) Otherwise:
%        a) Build a new charge ordering: [main_ch, the rest...]
%        b) For each position in that ordering, find where that charge currently
%           sits in row 1, and record the source column index → vector 'x'.
%        c) Reindex BOTH matrices (.pep_mz and .pep_ch) by 'x' for ALL rows.
%
% Complexity
%   O(ncharge^2 + npep*ncharge), dominated by the final column permutation.
%
% Caveats
%   - No bounds checking is performed on 'main_ch_idx'.
%   - If duplicated charge labels existed across columns (unexpected), 'find'
%     would return the first match.
% ------------------------------------------------------------------------------

    % Read the charge value at the requested "main" column from the first row
    main_ch = His.pep_ch(1,main_ch_idx);

    % If the main charge is already the first column, nothing to do
    if main_ch~=His.pep_ch(1,1)

        % Dimensions: npep = number of peptide rows, ncharge = number of cols
        [npep,ncharge] = size(His.pep_mz);

        % Build the desired charge order:
        %   - first: the selected main charge
        %   - then : all other charges in their current order (set difference)
        new_ch = [main_ch, setdiff(His.pep_ch(1,:), main_ch)];

        % Map each desired charge in 'new_ch' to its CURRENT column index 'x'
        % by searching on the first row of 'pep_ch'
        x = zeros([1,ncharge]);
        for ino = 1:ncharge
            x(ino) = find(His.pep_ch(1,:) == new_ch(ino));
        end;

        % Reorder columns for ALL rows using the computed permutation 'x'
        tune = 1:npep;  % row index vector (all rows)
        His.pep_mz(tune,:) = His.pep_mz(tune, x);
        His.pep_ch(tune,:) = His.pep_ch(tune, x);
    end;
end
