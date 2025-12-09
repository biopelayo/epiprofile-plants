function [top1_rt,top1_inten_sum] = get_rts0(MS1_index,MS1_peaks,ptol,His,hno)
%GET_RTS0  Estimate the dominant retention time (RT) of a target peptide from MS1
% ------------------------------------------------------------------------------
% Plain-English purpose (no proteomics jargon needed):
%   From a time-ordered list of spectra (MS1), this function builds the
%   intensity-vs-time profile for the peptide of interest and returns:
%     - top1_rt         : the time (in minutes) where the peptide signal peaks
%     - top1_inten_sum  : the summed intensity around that peak (a proxy for area)
%
% How it works (conceptually):
%   1) Compute the expected m/z of the peptide (and its natural isotope neighbors)
%      using its charge.
%   2) Extract a time profile (chromatogram) for the monoisotopic channel
%      by calling GetProfiles(...).
%   3) Find the main peak (time window where the curve is strongest)
%      using GetTopBottom(...).
%   4) Return the peak time and its summed intensity.
%
% Special case:
%   For a short list of specific peptide sequences (see 'str' below),
%   we pre-scan with a small +14.01565/charge m/z offset to guess an
%   upper scan boundary ('new_pos'). This shrinks the search range and
%   can stabilize peak picking in those cases.
%
% Inputs
%   MS1_index : [N×2] double
%               Column 1 = scan number (integer, increasing with time)
%               Column 2 = retention time in minutes
%   MS1_peaks : packed peak list used by GetProfiles (opaque here)
%   ptol      : mass tolerance (ppm). Convention in this code:
%                 - if ptol == 100 → use 10 ppm (historical, kept for compat.)
%                 - if ptol > 100 and charge >= 3 → enable C13-picking helper
%   His       : struct with fields used here:
%                 .pep_seq  {K×1}  cell, peptide sequences (strings)
%                 .pep_mz   [K×1]  double, precursor m/z for each peptide
%                 .pep_ch   [K×1]  double/int, charge state for each peptide
%   hno       : index of the peptide within His.* arrays (1-based)
%
% Outputs
%   top1_rt        : dominant peak retention time (min). If not found or < 4 min → 0
%   top1_inten_sum : summed intensity around the dominant local maximum (arbitrary units)
%
% Assumptions / relied-upon helpers
%   GetProfiles(MS1_index,MS1_peaks,isotope_mzs,charge,ptol,nC13,scan_range)
%     → [times, iso_intensities] where iso_intensities(:,2) is the monoisotopic channel
%   GetTopBottom(mono_trace)
%     → [nt, nb, top1_idx, inten_sum] where:
%         - nt are indices of detected local maxima
%         - nb are indices of local minima (valleys)
%         - top1_idx selects the "best" maximum among nt
%         - inten_sum(top1_idx) is an area-like intensity around that maximum
%
% Edge-case handling
%   - If nothing convincing is found or peak is very early (< 4 min),
%     outputs are set to zero.
%
% ------------------------------------------------------------------------------
% Author notes:
%   - This is a light-touch, documented version of the original logic.
%   - Numerical constants are preserved (e.g., 1.0032, 14.01565).
% ------------------------------------------------------------------------------

    % Natural isotopic spacing in m/z units (approx. 13C - 12C mass diff)
    unitdiff = 1.0032;

    % Total number of MS1 entries (time points/scans)
    num_MS1 = size(MS1_index,1);

    % Normalize ptol as in original code (keep historical behaviors)
    if ptol==100
        ptol = 10;
    end

    % Enable the "nC13" helper when tolerance is very wide and charge >= 3
    % (used inside GetProfiles to handle isotope picking)
    if ptol>100 && His.pep_ch(hno)>=3
        nC13 = 1;
    else
        nC13 = 0;
    end

    % --------------------------------------------------------------------------
    % 1) Optional pre-scan for a few specific sequences
    %    Rationale: for these sequences, look first at a slightly shifted m/z
    %    (+14.01565/charge) to tentatively find an upper scan boundary (new_pos).
    % --------------------------------------------------------------------------
    special_list = {'TKQTAR','KSTGGKAPR','KSAPATGGVKKPHR','KVLR'};
    if ismember(His.pep_seq{hno,1}, special_list)
        c_ch = His.pep_ch(hno);
        c_mz = His.pep_mz(hno) + 14.01565/c_ch;  % small offset pre-screen

        % Build the 4-isotope m/z list around the candidate center
        c_ref_isomzs = [c_mz-unitdiff/c_ch, c_mz, c_mz+unitdiff/c_ch, c_mz+2*unitdiff/c_ch];

        % Extract time trace for all scans
        [c_isorts, c_ref_isointens] = GetProfiles( ...
            MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, 1:num_MS1);

        % Use monoisotopic channel (column 2 by convention)
        c_mono_isointens = c_ref_isointens(:,2);

        % Find dominant maximum
        [nt, ~, top1_idx, inten_sum] = GetTopBottom(c_mono_isointens); %#ok<ASGLU>

        % Guard against empty peak lists (rare but safer)
        if isempty(nt) || isempty(top1_idx)
            top1_rt = 0;
        else
            top1_rt = c_isorts(nt(top1_idx));
        end

        % Discard ultra-early peaks (< 4 min) as unreliable
        if isempty(top1_rt) || top1_rt < 4
            top1_rt = 0;
        end

        % Translate the chosen RT into a scan boundary 'new_pos'
        if top1_rt == 0
            new_pos = num_MS1;
        else
            p = find(MS1_index(:,2) <= top1_rt);
            if isempty(p)
                new_pos = num_MS1;
            else
                new_pos = p(end);
            end
        end
    else
        % Default: search the full range
        new_pos = num_MS1;
    end

    % --------------------------------------------------------------------------
    % 2) Build the final MS1 profile at the *true* target m/z (no offset)
    % --------------------------------------------------------------------------
    c_mz = His.pep_mz(hno);
    c_ch = His.pep_ch(hno);

    % 4-isotope window centered at the target m/z
    c_ref_isomzs = [c_mz-unitdiff/c_ch, c_mz, c_mz+unitdiff/c_ch, c_mz+2*unitdiff/c_ch];

    % Extract up to new_pos (either full scan count or pre-screen boundary)
    [c_isorts, c_ref_isointens] = GetProfiles( ...
        MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, 1:new_pos);

    % Monoisotopic channel
    c_mono_isointens = c_ref_isointens(:,2);

    % Peak picking on the mono trace
    [nt, ~, top1_idx, inten_sum] = GetTopBottom(c_mono_isointens); %#ok<ASGLU>

    % Defensive checks
    if isempty(nt) || isempty(top1_idx)
        top1_rt = 0;
        top1_inten_sum = 0;
        return
    end

    % Pick dominant peak and its "area-like" intensity
    top1_rt = c_isorts(nt(top1_idx));
    top1_inten_sum = inten_sum(top1_idx);

    % Final early-elution guard
    if isempty(top1_rt) || top1_rt < 4
        top1_rt = 0;
        top1_inten_sum = 0;
    end
end
