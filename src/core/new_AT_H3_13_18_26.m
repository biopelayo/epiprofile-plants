function H3_13_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ============================= OVERVIEW (EN) =================================
% This entry point performs targeted quantification for a fixed panel of four
% H3 peptides (see init_histone). It:
%   1) Skips work if an output .mat already exists (idempotent behavior).
%   2) Builds the peptide panel specification (seq/mods/charges/mz/RT anchors).
%   3) Computes representative retention times (RTs) and intensities from MS1,
%      optionally using MS2 to refine anchors (depending on 'special').
%   4) Saves results to disk and produces diagnostic plots.
%   5) Optionally extracts PSM evidence in data-assisted mode (nDAmode==1).
%
% Inputs:
%   MS1_index   : [nMS1 x 2] table; column 2 typically stores the RT (in minutes).
%   MS1_peaks   : structure/container with MS1 peaks (format assumed downstream).
%   MS2_index   : MS2 index table (used by drawing/PSM routines and DA refinements).
%   MS2_peaks   : structure/container with MS2 peaks.
%   ptol        : mass tolerance (ppm or Th, depending on downstream functions).
%   cur_outpath : output folder for <out_filename>.mat and figures.
%   special     : control struct. Common fields: ndebug, nDAmode, nhmass, raw_path.
%
% Side effects:
%   - Writes <cur_outpath>/H3_13_18_26.mat via output_histone().
%   - Generates figures via draw_layout().
%   - May write PSM evidence via GetPSM2() when special.nDAmode==1.

% ---------------------------------- check -----------------------------------
out_filename = 'H3_13_18_26';                                  % fixed stem used for all outputs
% fprintf(1,'%s..',out_filename);                              % optional console log (kept silent)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);       % expected output file path
if 0~=exist(out_file0,'file')                                   % if .mat already exists on disk...
    return;                                                     % ...quit early to avoid recomputation
end;

% ----------------------------------- init -----------------------------------
His = init_histone();                                           % define peptide panel:
                                                                %   - mod_short/mod_type encodings
                                                                %   - candidate charge states (pep_ch)
                                                                %   - theoretical m/z grid (pep_mz)
                                                                %   - RT anchors (rt_ref) and display flags

% -------------------------------- calculate ---------------------------------
unitdiff = 1.0032;                                              % ~C13-C12 isotopic mass difference
[pep_rts,pep_intens,mono_isointens] = ...                       % core quantification (RT + intensity)
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);                % also returns monoisotopic traces

% --------------------------------- output -----------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);% serialize results (His + matrices)

% ----------------------------------- draw -----------------------------------
num_MS1 = size(MS1_index,1);                                    % number of MS1 scans (rows)
isorts   = MS1_index(1:num_MS1,2);                              % RT vector used for x-axis/sorting
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...% multi-panel diagnostics/QC
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ---------------------------------- PSMs ------------------------------------
if 1==special.nDAmode                                           % when DA mode==1, extract PSM-level support
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
             isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

% ============================= LOCAL: init_histone ===========================
function His = init_histone()
%%
% Purpose:
%   Construct the histone-peptide panel: short labels, modification map,
%   charge states, theoretical m/z table, reference RTs, and display flags.

His.pep_seq = 'unmod';                                          % descriptive tag for the panel
His.mod_short = {'KQLATKAAR';                                   % four target peptides (short labels)
    'KELATKAAR';
    'TLLATKAAR';
    'KQLAPKAAR'};
His.mod_type = {'0,pr;1,pr;6,pr;';                              % modification encoding per peptide:
    '0,pr;1,pr;6,pr;';                                          %   - index,value pairs: e.g., '0,pr' → N-term
    '0,pr;6,pr;';                                               %   - 'pr' typically denotes propionylation
    '0,pr;1,pr;6,pr;'};                                         %     (see README note)

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);            % charges considered (z=1,2,3), same for all
%{
His.pep_mz = [900.5302  450.7687  300.8482                      % (example) hard-coded m/z table (commented)
    842.4955  421.7514  281.5034
    757.4315  379.2194  253.1487];
%}
His.pep_mz = calculate_pepmz(His);                               % derive m/z per charge from seq+mods
His.rt_ref = [39.63                                              % reference RTs (minutes) used as anchors
    17.36
    22
    22];
His.display = ones(length(His.mod_type),1);                      % 1 = include each peptide in outputs/plots

% Ensure the "main" charge appears in the first column of pep_mz/pep_ch
main_ch = His.pep_ch(1,2);                                      % here we pick z=2 as main (display preference)
if main_ch~=His.pep_ch(1,1)                                     % if main charge isn't the leftmost column...
    [npep,ncharge] = size(His.pep_mz);                          % ...reorder columns for all peptides
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];        % new order: main first, then the others
    x = zeros([1,ncharge]);                                     % mapping old→new column indices
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));            % locate the source column for each target slot
    end;
    tune = 1:npep;                                              % apply the same permutation to all peptide rows
    His.pep_mz(tune,:) = His.pep_mz(tune,x);                    % reorder m/z columns
    His.pep_ch(tune,:) = His.pep_ch(tune,x);                    % reorder charge-ID columns
end;

% ============================= LOCAL: calculate_layout =======================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% Purpose:
%   For each peptide and charge state:
%     - Determine a representative retention time (RT).
%     - Integrate its intensity (MS1-driven).
%     - Build a monoisotopic trace across the RT grid.
%   It also recalibrates RT anchors unless debug mode is requested.

[npep,ncharge] = size(His.pep_mz);                               % dimensions derived from the m/z table
num_MS1 = size(MS1_index,1);                                     % number of MS1 scans (rows)
pep_rts = zeros([npep,ncharge]);                                 % RT matrix (peptide × charge)
pep_intens = zeros([npep,ncharge]);                              % intensity matrix (peptide × charge)
mono_isointens = zeros([num_MS1,npep]);                          % monoisotopic trace: [scan × peptide]

% ------------------------ RT-anchor (rt_ref) calibration --------------------
His.rt_unmod_orig = His.rt_ref(1);                               % keep the original RT anchor (first peptide)
if 1~=special.ndebug                                             % allow recalibration unless in debug mode
    if 2~=special.nDAmode                                        % standard confirmation (nDAmode != 2)
        for hno=1:4                                              % iterate over the 4 peptides in this panel
            [His.rt_ref(hno),special.ndebug] = ...               % check/adjust anchor based on raw path
                check_ref(special.raw_path, ...
                          [His.mod_short{hno},His.mod_type{hno}], ...
                          His.rt_ref(hno), ...
                          special.ndebug);
        end;
    else                                                         % DA mode (nDAmode == 2): flexible refinement
        nhmass = special.nhmass;                                 % helper mass for get_rts2 (if used downstream)
        for hno=1:4
            rt_unmod_orig = His.rt_ref(hno);                     % store current candidate anchor
            His.rt_ref(hno) = check_ref(special.raw_path, ...
                                         [His.mod_short{hno},His.mod_type{hno}], ...
                                         His.rt_ref(hno), ...
                                         special.ndebug);        % try to confirm/update anchor
            if rt_unmod_orig==His.rt_ref(hno)                    % if anchor did not change (no clear signal)...
                t1 = 0;                                          % ...search the whole RT range
                t2 = MS1_index(num_MS1,2);                       % end at the last observed RT
            else                                                 % if anchor updated...
                delta = 5;                                       % ...focus search within ±5 minutes
                t1 = His.rt_ref(hno)-delta;
                t2 = His.rt_ref(hno)+delta;
            end;
            [rts1,top1_rt1] = ...                                % refine via RT profiling in the window
                get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                         ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok  % returned rts1 unused; top1_rt1 is key
            if 0==isempty(top1_rt1)                              % if a best-RT candidate exists...
                His.rt_ref(hno) = top1_rt1;                      % ...adopt it as the new anchor
            end;
        end;
        special.ndebug = 1;                                      % lock anchors after DA-mode refinement
    end;
end;

% Legacy notes (kept from original source):
% 64-69KLPFQR
% 129-134RIRGER
% 130-135IRGERA

% --------------------------- Per-peptide extraction -------------------------
for hno=1:4                                                      % loop over each of the 4 peptides
    [cur_rts,cur_intens,cur_mono_isointens] = ...               % get RTs, intensities, monoisotopic trace
        get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);
    if cur_rts(1)>0                                             % keep only if a valid RT was found
        pep_rts(hno,1:ncharge) = cur_rts;                       % store RTs for all charges
        pep_intens(hno,1:ncharge) = cur_intens;                 % store intensities for all charges
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;     % store monoisotopic trace (scan-wise)
    end;
end;
