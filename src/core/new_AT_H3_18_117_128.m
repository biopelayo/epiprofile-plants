function H3_18_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ============================== OVERVIEW (EN) ================================
% Targeted quantification for a panel of three H3 peptides (see init_histone),
% across the charge states defined in His.pep_ch. The function:
%   1) Skips computation if an output .mat already exists (idempotent behavior).
%   2) Initializes the peptide panel (labels, derivatization map, charges, m/z, RT anchors).
%   3) Computes representative retention times (RTs), integrates intensities, and
%      builds monoisotopic traces from MS1 (with optional MS2-assisted refinement).
%   4) Serializes results to disk and draws diagnostic plots.
%   5) Optionally extracts PSM evidence when special.nDAmode==1.
%
% Inputs:
%   MS1_index   : [nMS1×2] MS1 index table; column 2 typically stores RT (minutes).
%   MS1_peaks   : MS1 peaks container/struct (format consumed by downstream functions).
%   MS2_index   : MS2 index table (used for DA/PSM/plotting).
%   MS2_peaks   : MS2 peaks container/struct.
%   ptol        : mass tolerance (ppm or Th, depending on downstream expectations).
%   cur_outpath : output folder for <out_filename>.mat and figures.
%   special     : control struct (common fields: ndebug, nDAmode, nhmass, raw_path).
%
% Disk side-effects:
%   - <cur_outpath>/H3_18_117_128.mat via output_histone().
%   - Figures via draw_layout().
%   - Optional PSM evidence via GetPSM2() when special.nDAmode==1.

% ------------------------------------ check ----------------------------------
out_filename = 'H3_18_117_128';                                 % fixed stem to name outputs consistently
% fprintf(1,'%s..',out_filename);                               % optional console log (left muted)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);        % absolute path to the expected .mat result
if 0~=exist(out_file0,'file')                                    % idempotency guard: if .mat exists already...
    return;                                                      % ...exit early to avoid recomputation
end;

% ------------------------------------- init ----------------------------------
His = init_histone();                                            % build peptide panel metadata:
                                                                 %   - mod_short/mod_type encodings
                                                                 %   - pep_ch (candidate charge states)
                                                                 %   - pep_mz (m/z per charge; computed)
                                                                 %   - rt_ref (reference RT anchors)
                                                                 %   - display flags

% ---------------------------------- calculate --------------------------------
unitdiff = 1.0032;                                               % ~C13-C12 isotopic mass difference (Da)
[pep_rts,pep_intens,mono_isointens] = ...                        % core quantification:
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...%   - RT detection & integration
                     ptol,unitdiff,His,special);                 %   - monoisotopic trace extraction

% ----------------------------------- output ----------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts); % persist His + quant matrices to <stem>.mat

% ------------------------------------ draw -----------------------------------
num_MS1 = size(MS1_index,1);                                     % number of MS1 scans (rows)
isorts = MS1_index(1:num_MS1,2);                                 % RT vector used as x-axis/sort reference
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ... % diagnostic/QC figures
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% --------------------------------- Get PSM -----------------------------------
if 1==special.nDAmode                                            % if data-assisted (DA) mode is enabled (==1)
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens, ... % extract PSM-level evidence for confidence
             isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

% ============================== LOCAL: init_histone ==========================
function His = init_histone()
%%
% Purpose:
%   Define the targeted H3 peptide set: short labels, modification encodings,
%   charge states, theoretical m/z per charge, reference RT anchors, and display flags.

His.pep_seq = 'unmod';                                           % panel descriptor (unmodified backbone label)
His.mod_short = {'VTIMPKDIQLAR';                                 % three target peptides (short labels)
    'VTIMPKDVQLAR';
    'VTIMPKEIQLAR'};
His.mod_type = {'0,pr;6,pr;';                                    % modification encodings per peptide:
    '0,pr;6,pr;';                                                %   indices refer to peptide coords (0 often N-term)
    '0,pr;6,pr;'};                                               %   'pr' typically denotes propionylation (see README)

His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);           % consider charges z = 1,2,3,4 for all peptides
%{
His.pep_mz = [900.5302  450.7687  300.8482                       % OPTIONAL example of a hard-coded m/z table
    842.4955  421.7514  281.5034
    757.4315  379.2194  253.1487];
%}
His.pep_mz = calculate_pepmz(His);                                % compute theoretical m/z per charge from seq+mods
His.rt_ref = [39.63                                              % reference RT anchors (minutes)
    17.36
    22];
His.display = ones(length(His.mod_type),1);                       % 1 = include in outputs/plots

% main ch
main_ch = His.pep_ch(1,2);                                       % choose z=2 as the "main" display charge
if main_ch~=His.pep_ch(1,1)                                      % if z=2 is not already the first column...
    [npep,ncharge] = size(His.pep_mz);                           % ...reorder columns for all peptides so that
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];         %  the main charge appears first
    x = zeros([1,ncharge]);                                      % mapping from old→new column indices
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));             % locate source column for each target slot
    end;
    tune = 1:npep;                                               % apply column permutation to all peptide rows
    His.pep_mz(tune,:) = His.pep_mz(tune,x);                     % reorder m/z columns accordingly
    His.pep_ch(tune,:) = His.pep_ch(tune,x);                     % reorder charge-ID columns accordingly
end;

% ============================ LOCAL: calculate_layout ========================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% Purpose:
%   For each peptide and for each considered charge:
%     - Determine a representative retention time (pep_rts).
%     - Integrate MS1 intensity at that RT (pep_intens).
%     - Build monoisotopic intensity traces along the RT axis (mono_isointens).
%   Also recalibrates RT anchors unless debug mode is requested.

[npep,ncharge] = size(His.pep_mz);                                % dimensions taken from m/z grid
num_MS1 = size(MS1_index,1);                                      % number of MS1 scans
pep_rts = zeros([npep,ncharge]);                                  % preallocate RTs [peptide × charge]
pep_intens = zeros([npep,ncharge]);                               % preallocate intensities [peptide × charge]
mono_isointens = zeros([num_MS1,npep]);                           % preallocate monoisotopic traces [scan × peptide]

% --------------------------- RT-anchor (rt_ref) calibration ------------------
His.rt_unmod_orig = His.rt_ref(1);                                % store original RT of the first peptide
if 1~=special.ndebug                                              % allow RT recalibration unless in debug mode
    if 2~=special.nDAmode                                         % standard confirmation path (nDAmode ≠ 2)
        for hno=1:3                                               % iterate over all three peptides in this panel
            [His.rt_ref(hno),special.ndebug] = ...
                check_ref(special.raw_path, ...
                          [His.mod_short{hno},His.mod_type{hno}], ...
                          His.rt_ref(hno), ...
                          special.ndebug);                        % verify/adjust anchors using raw_path hints
        end;
    else                                                          % DA mode (nDAmode == 2): flexible refinement
        nhmass = special.nhmass;                                  % helper mass parameter for get_rts2 (if used)
        for hno=1:3
            rt_unmod_orig = His.rt_ref(hno);                      % keep previous anchor for comparison
            His.rt_ref(hno) = check_ref(special.raw_path, ...
                                         [His.mod_short{hno},His.mod_type{hno}], ...
                                         His.rt_ref(hno), ...
                                         special.ndebug);         % attempt to confirm/update the anchor
            if rt_unmod_orig==His.rt_ref(hno)                     % unchanged → no strong evidence found
                t1 = 0;                                           % broaden search to entire RT range
                t2 = MS1_index(num_MS1,2);                        % up to the last RT observed in MS1_index
            else                                                  % updated anchor present
                delta = 5;                                        % refine within a ±5 min window around anchor
                t1 = His.rt_ref(hno)-delta;
                t2 = His.rt_ref(hno)+delta;
            end;
            [rts1,top1_rt1] = ...
                get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                         ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok % rts1 unused; top1_rt1 is the best RT
            if 0==isempty(top1_rt1)                               % if a best-RT candidate is returned
                His.rt_ref(hno) = top1_rt1;                       % adopt it as the new anchor
            end;
        end;
        special.ndebug = 1;                                       % lock anchors after DA refinement (freeze)
    end;
end;

% Legacy notes retained from original source:
% 64-69KLPFQR
% 129-134RIRGER
% 130-135IRGERA

% ----------------------------- Per-peptide extraction ------------------------
for hno=1:3                                                       % iterate over the three targeted peptides
    [cur_rts,cur_intens,cur_mono_isointens] = ...                % extract RTs, intensities, monoisotopic trace
        get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);
    if cur_rts(1)>0                                              % guard: only keep values if a valid RT was found
        pep_rts(hno,1:ncharge) = cur_rts;                        % store RTs for all considered charges
        pep_intens(hno,1:ncharge) = cur_intens;                  % store intensities for all considered charges
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;      % store monoisotopic trace (scan-wise)
    end;
end;
