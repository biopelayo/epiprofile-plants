function HH2A_01u_1_7(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ============================== OVERVIEW (EN) ================================
% Targeted quantification for an H2A peptide panel (8 targets; see init_histone),
% across charge states defined in His.pep_ch. The function:
%   1) Skips work if an output .mat already exists (idempotent behavior).
%   2) Builds the peptide panel (labels, modification map, charge states, m/z, RT anchors).
%   3) Computes representative retention times (RTs), integrates intensities, and
%      extracts monoisotopic traces from MS1 (optionally refined with MS2).
%   4) Serializes results and draws diagnostic plots.
%   5) Optionally extracts PSM evidence when special.nDAmode==1.
%
% Inputs:
%   MS1_index   : [nMS1×2] index; column 2 typically stores RT (minutes).
%   MS1_peaks   : container/struct with MS1 peaks (downstream functions know the format).
%   MS2_index   : MS2 index table (used in DA/PSM and plotting).
%   MS2_peaks   : container/struct with MS2 peaks.
%   ptol        : mass tolerance (ppm or Th; whichever is expected downstream).
%   cur_outpath : output directory for <out_filename>.mat and figures.
%   special     : control struct; common fields: ndebug, nDAmode, nhmass, raw_path.
%
% Disk side-effects:
%   - Writes <cur_outpath>/HH2A_01u_1_7.mat via output_histone().
%   - Generates figures via draw_layout().
%   - May write PSM evidence via GetPSM2() when special.nDAmode==1.

% ----------------------------------- check -----------------------------------
out_filename = 'HH2A_01u_1_7';                                  % fixed base name for outputs
% fprintf(1,'%s..',out_filename(2:end));                        % optional console log (kept muted)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);        % expected path to serialized results
if 0~=exist(out_file0,'file')                                    % idempotency: if results already exist...
    return;                                                      % ...quit early to avoid recomputation
end;

% ------------------------------------ init -----------------------------------
His = init_histone();                                            % define peptide panel metadata:
                                                                 %  - mod_short/mod_type encodings
                                                                 %  - pep_ch (candidate charge states)
                                                                 %  - pep_mz (computed from seq+mods)
                                                                 %  - rt_ref (reference RT anchors)
                                                                 %  - display flags

% --------------------------------- calculate ---------------------------------
unitdiff = 1.0032;                                               % ~C13-C12 isotopic mass difference (Da)
[pep_rts,pep_intens,mono_isointens] = ...                        % core quantification routine:
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...%  - RT detection & integration
                     ptol,unitdiff,His,special);                 %  - monoisotopic trace extraction

% ----------------------------------- output ----------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts); % persist His + quant matrices to <stem>.mat

% ------------------------------------ draw -----------------------------------
num_MS1 = size(MS1_index,1);                                     % number of MS1 scans (rows)
isorts = MS1_index(1:num_MS1,2);                                 % RT vector (x-axis / sort reference)
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ... % diagnostic/QC figures for all targets
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
%   Define the targeted H2A peptide set: short labels, modification encodings,
%   charge states, theoretical m/z per charge, reference RT anchors, display flags.

His.pep_seq = 'unmod';                                           % panel descriptor (unmodified backbone tag)
His.mod_short = {'DNKKSR';                                       % eight target peptides (short labels)
    'DNKKTR';
    'DNKKNR';
    'HLLLAIR';
    'HLQLAIR';
    'HLCLAIR';
    'HIQLAVR';
    'HVLLAVR'};
His.mod_type = {'0,pr;3,pr;4,pr;';                               % per-peptide modification encodings:
    '0,pr;3,pr;4,pr;';                                           %   indices refer to residue positions
    '0,pr;3,pr;4,pr;';                                           %   'pr' typically denotes propionylation (see README)
    '0,pr;';
    '0,pr;';
    '0,pr;';
    '0,pr;';
    '0,pr;'};

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);             % consider charges z = 1, 2, 3 for all peptides
%{
His.pep_mz = [906.5520  453.7796  302.8555                        % OPTIONAL example of hard-coded m/z values
    892.5363  446.7718  298.1836
    1457.8322 729.4197  486.6156
    1529.8719 765.4396  510.6288
    1172.7011 586.8542  391.5719
    1158.6854 579.8464  386.9000
    1412.7380 706.8726  471.5842
    1000.5574 500.7824  334.1907
    849.5417  425.2745  283.8521
    1372.8311 686.9192  458.2819
    1513.8525 757.4299  505.2890];
%}
His.pep_mz = calculate_pepmz(His);                                % compute theoretical m/z per charge from seq+mods
His.rt_ref = [34.47                                              % reference RT anchors (minutes) for each peptide
    32.31
    30.92
    30
    38.84
    22.07
    44.64
    39.95];
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
    tune = 1:npep;                                               % apply the same permutation to all peptide rows
    His.pep_mz(tune,:) = His.pep_mz(tune,x);                     % reorder m/z columns accordingly
    His.pep_ch(tune,:) = His.pep_ch(tune,x);                     % reorder charge-ID columns accordingly
end;

% ============================ LOCAL: calculate_layout ========================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% Purpose:
%   For each peptide and each considered charge:
%     - Determine a representative RT (pep_rts).
%     - Integrate MS1 intensity (pep_intens).
%     - Build the monoisotopic intensity trace along RT (mono_isointens).
%   Also recalibrates RT anchors unless debug mode is requested.

[npep,ncharge] = size(His.pep_mz);                                % dimensions derived from the m/z grid
num_MS1 = size(MS1_index,1);                                      % number of MS1 scans
pep_rts = zeros([npep,ncharge]);                                  % preallocate RTs [peptide × charge]
pep_intens = zeros([npep,ncharge]);                               % preallocate intensities [peptide × charge]
mono_isointens = zeros([num_MS1,npep]);                           % preallocate monoisotopic traces [scan × peptide]

% --------------------------- RT-anchor (rt_ref) calibration ------------------
His.rt_unmod_orig = His.rt_ref(1);                                % store original RT of the first peptide
if 1~=special.ndebug                                              % allow RT recalibration unless in debug mode
    if 2~=special.nDAmode                                         % standard confirmation path (nDAmode ≠ 2)
        for hno=1:8                                               % iterate over all eight peptides in this panel
            [His.rt_ref(hno),special.ndebug] = ...
                check_ref(special.raw_path, ...
                          [His.mod_short{hno},His.mod_type{hno}], ...
                          His.rt_ref(hno), ...
                          special.ndebug);                        % verify/adjust anchors using raw hints
        end;
    else                                                          % DA mode (nDAmode == 2): flexible refinement
        nhmass = special.nhmass;                                  % helper mass param for get_rts2 (if used)
        for hno=1:8
            rt_unmod_orig = His.rt_ref(hno);                      % keep previous anchor for comparison
            His.rt_ref(hno) = check_ref(special.raw_path, ...
                                         [His.mod_short{hno},His.mod_type{hno}], ...
                                         His.rt_ref(hno), ...
                                         special.ndebug);         % attempt to confirm/update the anchor
            if rt_unmod_orig==His.rt_ref(hno)                     % unchanged → no strong evidence found
                t1 = 0;                                           % broaden to entire RT range
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

% The following legacy notes (from the original source) keep context about
% canonical tryptic/HCD windows and neighboring peptides (not used programmatically):
% 82-88HLQLAIR
% 82-88HLQLAVR
% 4-14GGKKKSTKTSR
% 4-14SGKKKMSKLSR
% 32-39IHRHLKTR
% 32-39IHRHLKSR
% 89-99NDEELNKLLGR
% 21-29AGLQFPVGR

% ----------------------------- Per-peptide extraction ------------------------
for hno=1:8                                                       % iterate over all eight targeted peptides
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special); % extract RTs, intensities, monoisotopic trace
    if cur_rts(1)>0                                              % guard: only keep if a valid RT was found
        pep_rts(hno,1:ncharge) = cur_rts;                        % store RTs for all considered charges
        pep_intens(hno,1:ncharge) = cur_intens;                  % store intensities for all considered charges
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;      % store monoisotopic trace (scan-wise)
    end;
end;
