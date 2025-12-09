function H3_16_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ============================== OVERVIEW (EN) ================================
% Targeted quantification for a panel of two H3 peptides (init_histone.mod_short)
% across charge states defined in His.pep_ch. Steps:
%   1) Idempotency check: skip if output .mat already exists.
%   2) Initialize peptide panel (labels, mods, charge states, theoretical m/z, RT anchors).
%   3) Compute representative retention times (RTs), integrate intensities, and build
%      monoisotopic traces from MS1 (optionally refined by MS2).
%   4) Save results and draw diagnostic plots.
%   5) Optionally extract PSM evidence in data-assisted mode (special.nDAmode==1).
%
% Inputs:
%   MS1_index   : [nMS1×2], column 2 typically holds retention time (minutes).
%   MS1_peaks   : MS1 peaks container/struct (format consumed downstream).
%   MS2_index   : MS2 index table (for DA/PSM/plotting support).
%   MS2_peaks   : MS2 peaks container/struct.
%   ptol        : mass tolerance (ppm or Th depending on downstream conventions).
%   cur_outpath : output directory for <out_filename>.mat and figures.
%   special     : control struct; common fields: ndebug, nDAmode, nhmass, raw_path.
%
% Side effects on disk:
%   - <cur_outpath>/H3_16_53_63.mat (via output_histone).
%   - Figures (via draw_layout).
%   - Optional PSM evidence (via GetPSM2 when nDAmode==1).

% ----------------------------------- check -----------------------------------
out_filename = 'H3_16_53_63';                                  % fixed stem for naming outputs
% fprintf(1,'%s..',out_filename);                              % optional console log (kept muted)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);       % absolute path to expected .mat
if 0~=exist(out_file0,'file')                                   % idempotency: if .mat exists already...
    return;                                                     % ...skip the rest to avoid recomputation
end;

% ------------------------------------ init -----------------------------------
His = init_histone();                                           % define peptide panel metadata:
                                                                %  - mod_short/mod_type (derivatization map)
                                                                %  - pep_ch (candidate charge states)
                                                                %  - pep_mz (m/z per charge, computed)
                                                                %  - rt_ref (RT anchors), display flags

% --------------------------------- calculate ---------------------------------
unitdiff = 1.0032;                                              % ~C13-C12 isotopic mass difference (Da)
[pep_rts,pep_intens,mono_isointens] = ...                       % main quantification:
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...%  - RT detection + intensity integration
                     ptol,unitdiff,His,special);                %  - monoisotopic trace extraction

% ---------------------------------- output -----------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);% persist results (His + matrices) to .mat

% ----------------------------------- draw ------------------------------------
num_MS1 = size(MS1_index,1);                                    % number of MS1 scans/rows
isorts = MS1_index(1:num_MS1,2);                                % x-axis: RT vector taken from column 2
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...% diagnostic/QC figures
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------ PSM (optional) -------------------------------
if 1==special.nDAmode                                           % if data-assisted mode is enabled (==1)
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
             isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff); % extract PSM support
end;

% ============================== LOCAL: init_histone ==========================
function His = init_histone()
%%
% Purpose:
%   Define the targeted peptide set: short labels, modification encodings,
%   charge states, theoretical m/z per charge, reference RT anchors, display flags.

His.pep_seq = 'unmod';                                          % panel descriptor (unmodified backbone)
His.mod_short = {'KYQKSTELLIR';                                 % two peptide short labels
    'KYQKSTELLNR'};
His.mod_type = {'0,pr;1,pr;4,pr;';                              % modification encodings per peptide:
    '0,pr;1,pr;4,pr;'};                                         %  'pr' typically denotes propionylation (see README)

His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);          % consider charges z = 1,2,3,4 for all peptides
%{
His.pep_mz = [900.5302  450.7687  300.8482                      % OPTIONAL hard-coded m/z table (commented)
    842.4955  421.7514  281.5034
    757.4315  379.2194  253.1487];
%}
His.pep_mz = calculate_pepmz(His);                               % compute theoretical m/z per charge from seq+mods
His.rt_ref = [39.63                                             % initial reference RT anchors (minutes)
    17.36];
His.display = ones(length(His.mod_type),1);                      % 1 = include in outputs/plots

% main ch
main_ch = His.pep_ch(1,2);                                      % choose z=2 as the "main" display charge
if main_ch~=His.pep_ch(1,1)                                     % if z=2 is not in the first column...
    [npep,ncharge] = size(His.pep_mz);                          % ...reorder columns for all peptides so that
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];        %  the main charge appears first
    x = zeros([1,ncharge]);                                     % mapping from old→new column indices
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));            % locate source column for each target slot
    end;
    tune = 1:npep;                                              % apply column permutation to all peptide rows
    His.pep_mz(tune,:) = His.pep_mz(tune,x);                    % reorder m/z columns accordingly
    His.pep_ch(tune,:) = His.pep_ch(tune,x);                    % reorder charge-ID columns accordingly
end;

% ============================ LOCAL: calculate_layout ========================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% Purpose:
%   For each peptide and charge:
%     - Determine a representative RT (pep_rts).
%     - Integrate MS1 intensity (pep_intens).
%     - Build monoisotopic intensity traces along RT (mono_isointens).
%   Also recalibrate RT anchors unless debug mode is requested.

[npep,ncharge] = size(His.pep_mz);                               % dimensions taken from m/z grid
num_MS1 = size(MS1_index,1);                                     % number of MS1 scans
pep_rts = zeros([npep,ncharge]);                                 % preallocate RTs [peptide × charge]
pep_intens = zeros([npep,ncharge]);                              % preallocate intensities [peptide × charge]
mono_isointens = zeros([num_MS1,npep]);                          % preallocate monoisotopic traces [scan × peptide]

% -------------------------- RT-anchor (rt_ref) calibration -------------------
His.rt_unmod_orig = His.rt_ref(1);                               % store original RT of first peptide
if 1~=special.ndebug                                             % allow RT recalibration unless in debug
    if 2~=special.nDAmode                                        % standard confirm (nDAmode≠2)
        for hno=1:2                                              % loop over both peptides in the panel
            [His.rt_ref(hno),special.ndebug] = ...
                check_ref(special.raw_path, ...
                          [His.mod_short{hno},His.mod_type{hno}], ...
                          His.rt_ref(hno), ...
                          special.ndebug);                       % verify/adjust anchors using raw_path
        end;
    else                                                         % DA mode (nDAmode==2): flexible refinement
        nhmass = special.nhmass;                                 % helper mass for get_rts2 (if used)
        for hno=1:2
            rt_unmod_orig = His.rt_ref(hno);                     % keep previous anchor
            His.rt_ref(hno) = check_ref(special.raw_path, ...
                                         [His.mod_short{hno},His.mod_type{hno}], ...
                                         His.rt_ref(hno), ...
                                         special.ndebug);        % try to confirm/update anchor
            if rt_unmod_orig==His.rt_ref(hno)                    % unchanged → no clear evidence found
                t1 = 0;                                          % search whole RT range
                t2 = MS1_index(num_MS1,2);                       % from 0 to last RT in MS1_index
            else                                                 % updated anchor
                delta = 5;                                       % refine in a ±5 min window around anchor
                t1 = His.rt_ref(hno)-delta;
                t2 = His.rt_ref(hno)+delta;
            end;
            [rts1,top1_rt1] = ...
                get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                         ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok % rts1 unused; top1_rt1 is the best RT
            if 0==isempty(top1_rt1)                              % if a best-RT candidate is returned
                His.rt_ref(hno) = top1_rt1;                      % adopt it as the new anchor
            end;
        end;
        special.ndebug = 1;                                      % lock anchors after DA refinement (no more changes)
    end;
end;

% Legacy notes (kept from original source):
% 64-69KLPFQR
% 129-134RIRGER
% 130-135IRGERA

% ---------------------------- Per-peptide extraction -------------------------
for hno=1:2                                                      % iterate over each peptide
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special); % extract RTs, intensities, trace
    if cur_rts(1)>0                                             % guard: only keep if RT is valid
        pep_rts(hno,1:ncharge) = cur_rts;                       % store RTs for all charges
        pep_intens(hno,1:ncharge) = cur_intens;                 % store intensities for all charges
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;     % store monoisotopic trace (scan-wise)
    end;
end;
