function H3_14_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% === Purpose (English) =======================================================
% Targeted quantification for a fixed panel of three H3 peptides (see
% init_histone.mod_short) across charge states defined in His.pep_ch.
% Workflow:
%   1) Skip if results already exist (idempotent execution).
%   2) Initialize the histone-peptide panel (seq/mods/charges/mz/RT anchors).
%   3) Compute per-peptide retention times (RTs), intensities, and mono-isotopic
%      traces from MS1 (with optional MS2-assisted refinements).
%   4) Save results to disk and draw diagnostic figures.
%   5) Optionally extract PSM-level evidence when special.nDAmode==1.
%
% Inputs:
%   MS1_index  : [nMS1×2] MS1 index; column 2 typically holds RT (minutes).
%   MS1_peaks  : container/struct of MS1 peaks (format handled downstream).
%   MS2_index  : MS2 index table (used in plotting/DA/PSM routines).
%   MS2_peaks  : container/struct of MS2 peaks.
%   ptol       : mass tolerance (ppm or Th, depending on downstream expectations).
%   cur_outpath: output directory for <out_filename>.mat and figures.
%   special    : control struct (common fields: ndebug, nDAmode, nhmass, raw_path).
%
% Disk outputs:
%   - <cur_outpath>/H3_14_27_40.mat via output_histone().
%   - Figures via draw_layout().
%   - Optional PSM evidence via GetPSM2() if special.nDAmode==1.

% ------------------------------- check --------------------------------------
out_filename = 'H3_14_27_40';                                  % fixed stem for output naming
% fprintf(1,'%s..',out_filename);                              % optional console log (muted)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);       % expected path to the .mat result
if 0~=exist(out_file0,'file')                                   % if the result already exists...
    return;                                                     % ...exit early (idempotent behavior)
end;

% ------------------------------- init ---------------------------------------
His = init_histone();                                           % build peptide panel metadata:
                                                                %   - mod_short/mod_type encodings
                                                                %   - pep_ch (charge states)
                                                                %   - pep_mz (m/z per charge, computed)
                                                                %   - rt_ref (reference RT anchors)
                                                                %   - display flags

% ----------------------------- calculate ------------------------------------
unitdiff = 1.0032;                                              % ~C13-C12 isotopic mass difference (Da)
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);                % main quantification:
                                                                %   - RT detection/integration
                                                                %   - mono-isotopic traces

% ------------------------------ output --------------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);% persist His + quant matrices to .mat

% ------------------------------- draw ---------------------------------------
num_MS1 = size(MS1_index,1);                                    % number of MS1 scans/rows
isorts = MS1_index(1:num_MS1,2);                                % RT vector used as x-axis/sort reference
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...% diagnostic plots/QC panels
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------ PSM (optional) ------------------------------
if 1==special.nDAmode                                           % if data-assisted mode is enabled (==1)
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
             isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff); % extract PSM evidence
end;

% ============================= LOCAL: init_histone ===========================
function His = init_histone()
%%
% Purpose:
%   Define the targeted peptide set: short labels, modification encodings,
%   charge states, theoretical m/z values per charge, reference RT anchors,
%   and display flags.

His.pep_seq = 'unmod';                                          % panel descriptor (unmodified backbone)
His.mod_short = {'KSAPATGGVKKPHR';                              % three short peptide labels
    'KSAPTTGGVKKPHR';
    'QSAPATGGVKKPHR'};
His.mod_type = {'0,pr;1,pr;10,pr;11,pr;';                       % modification encoding per peptide:
    '0,pr;1,pr;10,pr;11,pr;';                                   %  indices use peptide coords (0 often N-term)
    '0,pr;10,pr;11,pr;'};                                       %  'pr' typically denotes propionylation (see README)

His.pep_ch = repmat([2 3 4],length(His.mod_type),1);            % charges considered: z = 2, 3, 4
%{
His.pep_mz = [900.5302  450.7687  300.8482                      % OPTIONAL example m/z table (commented)
    842.4955  421.7514  281.5034
    757.4315  379.2194  253.1487];
%}
His.pep_mz = calculate_pepmz(His);                               % compute theoretical m/z per charge from seq+mods
His.rt_ref = [39.63                                             % reference RT anchors (minutes)
    17.36
    22];
His.display = ones(length(His.mod_type),1);                      % 1 = display peptide in outputs/plots

% main ch
main_ch = His.pep_ch(1,1);                                      % choose column 1 as "main" charge
if main_ch~=His.pep_ch(1,1)                                     % condition is a no-op as written (always false);
    [npep,ncharge] = size(His.pep_mz);                          % if logic changes in the future, this block
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];        % would reorder columns to place main_ch first
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = 1:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end;

% =========================== LOCAL: calculate_layout =========================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% Purpose:
%   For each peptide and charge state:
%     - Determine a representative RT (pep_rts).
%     - Integrate intensity (pep_intens).
%     - Build mono-isotopic intensity traces across the RT grid (mono_isointens).
%   Recalibrate RT anchors unless debug mode is on.

[npep,ncharge] = size(His.pep_mz);                               % dimensions from the m/z table
num_MS1 = size(MS1_index,1);                                     % number of MS1 scans/rows
pep_rts = zeros([npep,ncharge]);                                 % preallocate RT matrix
pep_intens = zeros([npep,ncharge]);                              % preallocate intensity matrix
mono_isointens = zeros([num_MS1,npep]);                          % preallocate mono-isotopic traces

% ------------------------ RT-anchor (rt_ref) calibration --------------------
His.rt_unmod_orig = His.rt_ref(1);                               % store original RT of peptide #1
if 1~=special.ndebug                                             % allow recalibration unless in debug mode
    if 2~=special.nDAmode                                        % standard confirmation (nDAmode ~= 2)
        for hno=1:3                                              % loop over the three peptides in this panel
            [His.rt_ref(hno),special.ndebug] = ...
                check_ref(special.raw_path, ...
                          [His.mod_short{hno},His.mod_type{hno}], ...
                          His.rt_ref(hno), ...
                          special.ndebug);                       % verify/adjust RT anchors using raw path
        end;
    else                                                         % DA mode (nDAmode == 2): flexible refinement
        nhmass = special.nhmass;                                 % helper mass parameter for get_rts2 (if used)
        for hno=1:3
            rt_unmod_orig = His.rt_ref(hno);                     % keep current anchor
            His.rt_ref(hno) = check_ref(special.raw_path, ...
                                         [His.mod_short{hno},His.mod_type{hno}], ...
                                         His.rt_ref(hno), ...
                                         special.ndebug);        % try to confirm/update the anchor
            if rt_unmod_orig==His.rt_ref(hno)                    % unchanged → no clear evidence
                t1 = 0;                                          % search entire RT range...
                t2 = MS1_index(num_MS1,2);                       % ...from 0 to the last RT in MS1_index
            else                                                 % updated anchor found
                delta = 5;                                       % refine search within ±5 minutes around anchor
                t1 = His.rt_ref(hno)-delta;
                t2 = His.rt_ref(hno)+delta;
            end;
            [rts1,top1_rt1] = ...
                get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                         ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok  % rts1 unused; top1_rt1 provides best RT
            if 0==isempty(top1_rt1)                              % if a best-RT candidate is available
                His.rt_ref(hno) = top1_rt1;                      % adopt it as the new anchor
            end;
        end;
        special.ndebug = 1;                                      % lock anchors after DA-mode refinement
    end;
end;

% Legacy notes retained from source:
% 64-69KLPFQR
% 129-134RIRGER
% 130-135IRGERA

% --------------------------- Per-peptide extraction --------------------------
for hno=1:3                                                      % iterate over all 3 peptides
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special); % extract RTs, intensities, trace
    if cur_rts(1)>0                                             % guard: only store if RT is valid
        pep_rts(hno,1:ncharge) = cur_rts;                       % store RTs for all considered charges
        pep_intens(hno,1:ncharge) = cur_intens;                 % store intensities per charge
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;     % store mono-isotopic trace (scan-wise)
    end;
end;
