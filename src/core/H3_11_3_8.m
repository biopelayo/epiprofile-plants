function H3_11_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% === Purpose (English) =======================================================
% Entry point to quantify and visualize three H3 tryptic/propionylated peptides
% (short labels in init_histone.mod_short) across available charge states.
% It: (1) checks if results already exist, (2) initializes histone metadata,
% (3) computes peptide RTs and intensities from MS1/MS2, (4) writes outputs,
% (5) draws diagnostic plots, and (6) optionally extracts PSM-level evidence
% in data-assisted mode (special.nDAmode==1).
%
% Inputs:
%   MS1_index  : [nMS1×2] index table for MS1 scans (col2 usually retention time).
%   MS1_peaks  : MS1 centroid/peak list container (format handled downstream).
%   MS2_index  : MS2 index table (used for MS2-based checks/validation).
%   MS2_peaks  : MS2 peaks container.
%   ptol       : mass tolerance (ppm or Th depending on downstream functions).
%   cur_outpath: output folder for .mat and figures.
%   special    : struct of flags/paths (e.g., ndebug, nDAmode, nhmass, raw_path).
%
% Outputs (on disk):
%   - <cur_outpath>/H3_11_3_8.mat with the His structure and quantification.
%   - Figures via draw_layout().
%   - Optional PSM evidence via GetPSM2() when nDAmode==1.

% check
out_filename = 'H3_11_3_8';                                 % fixed stem for output naming
% fprintf(1,'%s..',out_filename);                            % optional console log (muted)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);    % expected output .mat path
if 0~=exist(out_file0,'file')                                % if results already exist…
    return;                                                  % …skip work to save time
end;

% init
His = init_histone();                                        % build metadata for H3 peptides:
                                                             % sequences, mods, charge states,
                                                             % m/z, reference RTs, display flags

% calculate
unitdiff = 1.0032;                                           % ~C13-C12 isotopic mass difference
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,...
                     ptol,unitdiff,His,special);             % main quantification routine:
                                                             % RT detection + intensity integration
                                                             % + monoisotopic profiles

% output
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts); % persist quant results to .mat

% draw
num_MS1 = size(MS1_index,1);                                 % count of MS1 scans/rows
isorts = MS1_index(1:num_MS1,2);                             % extract RT column as sort/reference
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,...
            isorts,mono_isointens,MS2_index,MS2_peaks,special); % produce diagnostic plots

% Get PSM
if 1==special.nDAmode                                         % if data-assisted mode is enabled
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens,...
             isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff); % get PSM evidence
end;

% ============================================================================
% Local function: Define the histone-peptide panel and static metadata
% ============================================================================
function His = init_histone()
%%
% === Purpose ================================================================
% Define the peptide list, modification encodings, charge states, m/z,
% reference RTs and display flags for this H3 panel. This is the “spec”
% that the rest of the pipeline uses to look for signals.

His.pep_seq = 'unmod';                                       % panel label (here: unmodified backbone)
His.mod_short = {'TKQTAR';                                   % short peptide labels (3 entries)
    'TKQSAR';
    'SNQTAR'};
His.mod_type = {'0,pr;2,pr;';                                % per-peptide modification encoding:
    '0,pr;2,pr;';                                            % 'pr' typically stands for propionylation
    '0,pr;'};                                                % (N-terminus and Lys residues derivatized)

His.pep_ch = repmat([1 2],length(His.mod_type),1);           % track charge states per peptide (z=1,2)
%{
His.pep_mz = [900.5302	450.7687	300.8482                   % OPTIONAL hard-coded m/z table (commented):
    842.4955	421.7514	281.5034                            % columns correspond to charge states
    757.4315	379.2194	253.1487];
%}
His.pep_mz = calculate_pepmz(His);                            % compute m/z from sequence+mods+charges
His.rt_ref = [39.63                                           % initial reference RTs (minutes)
    17.36
    22];
His.display = ones(length(His.mod_type),1);                   % 1=display peptide in plots/tables

% main ch
main_ch = His.pep_ch(1,1);                                   % set "main" charge as the first of row1
if main_ch~=His.pep_ch(1,1)                                  % safeguard (no-op with current init)
    [npep,ncharge] = size(His.pep_mz);                       % if needed, re-order columns so that
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];     % the first column is the main charge
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));         % map old → new column order
    end;
    tune = 1:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);                 % reorder m/z columns
    His.pep_ch(tune,:) = His.pep_ch(tune,x);                 % reorder charge columns
end;

% ============================================================================
% Local function: Quantify RT/intensity profiles for the peptide panel
% ============================================================================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% === Purpose ================================================================
% For each peptide and charge state listed in His, derive:
%  - pep_rts(hno, ch)       : representative RTs (per peptide and charge)
%  - pep_intens(hno, ch)    : integrated intensities (per peptide and charge)
%  - mono_isointens(scan,hno): monoisotopic MS1 intensity trace per peptide
%
% Strategy:
%   1) Initialize output matrices sized by number of peptides (npep),
%      number of charge states (ncharge), and number of MS1 scans (num_MS1).
%   2) Calibrate reference RTs (His.rt_ref) by searching near prior anchors
%      or scanning the full RT range (when needed / DA mode).
%   3) For each peptide, call get_histone0() to extract RTs/intensities and
%      the monoisotopic profile within tolerance ptol and isotopic spacing
%      unitdiff, using MS1/MS2 info as available.

[npep,ncharge] = size(His.pep_mz);                            % matrix size from precomputed m/z table
num_MS1 = size(MS1_index,1);                                  % number of MS1 scans
pep_rts = zeros([npep,ncharge]);                              % preallocate RTs (peptide × charge)
pep_intens = zeros([npep,ncharge]);                           % preallocate intensities
mono_isointens = zeros([num_MS1,npep]);                       % MS1 monoisotopic traces per peptide

% calibrate the rt_ref
His.rt_unmod_orig = His.rt_ref(1);                            % keep a copy of initial RT for peptide #1
if 1~=special.ndebug                                          % allow RT recalibration unless in debug
    if 2~=special.nDAmode                                     % standard mode (not “2”/DA scan)
        for hno=1:3                                           % three peptides defined in this panel
            [His.rt_ref(hno),special.ndebug] = ...
                check_ref(special.raw_path, ...
                          [His.mod_short{hno},His.mod_type{hno}], ...
                          His.rt_ref(hno), ...
                          special.ndebug);                    % validate/adjust RT anchors from raw
        end;
    else                                                      % DA mode: refine anchors more flexibly
        nhmass = special.nhmass;                              % neutral loss / helper mass (DA tools)
        for hno=1:3
            rt_unmod_orig = His.rt_ref(hno);                  % store current anchor
            His.rt_ref(hno) = check_ref(special.raw_path, ...
                                         [His.mod_short{hno},His.mod_type{hno}], ...
                                         His.rt_ref(hno), ...
                                         special.ndebug);     % try to confirm anchor from raw
            if rt_unmod_orig==His.rt_ref(hno)                 % if unchanged (no strong evidence)…
                t1 = 0;                                       % …search broadly: start at 0 min
                t2 = MS1_index(num_MS1,2);                    % …end at last RT from MS1 index
            else                                              % if adjusted anchor was found…
                delta = 5;                                    % …search in a ±5 min window
                t1 = His.rt_ref(hno)-delta;
                t2 = His.rt_ref(hno)+delta;
            end;
            [rts1,top1_rt1] = ...
                get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                         ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok  % refine with RT profiling
            if 0==isempty(top1_rt1)                           % if a top-RT candidate exists
                His.rt_ref(hno) = top1_rt1;                   % update the anchor to the best RT
            end;
        end;
        special.ndebug = 1;                                   % lock anchors after DA refinement
    end;
end;

% 64-69KLPFQR
% 129-134RIRGER
% 130-135IRGERA
% (Comments above: historical peptide indexing notes for H3 fragments.)

for hno=1:3                                                   % iterate over the three H3 peptides
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special); % extract per-peptide traces
    if cur_rts(1)>0                                           % if a valid RT was found (guard)
        pep_rts(hno,1:ncharge) = cur_rts;                     % store RTs for all considered charges
        pep_intens(hno,1:ncharge) = cur_intens;               % store integrated intensities
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;   % store monoisotopic MS1 profile
    end;
end;
