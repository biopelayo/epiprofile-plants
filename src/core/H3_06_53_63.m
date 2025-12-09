function H3_06_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%% H3_06_53_63 — Orchestrates the quant/plot pipeline for H3 peptide 53–63 (K56 PTM panel)
% PURPOSE
%   This wrapper prepares a histone peptide panel (H3 53–63) with multiple K56
%   modification states (unmodified, me1, me2, me3, ac), computes their m/z
%   for several charge states, calibrates RT using the unmodified form, and
%   quantifies each PTM variant. It also writes outputs and figures and, in
%   DDA-assist mode, extracts PSM-level evidence.
%
% INPUTS
%   MS1_index  : [N x 2] (or similar) index table for MS1 scans; col 2 must hold RTs (minutes).
%   MS1_peaks  : structure/array of MS1 peaks compatible with GetProfiles/get_area downstream.
%   MS2_index  : index table for MS2 scans (used by get_rts2 when special.nDAmode==2).
%   MS2_peaks  : container with MS2 peaks/fragments (used in DDA-assisted workflows).
%   ptol       : mass tolerance for matching (unit depends on called functions; legacy behavior preserved).
%   cur_outpath: filesystem path for outputs (MAT-files, figures).
%   special    : configuration struct (fields used: nDAmode, ndebug, raw_path, nhmass).
%
% SIDE EFFECTS
%   - Creates a MAT file with quantification for this panel.
%   - Writes PDF figures (layout/heatmap) and a "Figure Legends.txt".
%   - May extract PSMs when DDA-assist mode is enabled (special.nDAmode==1).
%
% DEPENDENCIES (external to this file)
%   GetMods, calculate_layout, output_histone, draw_layout, GetPSM,
%   calculate_pepmz, check_ref, get_rts2, get_histone0, relocateD,
%   relocate, relocate2, get_histone1, get_rts, find_pair.

%%

% check
out_filename = 'H3_06_53_63';                    % Fixed panel identifier (used for file naming).
fprintf(1,'%s..',out_filename);                  % Console progress marker.
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']); % Target .mat result for this panel.
if 0~=exist(out_file0,'file')                    % If output already exists…
    return;                                      % …skip computation to be idempotent / fast.
end;

% init
His = init_histone(cur_outpath,out_filename);    % Build peptide library for this panel (seq, PTMs, charges, m/z, RT refs).

% calculate
unitdiff = 1.0032;                               % Nominal 13C–12C mass difference (~1.0032 Da).
Mods = GetMods();                                % Load modification mass table (external).
% Core quantification: computes RTs and intensities for all PTM variants,
% possibly using both MS1 and MS2 (get_rts2) when DDA/nhmass configured in special.
[pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

% output
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts); % Persist peptide-level results (AUC, RT).

% draw
num_MS1 = size(MS1_index,1);                     % Number of MS1 scans for building time axis.
isorts = MS1_index(1:num_MS1,2);                 % Extract RTs (assumed sorted/monotonic).
% Render figures for layout and optionally overlay MS2 info depending on 'special'.
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS2_index,MS2_peaks,special);

% Get PSM
% Optional: if DA-mode is enabled (value==1 here), collect PSM-level evidence for the panel.
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

function His = init_histone(cur_outpath,out_filename)
%% INIT_HISTONE — Define sequence/PTM states/charges for H3 53–63 (K56 panel)
% OUTPUT
%   His: structure containing:
%        - pep_seq    : base peptide sequence (string).
%        - mod_short  : cellstr of modification labels (display-friendly).
%        - mod_type   : cellstr PTM encodings (position,type) e.g. '4,me1;'.
%        - pep_ch     : [nPTM x nCharge] charge states for each PTM variant.
%        - pep_mz     : [nPTM x nCharge] computed precursor m/z for each (PTM,charge).
%        - rt_ref     : column vector of library RTs (min), one per PTM state.
%        - display    : mask to mark visible items in plots/tables.
%        - outpath, outfile: for downstream I/O.
%
% NOTES
%   - The PTM encoding uses indexes relative to the peptide substring, not full protein numbering.
%   - 'pr' denotes propionylation (derivatization) used in histone workflows.
%   - Charges (1–4) are provisioned for robustness; main channel is enforced below.

%%

His.pep_seq = 'KYQKSTELLIR';                     % H3 53–63 peptide (K56 within; 0-based PTM index 4 maps to K56 here).
His.mod_short = {'unmod';                         % Human-readable labels for each state:
    'K56me1';                                     %   unmodified, mono/di/tri-methylation at K56, and acetylation.
    'K56me2';
    'K56me3';
    'K56ac'};
His.mod_type = {'0,pr;1,pr;4,pr;';                % PTM grammar: 'pos,type;' where pos is 0-based within pep_seq.
    '0,pr;1,pr;4,me1;';                           % 'pr' marks propionylation at N-terminus and Lys positions.
    '0,pr;1,pr;4,me2;';
    '0,pr;1,pr;4,me3;';
    '0,pr;1,pr;4,ac;'};

His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1); % Allow charges z=1..4 for each PTM row (same z-grid across states).
%{
His.pep_mz = [1362.7627	681.885	454.9258	341.4461
    1376.7784	688.8928	459.5976	344.9501
    1334.7678	667.8876	445.5941	334.4474
    1348.7835	674.8954	450.2660	337.9513
    1348.7471	674.8772	450.2539	337.9422
    1418.7890	709.8981	473.6012	355.4527
    1418.7890	709.8981	473.6012	355.4527];
%}
His.pep_mz = calculate_pepmz(His);               % Compute m/z from sequence+PTMs+charges using current mass tables.
His.rt_ref = [40.37                                % Library RTs (min) for each PTM state (order aligned with mod lists).
    41.20
    33.92
    33.72
    39.17];
His.display = ones(length(His.mod_type),1);      % Flag all rows for display by default.

His.outpath = cur_outpath;                       % Persist paths/names for downstream helpers.
His.outfile = out_filename;

% main ch
% Ensure the "main" display/processing channel is the 2+ charge (historically robust for histone peptides).
main_ch = His.pep_ch(1,2);                       % Desired main charge = second column entry in first row (z=2).
if main_ch~=His.pep_ch(1,1)                      % If current first column is not the main charge…
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)]; % Build new column order putting main_ch first.
    x = zeros([1,ncharge]);                      % x maps old column positions -> new order.
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = 1:npep;                               % Reorder columns consistently for all PTM rows.
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end;

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special)
%% CALCULATE_LAYOUT — Core quantification logic for this H3 K56 panel
% OUTPUTS
%   pep_rts         : [nPTM x nCharge] selected RT for each (PTM,charge).
%   pep_intens      : [nPTM x nCharge] integrated AUC for each (PTM,charge).
%   mono_isointens  : [num_MS1 x nPTM] monoisotopic traces (first charge) for QC/plotting.
%
% STRATEGY
%   1) Calibrate RT using the unmodified species (row 1) — possibly refine via MS2 (DDA) if enabled.
%   2) Shift all PTM RT references by the measured delta (align panel).
%   3) For debug/DA modes, run relocation helpers to refine expected RTs.
%   4) Quantify each PTM variant with get_histone1 (MS1 profiles) using updated RT anchors.

%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);                 % Will store per-PTM, per-charge RTs.
pep_intens = zeros([npep,ncharge]);              % Will store per-PTM, per-charge integrated areas.
mono_isointens = zeros([num_MS1,npep]);          % MS1 monoisotopic trace kept per PTM (first charge only).

% unmod
His.rt_unmod_orig = His.rt_ref(1);               % Save original library RT for unmodified state.
if 1~=special.ndebug                              % If not in debug mode, attempt RT refinement.
    if 2~=special.nDAmode                         % Non-DA mode: check/adjust RT from external prior (file-based).
        [His.rt_ref(1),special.ndebug] = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
    else                                          % DA mode: allow MS2-informed RT refinement.
        nhmass = special.nhmass;                  % Neutral loss / mass setting for DA search (provided in special).
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)       % If check_ref did not change RT…
            t1 = 0;
            t2 = MS1_index(num_MS1,2);           % …search across the full gradient to find a better apex.
        else
            delta = 5;
            t1 = His.rt_ref(1)-delta;            % …or search in a ±5 min neighborhood around refined RT.
            t2 = His.rt_ref(1)+delta;
        end;
        hno = 1;% unmod
        % get_rts2 can leverage MS2 evidence (fragments) to sharpen the apex choice.
        [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,1,t1,t2,nhmass);%#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;            % Adopt refined apex RT if available.
        end;
    end;
end;

hno = 1;                                         % Index for unmodified state in this panel.
% get_histone0 extracts monoisotopic profiles, chooses apex near RT_ref, and integrates area per charge.
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% calibrate the rt_ref
if cur_rts(1)>0                                  % If the main charge produced a valid apex…
    His.rt_ref(1) = cur_rts(1);                  % …take its RT as new anchor,
    delta = cur_rts(1)-His.rt_unmod_orig;        % compute the shift w.r.t. the original RT,
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta; % and propagate the RT shift to all PTM states.
    pep_rts(hno,1:ncharge) = cur_rts;            % Persist measured RTs for the unmodified row.
    pep_intens(hno,1:ncharge) = cur_intens;      % Persist intensities for unmodified.
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens; % Store the EIC for QC/plots.
end;
% After calibration, (re)locate expected RTs depending on mode:
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);       % Debug relocation (diagnostic).
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);    % MS1-only relocation (rules below).
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass); % MS2-assisted relocation.
    end;
end;

% K56me1
% K56me2
% K56me3
% K56ac
for hno=2:5                                       % Loop over the 4 modified states following unmodified=1.
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0                               % Accept only if apex was found for main charge.
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%% RELOCATE — Heuristic RT re-anchoring for PTM states using MS1-only evidence
% For each PTM state, we define search windows relative to the unmodified RT,
% pick apex candidates with get_rts, and assign refined RT_ref values.
% This helps when PTM states elute systematically earlier/later.

%%

delta = 0.1;                                     % Narrow offset used for window definitions (minutes).
nsplit = 1;                                      % Legacy argument for get_rts (segmenting EIC if needed).

% K56me1 — search after the unmodified (later elution expected by rule)
hno = 2;
t1 = His.rt_ref(1)+delta;                        % Start just after unmod apex.
t2 = His.rt_ref(1)+14;                           % Allow up to +14 min (panel-specific heuristic).
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;                          % No apex found: mark as zero (downstream will skip quant).
else
    His.rt_ref(hno) = top1_rt2;                   % Use the top RT candidate.
end;

% K56me2 — search earlier than the unmodified
hno = 3;
t1 = 6;                                          % Lower absolute bound (min 6 min).
t2 = His.rt_ref(1)-3;                            % Up to (unmod - 3min).
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K56me3 — similarly earlier
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% Pairing rule: choose a consistent me2/me3 pair using both RT and area evidence.
[His.rt_ref(3),His.rt_ref(4)] = find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3,rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K56ac — search between unmod and me2 (typical acetyl elution rule)
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;            % Midpoint between unmod and me2.
t2 = His.rt_ref(1)-delta;                        % Up to just before unmod.
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;                          % Not found → mark zero.
else
    His.rt_ref(hno) = top1_rt5;                   % Assign RT of the best candidate.
end;

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%% RELOCATE2 — Same idea as relocate but leveraging MS2 (DDA) via get_rts2
% This variant uses fragment evidence (when available) to score/confirm
% apex candidates, increasing robustness when MS1 signal is noisy.

%%

delta = 0.1;
nsplit = 1;

% K56me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K56me2
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K56me3
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% Select consistent me2/me3 pair considering intensities and RT spacing.
[His.rt_ref(3),His.rt_ref(4)] = find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3,rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K56ac
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;
