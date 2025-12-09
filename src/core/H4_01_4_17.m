function H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% =========================================================================
% H4_01_4_17  —  Quantification panel for Histone H4 (aa 4–17) N-terminal peptide
% =========================================================================
%
% PURPOSE
%   Extract and quantify XICs for the H4 N-terminal peptide:
%       GKGGKGLGKGGAKR   (residues 4–17 in histone H4)
%   across its acetylation states: unmodified (chemically propionylated) and
%   all combinations of K5ac, K8ac, K12ac, K16ac (mono-, di-, tri-, tetra-).
%
% HOW THIS PANEL FITS IN
%   This is a building block in an EpiProfile-style targeted workflow for
%   histone PTMs. It assumes chemical propionylation (“pr”) in sample prep.
%   The unmodified derivatized peptide (“unmod”) is treated as the RT anchor.
%   We then relocate expected windows for modified forms around that anchor.
%
% INPUTS
%   MS1_index   : [nMS1 x 2] matrix, where column 1 is scan index (or any
%                 monotonically increasing identifier) and column 2 is
%                 retention time (in minutes).
%   MS1_peaks   : Core-specific container of centroided MS1 peaks
%                 (whatever your EpiProfile-like core expects).
%   MS2_index   : [nMS2 x 2] index-RT matrix for MS2 scans (used in DA mode).
%   MS2_peaks   : Core-specific container of centroided MS2 peaks (DA mode).
%   ptol        : Mass tolerance passed through to extraction routines. In
%                 ppm for Orbitrap/TOF-style centroided data is typical.
%   cur_outpath : Output directory where this panel writes .mat and figures.
%   special     : struct of runtime switches and hints:
%                 - ndebug   : if 1, diagnostic relocation (relocateD)
%                              with minimal logic (useful for debugging).
%                 - nDAmode  : relocation/PSM semantics (historical):
%                               * 2 → use MS2-assisted relocation (relocate2/get_rts2)
%                               * ≠2 → MS1-only relocation (relocate/get_rts)
%                               * 1 → additionally export PSM after quant
%                 - nhmass   : neutral-loss/precursor mass hint for get_rts2
%                               (if your core supports it; can be 0 if n/a).
%                 - raw_path : vendor path for check_ref (anchors unmod RT).
%
% OUTPUTS (written to disk)
%   <cur_outpath>/H4_01_4_17.mat  containing:
%     - His             : struct with peptide/charge/PTM definitions and RT refs
%     - pep_rts         : [16 x 4] RT per modification row x charge-state column
%     - pep_intens      : [16 x 4] integrated XIC intensities (aligned to pep_rts)
%     - mono_isointens  : [nMS1 x 16] monoisotopic XIC trace per modification row
%   plus figures from draw_layout(...) (format depends on your plotting code).
%   If special.nDAmode==1, the function also triggers GetPSM(...) to export
%   PSM tables from your core.
%
% IMPORTANT CONVENTIONS
%   - “pr” inside mod_type denotes chemical propionylation from sample prep,
%     not a biological PTM.
%   - Indexing inside mod_type is relative to the peptide, not to the full
%     histone sequence. For this panel, K positions map to local indices
%     2, 5, 9, 13 (see init_histone for details). Keep this consistent with
%     your EpiProfile-like core; do not mix mapping schemes between panels.
%
% DEPENDENCIES (must exist in your codebase)
%   calculate_pepmz, GetMods, output_histone, draw_layout, GetPSM,
%   check_ref, relocateD, get_rts, get_rts2, get_histone0, get_histone1,
%   get_histone4, get_histone6  (and any helpers those rely on).
%
% NOTE FOR PROTEOMICS USERS (NON-BIOINFORMATICS)
%   You do not need to edit the code to run the panel. You only need to:
%     1) Provide centroided MS1(/MS2) data as expected by your core,
%     2) Set ptol appropriately (e.g., 5–10 ppm Orbitrap; 10–20 ppm TOF),
%     3) Pass special.nDAmode=2 if you want RT relocation to leverage MS2,
%     4) Ensure raw_path is valid so check_ref can anchor the unmodified RT.
%   Everything else (m/z lists, charge states, RT heuristics) is encoded here.
%
% =========================================================================

% ----------------------------
% 1) Idempotency guard: if this panel was already computed, return.
% ----------------------------
out_filename = 'H4_01_4_17';
fprintf(1,'%s..',out_filename);           % progress print (non-blocking)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    % If output exists, we do not recompute. This keeps runs fast and
    % avoids overwriting results unless you delete the .mat yourself.
    return;
end;

% ----------------------------
% 2) Define peptide/PTMs/charges/seed RTs inside His struct
% ----------------------------
His = init_histone(cur_outpath,out_filename);

% ----------------------------
% 3) Core quantification and RT relocation
% ----------------------------
unitdiff = 1.0032;                       % 13C isotopic spacing for XIC tracing
Mods = GetMods();                        % PTM mass table (used by get_histone4/6)
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,Mods,His,special);

% The tri-acetyl block is reordered for readability in figures/reports.
[His,pep_rts,pep_intens,mono_isointens] = ...
    change_order(His,pep_rts,pep_intens,mono_isointens);

% ----------------------------
% 4) Persist results and draw summary plots
% ----------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% Prepare “isorts” (the RT vector) for plotting utilities
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);

% Plot panel layout and per-row XICs according to your core’s style
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ----------------------------
% 5) Optional PSM export (historical switch: nDAmode==1)
% ----------------------------
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

end % <-- end main function


% =========================================================================
% NESTED: init_histone
%   Build the “His” struct with peptide sequence, PTM rows, charge grid,
%   theoretical m/z per charge, seed RT values, and plotting flags.
% =========================================================================
function His = init_histone(cur_outpath,out_filename)
%%

% Peptide sequence (H4 aa 4–17)
His.pep_seq = 'GKGGKGLGKGGAKR';

% Rows in the panel (16 total), from unmodified to tetra-acetylated.
% mod_short is a human-readable short name used in figures/reports.
His.mod_short = {'unmod';
    'K5ac';
    'K8ac';
    'K12ac';
    'K16ac';
    'K5acK8ac';
    'K5acK12ac';
    'K5acK16ac';
    'K8acK12ac';
    'K8acK16ac';
    'K12acK16ac';
    'K8acK12acK16ac';
    'K5acK12acK16ac';
    'K5acK8acK16ac';
    'K5acK8acK12ac';
    'K5acK8acK12acK16ac'};

% mod_type encodes, for each row, which positions have “pr” vs “ac”.
% Indexing is RELATIVE TO THE PEPTIDE SEQUENCE, not full-length H4:
%   0   → peptide N-terminus
%   2   → K5 within this peptide
%   5   → K8 within this peptide
%   9   → K12 within this peptide
%   13  → K16 within this peptide
% For example, unmod has pr at 0,2,5,9,13 (derivatized, no acetylation).
His.mod_type = {'0,pr;2,pr;5,pr;9,pr;13,pr;';
    '0,pr;2,ac;5,pr;9,pr;13,pr;';
    '0,pr;2,pr;5,ac;9,pr;13,pr;';
    '0,pr;2,pr;5,pr;9,ac;13,pr;';
    '0,pr;2,pr;5,pr;9,pr;13,ac;';
    '0,pr;2,ac;5,ac;9,pr;13,pr;';
    '0,pr;2,ac;5,pr;9,ac;13,pr;';
    '0,pr;2,ac;5,pr;9,pr;13,ac;';
    '0,pr;2,pr;5,ac;9,ac;13,pr;';
    '0,pr;2,pr;5,ac;9,pr;13,ac;';
    '0,pr;2,pr;5,pr;9,ac;13,ac;';
    '0,pr;2,pr;5,ac;9,ac;13,ac;';
    '0,pr;2,ac;5,pr;9,ac;13,ac;';
    '0,pr;2,ac;5,ac;9,pr;13,ac;';
    '0,pr;2,ac;5,ac;9,ac;13,pr;';
    '0,pr;2,ac;5,ac;9,ac;13,ac;'};

% Charge grid (1+,2+,3+,4+) replicated for each row.
% Keep 1+ because this short, Lys-rich peptide can show 1+ signals.
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);

% Theoretical m/z per row x charge. We compute it to keep the panel portable
% across forks (instead of hard-coding). calculate_pepmz reads mod_type.
His.pep_mz = calculate_pepmz(His);

% Seed RTs (minutes) used to build relocation windows.
% The unmodified derivatized peptide is the anchor (row 1).
His.rt_ref = [30.08
    28.80
    28.80
    28.79
    29.11
    27.78
    27.72
    27.94
    27.78
    27.94
    27.88
    26.72
    26.65
    26.76
    26.55
    25.35];

% Display flags for plotting. All ones here (we want all shown).
His.display = ones(length(His.mod_type),1);

% Bookkeeping for downstream writers/plotters
His.outpath = cur_outpath;
His.outfile = out_filename;

% For visual consistency across panels, re-order charge columns so that
% 2+ appears first, followed by the remaining charges (order preserved).
% This does NOT change extraction logic — only the column order.
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end
    tune = 1:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end

end % init_histone


% =========================================================================
% NESTED: calculate_layout
%   Main extraction pipeline:
%     - Anchor unmodified RT via check_ref and get_histone0
%     - Calibrate drift and shift seed RTs of modified forms
%     - Relocate windows (MS1-only or MS2-assisted)
%     - Extract mono/di/tri/tetra acetylation blocks
% =========================================================================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,Mods,His,special)
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);

% Pre-allocate outputs:
%   - pep_rts/pep_intens are row-wise (modification rows) with per-charge columns
%   - mono_isointens stores the monoisotopic XIC trace for each row
pep_rts        = zeros([npep,ncharge]);
pep_intens     = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% ------------- 1) Anchor the unmodified derivatized peptide -------------
His.rt_unmod_orig = His.rt_ref(1); % keep the original seed to compute drift later

if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only relocation: check_ref may refine the unmod RT around raw_path
        [His.rt_ref(1),special.ndebug] = ...
            check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                      His.rt_ref(1),special.ndebug);
    else
        % DA relocation path: still call check_ref to propose a center
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                                  His.rt_ref(1),special.ndebug);
        % If the check_ref anchor equals the original seed, widen the probe;
        % else search narrowly around the refined RT (±5 min).
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
        else
            delta = 5;
            t1 = His.rt_ref(1)-delta;
            t2 = His.rt_ref(1)+delta;
        end
        % Probe unmodified with get_rts2 (MS2-assisted) to sharpen the anchor
        hno = 1; % unmod is row 1
        [rts1,top1_rt1] = ...
            get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% Finalize unmodified extraction via get_histone0 (returns per-charge RT/area)
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% If unmodified was found (cur_rts(1)>0), calibrate drift globally:
%   - set anchor RT to observed RT
%   - compute delta vs the original seed
%   - shift all remaining seed RTs by that delta
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    % Record results for the unmodified row (row 1)
    pep_rts(hno,1:ncharge)            = cur_rts;
    pep_intens(hno,1:ncharge)         = cur_intens;
    mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
end

% ------------- 2) Relocate windows for modified forms -------------
% Choose relocation strategy based on nDAmode (MS1-only vs MS2-assisted).
if 1==special.ndebug
    % Diagnostic relocation (minimal logic; useful when checking inputs)
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end
end

% ------------- 3) Extract blocks of acetylation states -------------

% MONO-ACETYLATIONS (K5ac/K8ac/K12ac/K16ac)
% We use a combinatorial helper (get_histone4) that understands the
% mutually related rows and can disambiguate partially coeluting isomers.
hno = 2; % block starts at row 2
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone4(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                 ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+3,1:ncharge)            = cur_rts(1:4,:);
    pep_intens(hno:hno+3,1:ncharge)         = cur_intens(1:4,:);
    mono_isointens(1:num_MS1,hno:hno+3)     = cur_mono_isointens(:,1:4);
end

% DI-ACETYLATIONS (6 combinations): rows 6..11
% Similar idea: a helper that extracts a 6-row combinatorial block.
hno = 6; % block starts at row 6
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone6(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                 ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+5,1:ncharge)            = cur_rts(1:6,:);
    pep_intens(hno:hno+5,1:ncharge)         = cur_intens(1:6,:);
    mono_isointens(1:num_MS1,hno:hno+5)     = cur_mono_isointens(:,1:6);
end

% TRI-ACETYLATIONS (4 combinations): rows 12..15
% Again leverage the 4-row helper; the order is adjusted later by change_order.
hno = 12; % block starts at row 12
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone4(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                 ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+3,1:ncharge)            = cur_rts(1:4,:);
    pep_intens(hno:hno+3,1:ncharge)         = cur_intens(1:4,:);
    mono_isointens(1:num_MS1,hno:hno+3)     = cur_mono_isointens(:,1:4);
end

% TETRA-ACETYLATION (all four Lys acetylated): row 16
% This is a single row; use the simple extractor.
hno = 16;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)              = cur_rts;
    pep_intens(hno,1:ncharge)           = cur_intens;
    mono_isointens(1:num_MS1,hno)       = cur_mono_isointens;
end

end % calculate_layout


% =========================================================================
% NESTED: relocate (MS1-only relocation)
%   Build RT search windows for target rows relative to the unmodified RT.
%   Heuristics reflect that higher acetylation generally elutes earlier.
% =========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1;   % small gap to avoid capturing anchor tailing
nsplit = 1;    % pass-through to get_rts (your core may split/merge windows)

% --- K5ac (row 2): search BEFORE unmodified anchor ---
% Window: [unmod - 11, unmod - delta]
hno = 2;
t1 = His.rt_ref(1)-11;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% Remember original seed; we will use it to propagate changes to K8/K12/K16
old_t = His.rt_ref(hno);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;          % not found under current window
else
    His.rt_ref(hno) = top1_rt2;   % place ref at the main RT peak
end

% --- K8ac/K12ac/K16ac (rows 3..5): dependent on K5ac outcome ---
% If K5ac failed, set them to zero (do not hunt widely).
% If K5ac ref moved away from its old seed, shift the trio together so the
% relative order is preserved visually and for subsequent extraction.
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:5) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:5) = His.rt_ref(hno-1);
end

% --- K5acK8ac (row 6): earlier than both K5ac and unmod ---
% Window: [unmod - 17, min(K5ac, unmod) - delta]
hno = 6;
t1 = His.rt_ref(1)-17;
if 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end
[rts6,top1_rt6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end

% --- Remaining di-acets (rows 7..11): cascade from K5acK8ac outcome ---
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:11) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:11) = His.rt_ref(hno-1);
end

% --- Tri-acetyl K8acK12acK16ac (row 12): even earlier ---
% Window: [unmod - 23, min(K5acK8ac, K5ac, unmod) - delta]
hno = 12;
t1 = His.rt_ref(1)-23;
if 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end
[rts12,top1_rt12] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

old_t = His.rt_ref(hno);
if 1==isempty(rts12)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt12;
end

% --- Remaining tri-acets (rows 13..15): cascade from row 12 outcome ---
hno = 13;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:15) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:15) = His.rt_ref(hno-1);
end

% --- Tetra-acetyl (row 16): earliest in this panel ---
% Window: [unmod - 28, min(tri-acet, di-acet, mono-acet, unmod) - delta]
hno = 16;
t1 = His.rt_ref(1)-28;
if 0<His.rt_ref(12)
    t2 = His.rt_ref(12)-delta;
elseif 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end
[rts16,top1_rt16] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts16)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt16;
end

end % relocate


% =========================================================================
% NESTED: relocate2 (MS2-assisted relocation; same logic, but using get_rts2)
% =========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K5ac
hno = 2;
t1 = His.rt_ref(1)-11;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% K8ac/K12ac/K16ac (propagate from K5ac result)
hno = 3;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:5) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:5) = His.rt_ref(hno-1);
end

% K5acK8ac
hno = 6;
t1 = His.rt_ref(1)-17;
if 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts6)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt6;
end

% Remaining di-acets
hno = 7;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:11) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:11) = His.rt_ref(hno-1);
end

% K8acK12acK16ac
hno = 12;
t1 = His.rt_ref(1)-23;
if 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end
[rts12,top1_rt12] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

old_t = His.rt_ref(hno);
if 1==isempty(rts12)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt12;
end

% Remaining tri-acets
hno = 13;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno:15) = 0;
elseif old_t~= His.rt_ref(hno-1);
    His.rt_ref(hno:15) = His.rt_ref(hno-1);
end

% Tetra-acetyl
hno = 16;
t1 = His.rt_ref(1)-28;
if 0<His.rt_ref(12)
    t2 = His.rt_ref(12)-delta;
elseif 0<His.rt_ref(6)
    t2 = His.rt_ref(6)-delta;
elseif 0<His.rt_ref(2)
    t2 = His.rt_ref(2)-delta;
else
    t2 = His.rt_ref(1)-delta;
end
[rts16,top1_rt16] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts16)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt16;
end

end % relocate2


% =========================================================================
% NESTED: change_order
%   Post-processing to reorder the 4 tri-acetyl rows (12..15 → 15..14..13..12)
%   for a clearer visual grouping in plots/reports. Quantitative values are
%   copied; this only changes display order.
% =========================================================================
function [His,pep_rts,pep_intens,mono_isointens] = ...
    change_order(His,pep_rts,pep_intens,mono_isointens)
%% 12,13,14,15->15,14,13,12

% Copy current values for rows 12..15
b1 = His.mod_short{15};
b2 = His.mod_short{14};
b3 = His.mod_short{13};
b4 = His.mod_short{12};

c1 = His.mod_type{15};
c2 = His.mod_type{14};
c3 = His.mod_type{13};
c4 = His.mod_type{12};

d1 = His.rt_ref(15);
d2 = His.rt_ref(14);
d3 = His.rt_ref(13);
d4 = His.rt_ref(12);

e1 = pep_rts(15,:);
e2 = pep_rts(14,:);
e3 = pep_rts(13,:);
e4 = pep_rts(12,:);

f1 = pep_intens(15,:);
f2 = pep_intens(14,:);
f3 = pep_intens(13,:);
f4 = pep_intens(12,:);

g1 = mono_isointens(:,15);
g2 = mono_isointens(:,14);
g3 = mono_isointens(:,13);
g4 = mono_isointens(:,12);

% Paste back in reversed order (15→12, 14→13, 13→14, 12→15)
His.mod_short{12} = b1;
His.mod_short{13} = b2;
His.mod_short{14} = b3;
His.mod_short{15} = b4;

His.mod_type{12} = c1;
His.mod_type{13} = c2;
His.mod_type{14} = c3;
His.mod_type{15} = c4;

His.rt_ref(12) = d1;
His.rt_ref(13) = d2;
His.rt_ref(14) = d3;
His.rt_ref(15) = d4;

pep_rts(12,:) = e1;
pep_rts(13,:) = e2;
pep_rts(14,:) = e3;
pep_rts(15,:) = e4;

pep_intens(12,:) = f1;
pep_intens(13,:) = f2;
pep_intens(14,:) = f3;
pep_intens(15,:) = f4;

mono_isointens(:,12) = g1;
mono_isointens(:,13) = g2;
mono_isointens(:,14) = g3;
mono_isointens(:,15) = g4;

end % change_order
