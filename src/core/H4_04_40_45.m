function H4_04_40_45(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ========================================================================
% H4_04_40_45 — Targeted panel for histone H4 peptide "RGGVKR" (aa 40–45)
%               States: unmodified (derivatized baseline) and K44ac.
% ------------------------------------------------------------------------
% Rationale:
% - Short H4 peptide containing the biologically relevant Lys at H4K44.
% - In derivatization workflows, "unmod" means N-term and internal Lys are
%   propionylated (baseline for chromatographic behavior).
% - K44ac is expected to elute EARLIER than its propionylated counterpart.
%
% Inputs:
%   MS1_index, MS1_peaks : centroided MS1 index (incl. RT in minutes) + peaks
%   MS2_index, MS2_peaks : centroided MS2 (for DA relocation / PSM export)
%   ptol                 : precursor mass tolerance in ppm
%   cur_outpath          : output folder
%   special              : struct with fields:
%                          - ndebug   : 1 -> diagnostic relocateD; else normal
%                          - nDAmode  : 0/other = MS1-only; 2 = MS2-assisted;
%                                       1 = export PSM after quantification
%                          - nhmass   : helper neutral mass (if used by get_rts2)
%                          - raw_path : raw path used by check_ref
%
% Outputs (to disk):
%   <cur_outpath>/H4_04_40_45.mat with His, pep_rts, pep_intens, mono_isointens
%   Layout figures via draw_layout; optional PSM via GetPSM if nDAmode==1.
% ========================================================================

% ------------------------------
% 1) Idempotent check: skip if .mat result exists
% ------------------------------
out_filename = 'H4_04_40_45';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return; % Respect original short-circuit behavior
end

% ------------------------------
% 2) Initialize panel definition (sequence, PTMs, charges, m/z, RT seeds)
% ------------------------------
His = init_histone();

% ------------------------------
% 3) Main calculation (anchor, relocate, extract)
% ------------------------------
unitdiff = 1.0032; % ~C13 spacing used by monoisotopic helpers
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);

% ------------------------------
% 4) Persist results
% ------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% ------------------------------
% 5) Draw XIC / layout for visual QC
% ------------------------------
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2); % RT vector (minutes)
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------
% 6) Optional: export PSM if requested (nDAmode==1)
% ------------------------------
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % H4_04_40_45


% ========================================================================
% init_histone — define peptide, modification rows, charges, m/z, RT seeds
% ========================================================================
function His = init_histone()
%%
% Peptide sequence for H4 (aa 40–45)
His.pep_seq = 'RGGVKR';

% Two PTM rows:
% - 'unmod' : derivatized baseline (N-term + Lys propionylated)
% - 'K44ac' : acetylated Lys (index 5 in this peptide), N-term still propionylated
His.mod_short = {'unmod';
    'K44ac'};

% PTM encoding (peptide-relative indices):
%   0,pr -> N-terminal propionylation
%   5,pr -> Lys propionylation (for "unmod")
%   5,ac -> Lys acetylation (for "K44ac")
His.mod_type = {'0,pr;5,pr;';
    '0,pr;5,ac;'};

% Charge grid considered by the panel (1+ and 2+). The utility functions
% will reorganize so that the "main" charge (2+) sits on the first column.
His.pep_ch = repmat([1 2],length(His.mod_type),1);

%{
% Historical fixed m/z (for debugging); superseded by calculate_pepmz(His):
His.pep_mz = [784.4788  392.7430
              770.4631  385.7352];
%}

% Compute theoretical m/z dynamically from sequence + PTM specification
His.pep_mz = calculate_pepmz(His);

% Seed retention times (minutes); K44ac typically elutes earlier than unmod
His.rt_ref = [20.47
    18.70];

% Display control (1=visible in composite plots)
His.display = ones(length(His.mod_type),1);

% Cosmetic reorganization: make 2+ the first column for alignment,
% without altering numeric values.
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


% ========================================================================
% calculate_layout — anchor unmod, drift-correct, relocate, extract K44ac
% ========================================================================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special)
%%
[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);

% Preallocate outputs
pep_rts        = zeros([npep,ncharge]);
pep_intens     = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% ---- 1) Anchor unmodified (derivatized baseline) using check_ref
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only anchoring
        [His.rt_ref(1),special.ndebug] = ...
            check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                      His.rt_ref(1),special.ndebug);
    else
        % MS2-assisted anchoring: refine around check_ref (±5 min) or scan all
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                                  His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0; t2 = MS1_index(num_MS1,2);
        else
            delta = 5; t1 = His.rt_ref(1)-delta; t2 = His.rt_ref(1)+delta;
        end
        hno = 1; % unmod
        [rts1,top1_rt1] = ...
            get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% ---- 2) Extract unmodified to compute delta and correct global drift
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

if cur_rts(1)>0
    % Update anchor and propagate delta to remaining seed RTs
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    % Persist extracted metrics for the unmodified row
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% ---- 3) Relocate K44ac search window (MS1-only or MS2-assisted)
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His); % minimal diag path
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His); % MS1-only
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                        ptol,unitdiff,His,nhmass); % DA path
    end
end

% ---- 4) Extract K44ac intensities/RT per charge (row 2)
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

end % calculate_layout


% ========================================================================
% relocate — MS1-only relocation for K44ac (earlier than unmod)
% ========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
delta  = 0.1; % small exclusion to avoid overlapping the anchor peak apex
nsplit = 1;   % split control passed to get_rts (panel convention)

% K44ac search window: earlier than unmodified (14 min before until just before)
hno = 2;
t1 = His.rt_ref(1)-14;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0; % mark as not found
else
    His.rt_ref(hno) = top1_rt2; % adopt top-peak RT as seed
end

end % relocate


% ========================================================================
% relocate2 — MS2-assisted relocation for K44ac (same window, get_rts2)
% ========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
delta  = 0.1;
nsplit = 1;

hno = 2;
t1 = His.rt_ref(1)-14;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

end % relocate2
