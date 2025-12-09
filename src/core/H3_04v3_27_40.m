function H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H3_04v3_27_40
% ------------------------------------------------------------------------------
% Purpose:
%   Layout builder for the tryptic H3 peptide 27–40 (sequence set below),
%   enumerating a panel of K27/K36 PTM states (mono/di/tri-methyl, acetyl, and
%   selected combinations), extracting MS1 profiles (and optionally MS2 checks),
%   calibrating RTs to the current run, outputting .mat results, and drawing
%   figures.
%
% Inputs:
%   MS1_index, MS1_peaks : MS1 scan table and corresponding peak container
%   MS2_index, MS2_peaks : MS2 scan table and corresponding peak container
%   ptol                 : mass tolerance selector (100 is special-cased to 10 ppm)
%   cur_outpath          : output directory for files and figures of this layout
%   special              : struct of run-time switches (fields used here:
%                          .nDAmode, .ndebug, .nhmass and others passed through)
%
% Outputs (files):
%   cur_outpath/out_filename.mat       : serialized 'His' + intensities/RTs
%   figures and ancillary outputs via draw_layout()/output_histone()
%
% Notes:
%   - This function delegates most heavy lifting to helper routines:
%       * init_histone(): defines the peptide, PTM list, charges, RT refs
%       * calculate_layout(): orchestrates extraction per PTM and calibration
%       * relocate()/relocate2(): optional re-alignment of RT references
%   - The peptide sequence used here is KSAPTTGGVKKPHR (A. thaliana H3.3 diag.)
% ------------------------------------------------------------------------------

%{
% Optional: force tighter tolerance for specific chemistries (kept as author note)
if ptol<100
    ptol = 10; % for K27me2K36me2 and K27me3K36me2
end
%}

% -----------------------------
% 1) Skip if output already exists
% -----------------------------
out_filename = 'H3_04v3_27_40';                  % layout name (used in filenames)
out_file0    = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    % Idempotency: do nothing if results for this layout are already present
    return;
end

% -----------------------------
% 2) Build histone/peptide definition and prior RTs
% -----------------------------
His = init_histone(cur_outpath,out_filename);

% -----------------------------
% 3) Compute: extract profiles, calibrate RTs, populate arrays
% -----------------------------
unitdiff = 1.0032;                 % isotopic mass spacing (approx. 1 Da)
Mods     = GetMods();              % PTM dictionary for MS2 checks (external)
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

% -----------------------------
% 4) Save per-peptide results
% -----------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% -----------------------------
% 5) Draw figures for report
% -----------------------------
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2);  % RT axis for MS1 scans
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
            mono_isointens,MS2_index,MS2_peaks,special);

% -----------------------------
% 6) (Optional) Extract PSMs in "DA mode"
% -----------------------------
if 1==special.nDAmode
    % DA mode: produce PSM-level evidence summary for the panel
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
           mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

% ==========================================================================
function His = init_histone(cur_outpath,out_filename)
%%
% Initialize peptide, PTM panel, charge states, theoretical m/z, RT priors,
% visibility flags, and normalize "main charge" column ordering if needed.

His.pep_seq = 'KSAPTTGGVKKPHR';  % H3(27–40) diagnostic peptide (A. thaliana H3.3)

% Short labels for each PTM state in this panel (n = 17 states)
His.mod_short = {'unmod';
    'K36me1';
    'K27me1';
    'K27me2';
    'K36me2';
    'K27me3';
    'K36me3';
    'K27me2K36me1';
    'K27me1K36me2';
    'K27me1K36me1';
    'K27me3K36me1';
    'K27me1K36me3';
    'K27me2K36me2';
    'K27me3K36me2';
    'K27me2K36me3';
    'K27me3K36me3';
    'K27ac'};

% PTM encoding per state using internal position IDs:
% Format example: '0,pr;1,me1;10,pr;11,pr;' → (positionID, PTM) pairs.
% 'pr' stands for "present/unmodified" at that position ID.
% NOTE [Inference]: position IDs {0,1,10,11} are the framework's internal map
% for residues within KSAPTTGGVKKPHR that correspond to K27/K36 loci.
His.mod_type = {'0,pr;1,pr;10,pr;11,pr;';
    '0,pr;1,pr;10,me1;11,pr;';
    '0,pr;1,me1;10,pr;11,pr;';
    '0,pr;1,me2;10,pr;11,pr;';
    '0,pr;1,pr;10,me2;11,pr;';
    '0,pr;1,me3;10,pr;11,pr;';
    '0,pr;1,pr;10,me3;11,pr;';
    '0,pr;1,me2;10,me1;11,pr;';
    '0,pr;1,me1;10,me2;11,pr;';
    '0,pr;1,me1;10,me1;11,pr;';
    '0,pr;1,me3;10,me1;11,pr;';
    '0,pr;1,me1;10,me3;11,pr;';
    '0,pr;1,me2;10,me2;11,pr;';
    '0,pr;1,me3;10,me2;11,pr;';
    '0,pr;1,me2;10,me3;11,pr;';
    '0,pr;1,me3;10,me3;11,pr;';
    '0,pr;1,ac;10,pr;11,pr;'};

% Charge states per PTM row → use z = 2,3,4 for all
His.pep_ch = repmat([2 3 4],length(His.mod_type),1);

% Theoretical m/z for each (PTM state × charge) → computed elsewhere to
% guarantee consistency with 'mod_type' and 'pep_seq'.
His.pep_mz = calculate_pepmz(His);

% Prior RTs (minutes) for each PTM state (empirical, method-specific).
% These act as initial anchors before calibration.
His.rt_ref = [26.08
    27.81
    27.94
    21.40
    22.68
    21.40
    22.66
    22.64
    24.44
    29.15
    22.64
    24.44
    18.29
    18.20
    18.34
    18.20
    24.99];

% Visibility mask for plotting/reporting (1=show, 0=hide)
His.display = ones(length(His.mod_type),1);
His.display([15 16]) = 0;  % hide K27me2K36me3 and K27me3K36me3 in figures by default

% Output routing
His.outpath  = cur_outpath;
His.outfile  = out_filename;

% Normalize "main charge" column across PTM rows:
% This ensures the preferred charge (column 2 of row 1) is the first column
% across a selected subset of PTM rows, keeping arrays consistent.
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz); %#ok
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];  % reorder columns so main_ch is first
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end
    tune = [4 5 6 7 8 9 11 12 13 14 15 16]; % rows to apply column reordering
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end

% ==========================================================================
function [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special)
%%
% Master orchestrator for the per-PTM extraction:
%   - Calibrate RT using the unmodified form as anchor
%   - Optionally relocate (MS1-only or MS1+MS2)
%   - Choose specialized extractors for close/co-eluting species
%   - Fill per-PTM RT/intensity matrices and monoisotopic MS1 traces

[npep,ncharge] = size(His.pep_mz);
num_MS1        = size(MS1_index,1);
pep_rts        = zeros([npep,ncharge]);
pep_intens     = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% --- 1) Unmodified anchor: extract and calibrate run RTs
His.rt_unmod_orig = His.rt_ref(1);       % keep original prior for reporting
nhmass   = special.nhmass;               % heavy-N toggle (used later)
hno      = 1;                            % row index: 'unmod'
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% If we found a valid RT for unmod (first charge column), shift all RT priors
% by the delta observed in this run. This "calibrates" the panel to the run.
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta; % global shift
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% --- 2) Optional RT relocation using empirical patterns
%   - relocateD: developer/debug routine (if special.ndebug==1)
%   - relocate : MS1-only re-alignment (default when nDAmode~=2)
%   - relocate2: MS1+MS2 re-alignment (if nDAmode==2)
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end
end

% --- 3) Close pair K36me1 / K27me1 → pick extractor depending on separation
% If the two RTs differ by >0.4 min, treat independently (MS1-only extractor).
% Otherwise, use a coupled extractor with MS2 checks to disentangle them.
if His.rt_ref(3)-His.rt_ref(2)>0.4
    for hno=2:3
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone10( ...
            MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge)            = cur_rts;
            pep_intens(hno,1:ncharge)         = cur_intens;
            mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
        end
    end
else
    hno = 2; % starting at K36me1: get_histone2 returns a 2-row block (K36me1/K27me1)
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
        MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge)            = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge)         = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1)     = cur_mono_isointens(:,1:2);
    end
end

% --- 4) K27me2, K27me3 (MS1-only extractor 1) and K36me2 (extractor 11)
for hno=[4 6] % K27me2 and K27me3
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1( ...
        MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)            = cur_rts;
        pep_intens(hno,1:ncharge)         = cur_intens;
        mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
    end
end
hno = 5; % K36me2 via extractor 11 (tighter RT logic)
[cur_rts,cur_intens,cur_mono_isointens] = get_histone11( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)            = cur_rts;
    pep_intens(hno,1:ncharge)         = cur_intens;
    mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
end

% --- 5) K36me3 + K27me2K36me1 → coupled extractor with MS2 verification
hno = 7; % starting at K36me3, returns 2-row block (7 and 8)
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge)        = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge)     = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end

% --- 6) Sequence of combinations (9..13) via MS1-only extractor 1
for hno=9:13
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1( ...
        MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)            = cur_rts;
        pep_intens(hno,1:ncharge)         = cur_intens;
        mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
    end
end

% --- 7) K27me3K36me2 and (optionally) K27me2K36me3
% Default path: K27me3K36me2 via MS1-only extractor 1
hno = 14;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)            = cur_rts;
    pep_intens(hno,1:ncharge)         = cur_intens;
    mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
end

% Alternative (commented by author): jointly with MS2 using get_histone2
%{
hno = 14;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge)        = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge)     = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end
%}

% --- 8) (Optional) K27me3K36me3 via MS1-only extractor 1 (author kept commented)
%{
hno = 16;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)            = cur_rts;
    pep_intens(hno,1:ncharge)         = cur_intens;
    mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
end
%}

% --- 9) K27ac via MS1-only extractor 1
hno = 17;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)            = cur_rts;
    pep_intens(hno,1:ncharge)         = cur_intens;
    mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
end

% ==========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
% Lightweight relocation based on an external RT reference table for the
% sister peptide KSAPATGGVKKPHR (H3.1) stored in 'H3_04_27_40.xls'.
% It estimates a global shift 'nshift' (depends on run length) and, for a
% selected set of PTM rows, re-assigns RTs by searching a narrow window
% [ref_rts(hno)-nshift±delta] with get_rts().

delta  = 0.5;
nsplit = 1;

% Estimate run length (end of MS1)
num_MS1 = size(MS1_index,1);
end_rt  = MS1_index(num_MS1,2);
if end_rt>90
    nshift = 1.6;  % larger global shift for very long gradients
else
    nshift = 0.3;  % smaller shift otherwise
end

% Read external RT list from file (if present)
ref_rts = get_KSAPATGGVKKPHR_rt(His);
if 1==isempty(ref_rts)
    return; % nothing to do if external refs are unavailable
end

% Reassign RTs for a subset of PTM states by local search around ref_rts-nshift
for hno = [2 3 4 5 6 7 8 9 10 11 12 13 14 17]
    t1 = ref_rts(hno)-nshift-delta;
    t2 = ref_rts(hno)-nshift+delta;
     % ===== DEBUG: inspección de ventana previa a get_rts =====
        % — Guard: ref_rts inválido —
    if ref_rts(hno) <= 0 || isnan(ref_rts(hno))
        His.rt_ref(hno) = 0;
        fprintf('[relocate INFO] hno=%d: ref_rts<=0 → skip\n', hno);
        continue;
    end

    % — Clamp al rango del run —
    minRT = MS1_index(1,2); maxRT = MS1_index(end,2);
    t1c = max(t1, minRT);
    t2c = min(t2, maxRT);
    if t1c >= t2c
        His.rt_ref(hno) = 0;
        fprintf('[relocate INFO] hno=%d: ventana degenerada → skip\n', hno);
        continue;
    end

    fprintf('[relocate CHECK] calling get_rts with [%.3f, %.3f]\n', t1c, t2c);

    minRT = MS1_index(1,2); 
    maxRT = MS1_index(end,2);
    n_ge_t1 = sum(MS1_index(:,2) >= t1);
    n_le_t2 = sum(MS1_index(:,2) <= t2);
    n_in    = sum(MS1_index(:,2) >= t1 & MS1_index(:,2) <= t2);
    fprintf(['[relocate DEBUG] hno=%d | ref=%.3f | nshift=%.3f | delta=%.3f | ', ...
             't1=%.3f t2=%.3f | RUN=[%.3f,%.3f] | n>=t1=%d n<=t2=%d n(in)=%d\n'], ...
            hno, ref_rts(hno), nshift, delta, t1, t2, minRT, maxRT, n_ge_t1, n_le_t2, n_in);
    % ===== FIN DEBUG =====
    %[rts1,top1_rt1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
    [rts1,top1_rt1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1c,t2c);

    if 1==isempty(rts1)
        His.rt_ref(hno) = 0;       % mark as not found
    else
        His.rt_ref(hno) = top1_rt1; % use the best local peak
    end
end

% ==========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
% Same idea as relocate(), but uses get_rts2() (MS1+MS2 confirmation),
% which is typically more robust in complex backgrounds.

delta  = 0.5;
nsplit = 1;

num_MS1 = size(MS1_index,1);
end_rt  = MS1_index(num_MS1,2);
if end_rt>90
    nshift = 1.6;
else
    nshift = 0.3;
end

ref_rts = get_KSAPATGGVKKPHR_rt(His);
if 1==isempty(ref_rts)
    return;
end

for hno = [2 3 4 5 6 7 8 9 10 11 12 13 14 17]
    t1 = ref_rts(hno)-nshift-delta;
    t2 = ref_rts(hno)-nshift+delta;
     % ===== DEBUG: inspección de ventana previa a get_rts =====
    minRT = MS1_index(1,2); 
    maxRT = MS1_index(end,2);
    n_ge_t1 = sum(MS1_index(:,2) >= t1);
    n_le_t2 = sum(MS1_index(:,2) <= t2);
    n_in    = sum(MS1_index(:,2) >= t1 & MS1_index(:,2) <= t2);
    fprintf(['[relocate DEBUG] hno=%d | ref=%.3f | nshift=%.3f | delta=%.3f | ', ...
             't1=%.3f t2=%.3f | RUN=[%.3f,%.3f] | n>=t1=%d n<=t2=%d n(in)=%d\n'], ...
            hno, ref_rts(hno), nshift, delta, t1, t2, minRT, maxRT, n_ge_t1, n_le_t2, n_in);
    % ===== FIN DEBUG =====
    [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1c,t2c,nhmass);
    if 1==isempty(rts1)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt1;
    end
end

% ==========================================================================
function ref_rts = get_KSAPATGGVKKPHR_rt(His)
%%
% Read an external RT table produced previously for the H3(27–40) sister
% peptide KSAPATGGVKKPHR (H3.1 isoform) from an .xls-like text file.
% Expected format:
%   ... lines ...
%   [rt]
%   peptide
%   <tab-separated rows where first two fields include the RT>
% Lines are scanned until '[rt]' then subsequent lines are parsed.

out_file0 = fullfile(His.outpath,'H3_04_27_40.xls');
if 0~=exist(out_file0,'file')
    fp = fopen(out_file0,'r');
    str = fgetl(fp);
    while 0==feof(fp) && 0==strcmp(str,'[rt]')
        str = fgetl(fp);
    end
    str = fgetl(fp); % skip "peptide" header
    no = 0;
    ref_rts = [];
    while 0==feof(fp)
        str = fgetl(fp);
        p   = strfind(str,'	');          % tab-separated
        c_rt = str2num( str(p(1)+1:p(2)-1) ); %#ok<ST2NM> % parse RT number
        no = no + 1;
        ref_rts(no) = c_rt; %#ok<AGROW>   % append
    end
    fclose(fp);
else
    ref_rts = []; % missing external RT table → relocation will be skipped
end

% --------------------------------------------------------------------------
% Below: the author kept a fuller relocate/relocate2 version (with pattern-
% based co-elution logic and find_triple()) as commented reference. The live
% relocate()/relocate2() above are the "short" variant that uses an external
% RT table as anchor and local search windows. Keeping this here documents
% the original, more elaborate approach for future extension.
%--------------------------------------------------------------------------
% { ... long commented relocate / relocate2 with find_triple ... }
%--------------------------------------------------------------------------
