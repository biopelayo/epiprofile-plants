function H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H3_08_117_128
% -------------------------------------------------------------------------
% Purpose:
%   Quantify H3 peptide VTIMPKDIQLAR (117–128) and its common PTMs:
%     - unmod  : derivatized termini (EpiProfile 'pr' convention)
%     - M120ox : methionine oxidation (typically elutes earlier)
%     - K122ac : lysine acetylation (elutes slightly before unmod)
%
% Strategy:
%   1) Refine unmodified anchor RT (check_ref; DA probing optional).
%   2) Extract unmodified, then calibrate global drift based on observed RT.
%   3) Relocate modified rows (MS1-only vs MS2-assisted (DA) windows).
%   4) Extract modified rows (get_histone1) and build XIC matrices.
%   5) Save .mat, draw layout, optionally export PSM tables.
%
% Inputs:
%   MS1_index, MS1_peaks : centroided MS1 (core-internal format).
%   MS2_index, MS2_peaks : centroided MS2 (required for DA relocation).
%   ptol                 : mass tolerance (ppm or instrument-specific units).
%   cur_outpath          : output directory for .mat and figures.
%   special              : struct with flags/hints:
%       .ndebug   -> 1: use relocateD (diagnostic mode)
%       .nDAmode  -> relocation/PSM semantics (historical):
%                     2 => DA relocation (relocate2 + get_rts2)
%                     1 => export PSMs after quantification
%       .nhmass   -> neutral-loss/precursor mass hint for get_rts2
%       .raw_path -> vendor raw path for check_ref
%
% Outputs (to disk):
%   <cur_outpath>/H3_08_117_128.mat  (pep_rts, pep_intens, mono_isointens)
%   Figures via draw_layout; PSM tables if special.nDAmode == 1
% -------------------------------------------------------------------------

% ------------------------------- check -----------------------------------
out_filename = 'H3_08_117_128';
fprintf(1,'%s..',out_filename);  % simple progress print, preserved
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    % Idempotency: skip if results already exist.
    return;
end

% -------------------------------- init -----------------------------------
His = init_histone();

% ------------------------------ calculate --------------------------------
unitdiff = 1.0032;  % 13C spacing for monoisotopic handling
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special);

% -------------------------------- output ---------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% --------------------------------- draw ----------------------------------
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2); % RT axis from MS1 index table
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
            mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------- Get PSM ---------------------------------
if 1==special.nDAmode
    % Historical: nDAmode==1 triggers PSM export post-quantification.
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
           mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % main

% =========================================================================
function His = init_histone()
%%
% Build the "His" struct: sequence, PTM table, charge grid, m/z table,
% seed RTs, and display flags. Core-dependent PTM indexing applies.

His.pep_seq = 'VTIMPKDIQLAR';

His.mod_short = {'unmod';
                 'M120ox';
                 'K122ac'};

% Compact PTM encodings (do NOT change unless your core uses different
% residue indices for this panel):
% - '0,pr' => peptide termini derivatization; '6,pr' marks terminal handling.
% - '4,ox' => methionine oxidation (core-dependent index for Met).
% - '6,ac' => acetylation on the Lys residue (core-dependent index).
His.mod_type = {'0,pr;6,pr;';
                '0,pr;4,ox;6,pr;';
                '0,pr;6,ac;'};

% Allow charge states 1+..4+, replicated per PTM row
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);

% Legacy hard-coded masses retained as comment (for cross-checking only):
%{
His.pep_mz = [1496.8505  748.9289  499.6217  374.9681;
              1512.8454  756.9264  504.9533  378.9668;
              1482.8349  741.9211  494.9498  371.4642];
%}

% Compute portable m/z values from composition and mod_type
His.pep_mz = calculate_pepmz(His);

% Seed RTs (minutes): unmod, M120ox (earlier), K122ac (slightly earlier)
His.rt_ref = [48.69;
              46.37;
              47.68];

% Display flags: show all in figures, but hide M120ox by default if desired
His.display = ones(length(His.mod_type),1);
His.display(2) = 0;

% ----------------------- charge column reordering ------------------------
% Panel-wide convention: put 2+ as first column for easier downstream merge
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
function [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% High-level orchestration: anchor refinement, unmod extraction + drift
% calibration, relocation (MS1-only vs DA), final extraction of modifieds.

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% ---------------------------- anchor: unmod -------------------------------
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        [His.rt_ref(1),special.ndebug] = check_ref( ...
            special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
    else
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref( ...
            special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        % If anchor remained unchanged, probe a broad symmetric window in DA
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
        else
            delta = 5;
            t1 = His.rt_ref(1)-delta;
            t2 = His.rt_ref(1)+delta;
        end
        hno = 1; % unmod
        [~,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                                 ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% Extract the unmodified row and calibrate drift using the observed RT
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% ------------------------------ relocation -------------------------------
if 1==special.ndebug
    % Diagnostic relocation (no window logic change)
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special.nhmass);
    end
end

% ----------------------------- final extraction --------------------------
% M120ox and K122ac (rows 2..3)
for hno=2:3
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1( ...
        MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

end % calculate_layout

% =========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
% MS1-only relocation windows anchored to the observed unmod RT:
%   - M120ox : [unmod - 10, unmod - 0.1]
%   - K122ac : [max(unmod - 8, M120ox + 0.1), unmod - 0.1]

delta  = 0.1;
nsplit = 1;

% -------------------------------- M120ox ---------------------------------
hno = 2;
t1 = His.rt_ref(1)-10;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% -------------------------------- K122ac ---------------------------------
hno = 3;
t1 = max([His.rt_ref(1)-8, His.rt_ref(2)+delta]); % avoid overlap with M120ox
t2 = His.rt_ref(1)-delta;
[rts3,top1_rt3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end

end % relocate

% =========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
% DA relocation (MS2-assisted): same windows as MS1-only, but calls get_rts2
% to leverage diagnostic fragments and reduce false positives.

delta  = 0.1;
nsplit = 1;

% -------------------------------- M120ox ---------------------------------
hno = 2;
t1 = His.rt_ref(1)-10;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% -------------------------------- K122ac ---------------------------------
hno = 3;
t1 = max([His.rt_ref(1)-8, His.rt_ref(2)+delta]);
t2 = His.rt_ref(1)-delta;
[rts3,top1_rt3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end

end % relocate2
