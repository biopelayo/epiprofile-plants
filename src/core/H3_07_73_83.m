function H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H3_07_73_83
% -------------------------------------------------------------------------
% Purpose:
%   Quantify the H3 peptide EIAQDFKTDLR (73–83) and its K79 variants:
%     - unmod (derivatized termini per EpiProfile conventions)
%     - K79me1 / K79me2 / K79me3
%     - K79ac
%   Strategy:
%     1) Refine the unmodified anchor RT.
%     2) Extract unmodified and apply global drift calibration.
%     3) Relocate modified rows (MS1-only or MS2-assisted (DA)).
%     4) Extract modifieds (get_histone1).
%     5) Save, draw, optionally export PSMs.
%
% Inputs:
%   MS1_index, MS1_peaks: centroided MS1 in core-internal format.
%   MS2_index, MS2_peaks: centroided MS2 (required for DA relocation).
%   ptol: mass tolerance (ppm or instrument-specific units).
%   cur_outpath: output directory for .mat and figures.
%   special: struct with flags/hints:
%       .ndebug   -> 1: use relocateD; 0: normal relocation.
%       .nDAmode  -> relocation/PSM semantics (historical):
%                      2 => DA relocation (relocate2 + get_rts2)
%                      1 => export PSMs after quantification
%       .nhmass   -> mass hint for get_rts2 (neutral loss / precursor).
%       .raw_path -> vendor raw path (used by check_ref).
%
% Outputs (to disk):
%   <cur_outpath>/H3_07_73_83.mat with pep_rts, pep_intens, mono_isointens
%   Figures via draw_layout; optional PSM tables if special.nDAmode == 1.
% -------------------------------------------------------------------------

% ------------------------------- check -----------------------------------
out_filename = 'H3_07_73_83';
fprintf(1,'%s..',out_filename);  % progress print (kept as in original)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    % Idempotent run: if results exist, skip.
    return;
end

% -------------------------------- init -----------------------------------
His = init_histone();

% ------------------------------ calculate --------------------------------
unitdiff = 1.0032;  % 13C spacing used for monoisotopic tracing
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special);

% -------------------------------- output ---------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% --------------------------------- draw ----------------------------------
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2); % RT axis from MS1 index table
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
            mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------- Get PSM ---------------------------------
if 1==special.nDAmode
    % Historical semantics: nDAmode==1 triggers PSM export after quant.
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
           mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % main

% =========================================================================
function His = init_histone()
%%
% Build the "His" structure with sequence, PTMs, charge grid, m/z table,
% reference RTs, display flags. Core-dependent indexing applies.

His.pep_seq = 'EIAQDFKTDLR';

His.mod_short = {'unmod';
                 'K79me1';
                 'K79me2';
                 'K79me3';
                 'K79ac'};

% Compact PTM encoding (do NOT change unless your core uses different
% residue indices for this panel):
His.mod_type = {'0,pr;7,pr;';
                '0,pr;7,me1;';
                '0,pr;7,me2;';
                '0,pr;7,me3;';
                '0,pr;7,ac;'};

% Allow charges 1+..4+ for this length; replicate per row
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);

% Legacy hard-coded masses kept as comment for reference:
%{
His.pep_mz = [1447.7427  724.3750  483.2524  362.6911;
              1461.7584  731.3828  487.9243  366.1951;
              1419.7478  710.3775  473.9208  355.6924;
              1433.7635  717.3854  478.5927  359.1963;
              1433.7271  717.3672  478.5805  359.1872];
%}

% Compute portable m/z values
His.pep_mz = calculate_pepmz(His);

% Seed RTs (minutes): unmod, me1, me2, me3, ac
His.rt_ref = [40.37;
              41.20;
              33.92;
              33.72;
              39.17];

% Show all rows in plots by default (quantification unaffected)
His.display = ones(length(His.mod_type),1);

% ----------------------- charge column reordering ------------------------
% Convention across panels: put 2+ as first column to harmonize merges
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
% Orchestrates:
%   - Anchor refinement (check_ref; DA probe optional).
%   - Unmod extraction -> drift calibration.
%   - Relocation (MS1-only vs DA).
%   - Extraction of modified rows with get_histone1.

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
        % If the anchor did not move, probe a broader symmetric window using DA
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
        else
            delta = 5;
            t1 = His.rt_ref(1)-delta;
            t2 = His.rt_ref(1)+delta;
        end
        hno = 1; % unmod
        % NOTE: original code calls get_rts22 (two '2'). If your core exposes
        % get_rts2 (common), adjust accordingly. We keep the original call.
        [~,top1_rt1] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                                  ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% Extract unmodified row and optionally calibrate seeds by observed drift
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
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special.nhmass);
    end
end

% ----------------------------- final extraction --------------------------
% K79me1, K79me2, K79me3, K79ac (rows 2..5)
for hno=2:5
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
% MS1-only relocation windows per K79 PTM:
%   - me1: [unmod+0.1, unmod+14]
%   - me2: [6, unmod-3]
%   - me3: [6, unmod-3] (paired with me2 via find_pair)
%   - ac : [(unmod+me2)/2, unmod-0.1]

delta = 0.1;
nsplit = 1;

% ------------------------------- K79me1 ----------------------------------
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% ------------------------------- K79me2 ----------------------------------
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% ------------------------------- K79me3 ----------------------------------
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% Resolve me2 vs me3 as a pair, using intensity evidence
[His.rt_ref(3),His.rt_ref(4)] = find_pair( ...
    rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
    rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% -------------------------------- K79ac ----------------------------------
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end

end % relocate

% =========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
% DA relocation (MS2-assisted). Same windows, calls get_rts2.

delta = 0.1;
nsplit = 1;

% ------------------------------- K79me1 ----------------------------------
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+14;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% ------------------------------- K79me2 ----------------------------------
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts2( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% ------------------------------- K79me3 ----------------------------------
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts2( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[His.rt_ref(3),His.rt_ref(4)] = find_pair( ...
    rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
    rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% -------------------------------- K79ac ----------------------------------
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end

end % relocate2
