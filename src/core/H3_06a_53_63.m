function H3_06a_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H3_06a_53_63
% -------------------------------------------------------------------------
% Purpose:
%   Quantify the H3 peptide RYQKSTELLIR (53–63) and two propionylated
%   variants (S57pr, Y54pr) using EpiProfile-style logic:
%     - Anchor the unmodified form (derivatized termini).
%     - Calibrate RT seeds by the observed drift.
%     - Relocate candidates in MS1-only or MS2-assisted (DA) modes.
%     - Disambiguate S57pr vs Y54pr either by paired-peak logic or by an
%       MS1/MS2 joint extractor when they are not well separated.
%
% Inputs:
%   MS1_index, MS1_peaks: centroided MS1 data in core-internal format.
%   MS2_index, MS2_peaks: centroided MS2 data (required for DA relocation).
%   ptol: mass tolerance (ppm or instrument-specific units).
%   cur_outpath: output folder path.
%   special: struct with control flags and method hints:
%       .ndebug   -> if 1, use relocateD (diagnostic relocation).
%       .nDAmode  -> relocation/PSM semantics (historical):
%                      - if 2, use DA relocation (relocate2 + get_rts2)
%                      - if 1, export PSMs after quantification
%       .nhmass   -> mass hint for get_rts2 (neutral loss / precursor).
%       .raw_path -> vendor path, used by check_ref to refine the unmod RT.
%
% Outputs (to disk):
%   - <cur_outpath>/H3_06a_53_63.mat with pep_rts, pep_intens, mono_isointens
%   - Figures drawn by draw_layout (layout + XIC summary)
%   - Optional PSM tables if special.nDAmode == 1
%
% Notes:
%   - This function assumes the core’s mod_type indexing. If your build
%     differs, update init_histone() accordingly.
%   - Charge columns are reordered to put 2+ as the first column to keep a
%     consistent convention across panels.
% -------------------------------------------------------------------------

% ------------------------------- check -----------------------------------
out_filename = 'H3_06a_53_63';
% fprintf(1,'%s..',out_filename);  % optional progress print
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    % If a .mat result already exists, skip computation (idempotent run).
    return;
end

% -------------------------------- init -----------------------------------
His = init_histone(cur_outpath,out_filename);

% ------------------------------ calculate --------------------------------
unitdiff = 1.0032;               % C13 isotopic spacing used by mono-tracing
Mods = GetMods();                % PTM mass table (core function, not here)
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

% -------------------------------- output ---------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% --------------------------------- draw ----------------------------------
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2); % RT axis for plotting (from index table)
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
            mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------- Get PSM ---------------------------------
if 1==special.nDAmode
    % Historical semantics: when nDAmode==1, export PSMs after quant.
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
           mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % main

% =========================================================================
function His = init_histone(cur_outpath,out_filename)
%%
% Build the "His" structure with sequence, PTM compact encoding, charge
% grid, m/z table, reference RTs, display flags, and output bookkeeping.

His.pep_seq = 'RYQKSTELLIR';

% Short PTM labels used downstream (rows correspond 1:1 to mod_type):
His.mod_short = {'unmod';
                 'S57pr';
                 'Y54pr'};

% Compact PTM encoding (core-dependent indexing!):
%  - 'pr' -> propionylation
%  - indices (0,2,4,5,...) refer to core’s internal residue map for this panel
His.mod_type = {'0,pr;4,pr;';
                '0,pr;4,pr;5,pr;';
                '0,pr;2,pr;4,pr;'};

% Charge grid replicated per row (allow 1+..4+ for this medium-length peptide)
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);

% (Legacy hard-coded masses kept as comment for reference)
%{
His.pep_mz = [1518.8639  759.9356  506.9595  380.4714;
              1574.8901  787.9487  525.6349  394.4780;
              1574.8901  787.9487  525.6349  394.4780];
%}

% Compute m/z from sequence + mod_type via core helper, keeping the panel portable
His.pep_mz = calculate_pepmz(His);

% Seed RTs (minutes). These are used for initial windows and then drift-calibrated.
His.rt_ref = [37;
              41.58;
              42];

% Display flags per row (0=hide in default plots, 1=show)
% This panel hides all by default to avoid clutter; quant still occurs.
His.display = zeros(length(His.mod_type),1);

% Bookkeeping for downstream writers/drawers
His.outpath = cur_outpath;
His.outfile = out_filename;

% ----------------------- charge column reordering ------------------------
% Convention: ensure that the first column corresponds to 2+ (main_ch).
% This harmonizes merging across panels and simplifies plotting.
main_ch = His.pep_ch(1,2); % 2+ sits in the middle for [1 2 3 4]
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
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special)
%%
% High-level orchestrator:
%   1) Refine the unmodified anchor RT (check_ref +/- DA probe).
%   2) Extract unmodified (get_histone0) and apply drift calibration.
%   3) Relocate modified rows (MS1-only or DA).
%   4) Extract modifieds:
%        - If seeds are well separated (>0.4 min), use get_histone10 per row.
%        - Else use get_histone2 (joint resolution of the pair).
%
% Returns:
%   pep_rts (nPTM x nCharge), pep_intens (nPTM x nCharge),
%   mono_isointens (nMS1 x nPTM)

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% ---------------------------- anchor: unmod -------------------------------
His.rt_unmod_orig = His.rt_ref(1);

if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only: check_ref may refine the anchor and also flip ndebug
        [His.rt_ref(1),special.ndebug] = check_ref( ...
            special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
    else
        % DA relocation mode: still call check_ref, but also consider probing with MS2
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref( ...
            special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);

        % If anchor did not move, broaden a symmetric window and probe with get_rts2
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
        else
            delta = 5;
            t1 = His.rt_ref(1)-delta;
            t2 = His.rt_ref(1)+delta;
        end
        hno = 1; % unmod row
        % get_rts2 returns candidate RTs and a top-RT; we keep the top if present
        [~,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                                 ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% Extract unmodified with the standard helper (MS1-focused)
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% -------------------------- drift calibration ----------------------------
if cur_rts(1)>0
    % If unmod was detected, adopt its observed RT for row 1 and propagate
    % the global shift to other reference RTs (robust to LC drift).
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    % Fill outputs for unmod
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% ------------------------------ relocation -------------------------------
if 1==special.ndebug
    % Diagnostic relocation (minimal changes, useful to inspect seeds)
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        % MS1-only relocation (get_rts)
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        % DA relocation (get_rts2 uses MS2 diagnostics)
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special.nhmass);
    end
end

% ----------------------- final extraction of modifieds --------------------
% If the seed separation between rows 3 (Y54pr) and 2 (S57pr) is large
% enough, quantify them independently; otherwise, use a joint resolver.
if His.rt_ref(3)-His.rt_ref(2)>0.4
    for hno=2:3
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone10( ...
            MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end
    end
else
    % Joint extraction (MS1/MS2-aware) to deconvolve co-eluting pair
    hno = 2;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
        MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end
end

end % calculate_layout

% =========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
% MS1-only relocation strategy:
%   - Search S57pr over [unmod+0.1, unmod+30] minutes.
%   - If two candidate maxima are present (2nd >= 1/30 of 1st; |ΔRT| < 5),
%     assign earlier RT to S57pr and later RT to Y54pr.
%   - Otherwise, fix S57pr to the top RT and adjust Y54pr by propagating
%     the observed delta with respect to its previous seed.

delta = 0.1;     % avoid capturing anchor tailing
nsplit = 0;      % no sub-splitting for this search

% ------------------------------ S57pr row --------------------------------
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts2,top1_rt2,inten_sum2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% Sort candidates by total intensity (descending)
[tmp_sum,ix] = sort(inten_sum2,'descend');
tmp_rts = rts2(ix);

% If a second candidate is reasonably strong and close in RT, treat the
% pair as S57pr (earlier) and Y54pr (later).
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<5
    His.rt_ref(hno)   = min([tmp_rts(2),tmp_rts(1)]); % S57pr earliest
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]); % Y54pr latest
else
    % Fallback: fix S57pr to the main candidate (if any) and propagate delta to Y54pr
    old_t = His.rt_ref(hno);
    if 1==isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end

    % ------------------------------ Y54pr row -----------------------------
    hno = 3;
    if 0==His.rt_ref(hno-1)
        % If S57pr vanished, Y54pr cannot be safely placed
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1)
        % If S57pr moved from its old seed, shift Y54pr by the same delta
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end
end

end % relocate

% =========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
% DA relocation strategy (uses get_rts2 for MS2-assisted RT confirmation).
% Logic mirrors relocate(), but leverages MS2 to improve robustness.

delta = 0.1;
nsplit = 0;

% ------------------------------ S57pr row --------------------------------
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts2,top1_rt2,inten_sum2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                                      ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[tmp_sum,ix] = sort(inten_sum2,'descend');
tmp_rts = rts2(ix);

if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<5
    His.rt_ref(hno)   = min([tmp_rts(2),tmp_rts(1)]); % S57pr earliest
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]); % Y54pr latest
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end

    % ------------------------------ Y54pr row -----------------------------
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1)
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end
end

end % relocate2
