function H4_02a_18_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ========================================================================
% H4_02a_18_23 — Targeted panel for histone H4 (aa 18–23) peptide "HRKVLR"
%                 and its K20 PTMs: me1, me2, me3, ac (propionylated prep).
% ========================================================================
% This function mirrors your original logic and structure but adds
% exhaustive comments for clarity. It anchors the unmodified (chemically
% derivatized) peptide, relocates RT windows for K20me1/2/3 and K20ac,
% extracts per-charge XICs, writes outputs, draws plots, and optionally
% exports PSMs.
%
% INPUTS
%   MS1_index   : [nMS1 x 2] matrix with scan index and retention time (min)
%   MS1_peaks   : core-specific MS1 centroided peaks container
%   MS2_index   : [nMS2 x 2] (only needed for DA relocation or PSM export)
%   MS2_peaks   : core-specific MS2 centroided peaks container
%   ptol        : mass tolerance passed to the core extractors (e.g., ppm)
%   cur_outpath : output directory for .mat and figures
%   special     : struct with switches/hints:
%                 .ndebug   -> 1: diagnostic relocation (relocateD)
%                 .nDAmode  -> 2: use MS2-assisted relocation (relocate2/get_rts2)
%                                1: export PSM after quant (historical semantic)
%                                else: MS1-only relocation (relocate/get_rts)
%                 .nhmass   -> neutral loss / precursor mass hint for get_rts2
%                 .raw_path -> raw vendor path for check_ref (unmod anchoring)
%
% OUTPUTS (to disk)
%   <cur_outpath>/H4_02a_18_23.mat with:
%     His, pep_rts [5x3], pep_intens [5x3], mono_isointens [nMS1 x 5]
%   plus figures from draw_layout(...)
%   plus PSM tables if special.nDAmode == 1 (GetPSM)
% ========================================================================

% ------------------------------
% Idempotency guard: skip if a previous result exists
% ------------------------------
out_filename = 'H4_02a_18_23';
% fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% ------------------------------
% Build the histone peptide/PTM definition (His struct)
% ------------------------------
His = init_histone();

% ------------------------------
% Main computation: anchor, relocate, extract, save & draw
% ------------------------------
unitdiff = 1.0032; % 13C spacing used by monoisotopic tracing in helpers
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);

% Persist quantitative matrices and His definition
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% Prepare RT axis for plotting utilities (often called "isorts")
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);

% Draw panel summary/XICs using your core’s plotting routine
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% Optional PSM export (historical switch semantics)
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

end % main function


% ========================================================================
% init_histone — Define peptide, PTM rows, charge grid, theoretical m/z,
%                seed RTs and display flags.
% ========================================================================
function His = init_histone()
%%

% Peptide sequence (H4 18–23, includes K20 as the third residue in "HRKVLR")
His.pep_seq = 'HRKVLR';

% PTM rows: unmodified (propionylated prep), K20me1, K20me2, K20me3, K20ac
His.mod_short = {'unmod';
    'K20me1';
    'K20me2';
    'K20me3';
    'K20ac'};

% Compact PTM encoding (relative to peptide positions):
%  - 0 = peptide N-term (pr from derivatization)
%  - 3 = K20 within "HRKVLR"
His.mod_type = {'0,pr;3,pr;';
    '0,pr;3,me1;';
    '0,pr;3,me2;';
    '0,pr;3,me3;';
    '0,pr;3,ac;'};

% Charge states: consider 1+, 2+, and 3+ (slightly longer peptide; 3+ possible)
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);

% Theoretical m/z per row × charge (computed; portable across forks)
% (Legacy hard-coded table is commented in your original; we compute instead)
His.pep_mz = calculate_pepmz(His);

% Seed RTs (minutes) used to build relocation windows later
His.rt_ref = [21.85
    24.73
    14.33
    14.18
    20.15];

% Display flags (all rows hidden by default in figures for this auxiliary panel)
His.display = zeros(length(His.mod_type),1);

% Optional: reorder charge columns to put 2+ first (cosmetic consistency)
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    % Reorder for ALL rows (tune = 1:npep) as per your original code
    tune = 1:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end;

end % init_histone


% ========================================================================
% calculate_layout — Anchor unmodified, correct drift, relocate windows,
%                    extract modified rows, collect matrices.
% ========================================================================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special)
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);

% Preallocate output containers
pep_rts        = zeros([npep,ncharge]);
pep_intens     = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% ------------------ 1) Anchor "unmod" with check_ref (+ optional DA probe)
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only relocation path
        [His.rt_ref(1),special.ndebug] = ...
            check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                      His.rt_ref(1),special.ndebug);
        % Special fallback: if check_ref signals ndebug==2, try to infer
        % unmod position by also scanning K20me1 and aligning unmod before it.
        if 2==special.ndebug
            % Try to find unmod across the whole gradient
            hno = 1; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = ...
                get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);
            % Also search K20me1 across the whole gradient
            hno = 2; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = ...
                get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok
            if 0==isempty(top1_rt2)
                if 1==isempty(top1_rt1)
                    % If unmod peak is not obvious, place it ~3 min before me1
                    His.rt_ref(1) = top1_rt2-3;
                else
                    % If unmod has candidates, favor the strongest within
                    % a window just before me1 (−18..0 min)
                    p = find(rts1>top1_rt2-18 & rts1<top1_rt2);
                    if 0==isempty(p)
                        [tmp,pp] = max(inten_sum1(p)); %#ok
                        His.rt_ref(1) = rts1(p(pp));
                    end;
                end;
            else
                % If me1 not found, but unmod is, use unmod directly
                if 0==isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt1;
                end;
            end;
        end;
    else
        % DA relocation path: refine with check_ref then probe with get_rts2
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                                  His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0; t2 = MS1_index(num_MS1,2);
        else
            delta = 5; t1 = His.rt_ref(1)-delta; t2 = His.rt_ref(1)+delta;
        end;
        hno = 1;
        [rts1,top1_rt1] = ...
            get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end;
    end;
end;

% Extract unmodified row using the core helper for anchors
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% If unmod is found on at least one charge, correct global drift and store
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% ------------------ 2) Relocate PTM windows (MS1-only vs DA) ------------------
if 1==special.ndebug
    % Minimal/diagnostic relocation (implementation in your core)
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end;
end;

% ------------------ 3) Extract modified rows (me1, me2, me3, ac) -------------
for hno=2:5
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)        = cur_rts;
        pep_intens(hno,1:ncharge)     = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;

end % calculate_layout


% ========================================================================
% relocate — MS1-only relocation windows relative to the unmodified anchor
% ========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1; % avoid tailing of the anchor
nsplit = 1;  % passed to get_rts; splitting behavior is core-dependent

% --- K20me1 (row 2): later than unmod
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+18;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% --- K20me2 (row 3): earlier than unmod
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% --- K20me3 (row 4): also earlier than unmod
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% Resolve me2 vs me3 jointly (considers intensity sums)
[His.rt_ref(3),His.rt_ref(4)] = ...
    find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
              rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% --- K20ac (row 5): between me2 and unmod
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;

end % relocate


% ========================================================================
% relocate2 — MS2-assisted relocation (same windows; get_rts2 under the hood)
% ========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K20me1 (slightly tighter upper bound vs MS1 path: +16 here)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end;

% K20me2
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K20me3
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% Resolve me2 vs me3 jointly
[His.rt_ref(3),His.rt_ref(4)] = ...
    find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
              rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K20ac
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts5)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt5;
end;
end % relocate2
