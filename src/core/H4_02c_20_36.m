function H4_02c_20_36(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ========================================================================
% H4_02c_20_36 — Targeted panel for histone H4 peptide "KVLRDNIQGITKPAIRR"
%                 (aa 20–36), monitoring K20 PTMs (me1, me2, me3, ac).
% ------------------------------------------------------------------------
% This function quantifies one long H4 peptide (20–36, Arg-extended) and
% its K20 modification states. It follows the core workflow used across
% your panels:
%   • Build panel (sequence, PTMs, charges, theoretical m/z, seed RTs)
%   • Anchor the unmodified (derivatized) peptide via check_ref
%   • Optionally refine anchor with MS2-assisted scan (DA mode)
%   • Extract unmodified (get_histone0) and correct global RT drift
%   • Relocate PTM windows (MS1-only or DA) and extract (get_histone1)
%   • Persist matrices, draw layouts, and optionally export PSMs
%
% Inputs:
%   MS1_index   [nMS1 x 2]  : scan index and retention time (min)
%   MS1_peaks                : centroided MS1 peaks (core-defined struct)
%   MS2_index, MS2_peaks     : centroided MS2 (for DA relocation/PSM)
%   ptol                     : mass tolerance (e.g., ppm)
%   cur_outpath              : output directory
%   special                  : struct with switches/paths:
%       .ndebug   -> 1 for minimal diagnostics (relocateD), otherwise relocate
%       .nDAmode  -> 2: use MS2-assisted relocation (get_rts2/relocate2)
%                      1: export PSM after quantification (historical)
%                      else: MS1-only relocation (get_rts/relocate)
%       .nhmass   -> neutral loss / precursor hint for DA (if applicable)
%       .raw_path -> path to raw file for check_ref anchoring
%
% Outputs to disk:
%   <cur_outpath>/H4_02c_20_36.mat with fields:
%     - His (panel definition)
%     - pep_rts [5 x 3]     : RT per PTM row x charge (2+/3+/4+)
%     - pep_intens [5 x 3]  : integrated intensities aligned to pep_rts
%     - mono_isointens      : monoisotopic traces (nMS1 x 5)
%   Plus layout figures and (optionally) PSM tables.
% ========================================================================

% ------------------------------
% 1) Idempotent check: skip if .mat already exists
% ------------------------------
out_filename = 'H4_02c_20_36';
% fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return; % respect original behavior: do nothing if results exist
end

% ------------------------------
% 2) Panel initialization (sequence, PTM rows, charges, m/z, RT seeds)
% ------------------------------
His = init_histone();

% ------------------------------
% 3) Compute layout: anchor, relocate, extract
% ------------------------------
unitdiff = 1.0032; % ~13C mass spacing, used by monoisotopic helpers
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);

% ------------------------------
% 4) Persist matrices to disk
% ------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% ------------------------------
% 5) Draw summary plots / XIC layouts for QC
% ------------------------------
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2); % RT axis
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------
% 6) Optional PSM export (historical semantic: nDAmode==1)
% ------------------------------
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % function H4_02c_20_36


% ========================================================================
% init_histone — Define peptide, PTM catalog, charge grid, m/z, seed RTs
% ========================================================================
function His = init_histone()
%% Panel definition (exactly as in your original code)

% Peptide (H4 aa 20–36), Arg-extended vs 20–35
His.pep_seq = 'KVLRDNIQGITKPAIRR';

% PTM states focused on K20 (first residue in this peptide)
His.mod_short = {'unmod';
    'K20me1';
    'K20me2';
    'K20me3';
    'K20ac'};

% Compact PTM encoding per peptide-relative indices:
%  - 0 : peptide N-terminus (propionylation from sample prep)
%  - 1 : K20 site (me1/me2/me3/ac)
%  - 12: internal Lys propionylated (per your core convention)
His.mod_type = {'0,pr;1,pr;12,pr;';
    '0,pr;1,me1;12,pr;';
    '0,pr;1,me2;12,pr;';
    '0,pr;1,me3;12,pr;';
    '0,pr;1,ac;12,pr;'};

% Charge states expected for this long/basic peptide (favoring multi-charge)
His.pep_ch = repmat([2 3 4],length(His.mod_type),1);

% Theoretical m/z (row x charge) computed from sequence + modifications
His.pep_mz = calculate_pepmz(His);

% Seed retention times (min) for relocation windows
His.rt_ref = [39.53
    40.44
    35.63
    34.87
    38.83];

% Hidden by default in plots (auxiliary confirmation panel)
His.display = zeros(length(His.mod_type),1);

% Cosmetic: ensure first charge column matches the "main" (here: 3+ at col 2)
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end
    tune = 1:npep;          % reorder all rows consistently
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end

end % init_histone


% ========================================================================
% calculate_layout — Anchor unmod, drift-correct, relocate, extract PTMs
% ========================================================================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special)
%% Follow original control-flow exactly, with explanatory comments

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);

% Preallocate output matrices (rows: PTM states; cols: charges)
pep_rts        = zeros([npep,ncharge]);
pep_intens     = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% ---- Anchor: unmodified (derivatized) peptide
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only anchoring via check_ref (no MS2 assistance here)
        [His.rt_ref(1),special.ndebug] = ...
            check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                      His.rt_ref(1),special.ndebug);

        % Fallback path if check_ref signals ndebug==2:
        % scan unmod and me1 across the full run to infer an anchor
        if 2==special.ndebug
            % Unmodified scan over the whole gradient
            hno = 1; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = ...
                get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);
            % K20me1 scan likewise
            hno = 2; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = ...
                get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok

            if 0==isempty(top1_rt2)
                % If me1 is present:
                if 1==isempty(top1_rt1)
                    % Place unmod about 3 min before me1 apex when unmod peak is unclear
                    His.rt_ref(1) = top1_rt2-3;
                else
                    % If unmod candidates exist, pick the strongest within (me1−18 .. me1)
                    p = find(rts1>top1_rt2-18 & rts1<top1_rt2);
                    if 0==isempty(p)
                        [~,pp] = max(inten_sum1(p)); %#ok
                        His.rt_ref(1) = rts1(p(pp));
                    end
                end
            else
                % If me1 is missing but unmod has a candidate, use that
                if 0==isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt1;
                end
            end
        end
    else
        % DA anchoring: check_ref then probe ±5 min (or full gradient if unchanged)
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                                  His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0; t2 = MS1_index(num_MS1,2);
        else
            delta = 5; t1 = His.rt_ref(1)-delta; t2 = His.rt_ref(1)+delta;
        end
        hno = 1;
        [rts1,top1_rt1] = ...
            get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% Extract unmodified row; also yields monoisotopic trace for diagnostics
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% If an unmod RT was found, correct global drift and store into outputs
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% ---- Relocate PTM windows according to mode
if 1==special.ndebug
    % Minimal diagnostic relocation (implementation in core)
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        % MS1-only relocation
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        % MS2-assisted relocation (DA mode)
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end
end

% ---- Extract modified rows sequentially (me1, me2, me3, ac)
for hno=2:5
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)        = cur_rts;
        pep_intens(hno,1:ncharge)     = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

end % calculate_layout


% ========================================================================
% relocate — MS1-only relocation windows for K20 PTMs
% ========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
delta  = 0.1; % small offset to avoid tailing around the anchor
nsplit = 1;   % passthrough flag for get_rts (core-specific behavior)

% --- K20me1 (row 2): expected after unmodified
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+18;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% --- K20me2 (row 3): earlier than unmodified
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% --- K20me3 (row 4): also earlier than unmodified
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% Resolve me2 vs me3 jointly using intensity patterns
[His.rt_ref(3),His.rt_ref(4)] = ...
    find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
              rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% --- K20ac (row 5): between me2 and unmodified
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


% ========================================================================
% relocate2 — MS2-assisted relocation (same windows, get_rts2 backend)
% ========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
delta  = 0.1;
nsplit = 1;

% --- K20me1 (row 2) — DA path (narrower upper bound +16)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% --- K20me2 (row 3)
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% --- K20me3 (row 4)
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% Joint resolution (me2 vs me3)
[His.rt_ref(3),His.rt_ref(4)] = ...
    find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
              rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% --- K20ac (row 5)
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
