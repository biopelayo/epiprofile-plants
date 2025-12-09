function H4_02b_20_35(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ========================================================================
% H4_02b_20_35 — Targeted panel for histone H4 peptide "KVLRDNIQGITKPAIR"
%                 (aa 20–35) focusing on K20 PTMs: me1, me2, me3, ac.
% ========================================================================
% Workflow:
%   1) Build panel definition (His) with sequence, PTM rows, charge grid,
%      theoretical m/z and seed RTs.
%   2) Anchor the unmodified (derivatized) peptide via check_ref; optionally
%      refine by scanning with get_rts2 in DA mode.
%   3) Extract unmodified with get_histone0 and correct global drift (shift
%      all PTM seed RTs by the observed delta).
%   4) Relocate PTM windows (MS1-only: relocate ; DA: relocate2).
%   5) Extract modified rows with get_histone1 in a loop.
%   6) Save matrices, draw summary plots, and optionally export PSMs.
%
% INPUTS
%   MS1_index   : [nMS1 x 2] scan index and RT (min)
%   MS1_peaks   : core-specific container of centroided MS1 peaks
%   MS2_index   : [nMS2 x 2] (used if DA relocation or PSM export)
%   MS2_peaks   : core-specific container of centroided MS2 peaks
%   ptol        : mass tolerance (e.g., ppm) passed to core extractors
%   cur_outpath : output directory
%   special     : struct with switches (ndebug, nDAmode, nhmass, raw_path)
%
% OUTPUTS (to disk)
%   <cur_outpath>/H4_02b_20_35.mat with His, pep_rts [5x3], pep_intens [5x3],
%   mono_isointens [nMS1 x 5]; plus figures; plus PSM tables if nDAmode==1.
% ========================================================================

% ------------------------------
% Idempotency guard
% ------------------------------
out_filename = 'H4_02b_20_35';
% fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% ------------------------------
% Panel definition
% ------------------------------
His = init_histone();

% ------------------------------
% Compute layout: anchor, relocate, extract
% ------------------------------
unitdiff = 1.0032; % 13C spacing for monoisotopic tracing utilities
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);

% ------------------------------
% Persist results
% ------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% ------------------------------
% Plotting (XIC/layout)
% ------------------------------
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------
% Optional PSM export (historical semantic: nDAmode==1)
% ------------------------------
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

end % main


% ========================================================================
% init_histone — Define peptide, PTMs, charge grid, theoretical m/z, RTs
% ========================================================================
function His = init_histone()
%%

% Long H4 peptide spanning aa 20–35 (K20 at first peptide position)
His.pep_seq = 'KVLRDNIQGITKPAIR';

% PTM rows on K20 plus derivatization marks:
His.mod_short = {'unmod';
    'K20me1';
    'K20me2';
    'K20me3';
    'K20ac'};

% Compact PTM encoding (peptide-relative indices):
%  - 0 = peptide N-term (propionylation from sample prep)
%  - 1 = Lys20 within this peptide (target PTM site)
%  - 12 = internal Lys in the peptide (propionylated per prep), per core spec
His.mod_type = {'0,pr;1,pr;12,pr;';
    '0,pr;1,me1;12,pr;';
    '0,pr;1,me2;12,pr;';
    '0,pr;1,me3;12,pr;';
    '0,pr;1,ac;12,pr;'};

% Charge states: 2+, 3+, 4+ (longer peptide; multi-charge dominates)
His.pep_ch = repmat([2 3 4],length(His.mod_type),1);

% Theoretical m/z values for each row×charge (computed from sequence+mods)
His.pep_mz = calculate_pepmz(His);

% Seed RTs (minutes) — starting points for relocation windows
His.rt_ref = [42.43
    43.40
    38.03
    37.28
    41.13];

% Display flags (0 = hidden in plots); this is an auxiliary/confirmatory panel
His.display = zeros(length(His.mod_type),1);

% Cosmetic: ensure the first charge column corresponds to 2+
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = 1:npep; % reorder ALL rows, as in original logic
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end;

end % init_histone


% ========================================================================
% calculate_layout — Anchor unmod, correct drift, relocate, extract PTMs
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

% ---- 1) Anchor unmodified with check_ref (+ DA scan if applicable)
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only anchoring
        [His.rt_ref(1),special.ndebug] = ...
            check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                      His.rt_ref(1),special.ndebug);
        % Fallback strategy if check_ref signals ndebug==2
        if 2==special.ndebug
            % Scan unmod across the whole gradient
            hno = 1; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = ...
                get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);
            % Scan K20me1 across the whole gradient (reference to place unmod)
            hno = 2; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = ...
                get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok
            if 0==isempty(top1_rt2)
                if 1==isempty(top1_rt1)
                    % If unmod is unclear, place it ~3 min before me1 apex
                    His.rt_ref(1) = top1_rt2-3;
                else
                    % If candidate unmod peaks exist, pick strongest within
                    % (me1−18 .. me1) window
                    p = find(rts1>top1_rt2-18 & rts1<top1_rt2);
                    if 0==isempty(p)
                        [tmp,pp] = max(inten_sum1(p)); %#ok
                        His.rt_ref(1) = rts1(p(pp));
                    end;
                end;
            else
                % If me1 is missing but unmod has a candidate, use it
                if 0==isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt1;
                end;
            end;
        end;
    else
        % DA anchoring: refine with check_ref then probe ±5 min (or full run)
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

% Extract unmodified row (also produces monoisotopic trace for diagnostics)
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% If unmod is found, correct global drift and store matrices
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;

% ---- 2) Relocate PTM windows (MS1-only vs DA)
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end;
end;

% ---- 3) Extract modified rows sequentially (me1, me2, me3, ac)
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
% relocate — MS1-only relocation strategy for K20 PTMs
% ========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%

delta = 0.1; % avoid capturing tailing of the anchor
nsplit = 1;  % passthrough flag to get_rts (core-specific behavior)

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

% Resolve me2 vs me3 together using intensity patterns
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
% relocate2 — MS2-assisted relocation (same windows, uses get_rts2)
% ========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%

delta = 0.1;
nsplit = 1;

% K20me1 (row 2) — DA path slightly tighter upper bound (+16)
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

% K20me2 (row 3)
hno = 3;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K20me3 (row 4)
hno = 4;
t1 = 6;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% Joint resolution of me2 vs me3
[His.rt_ref(3),His.rt_ref(4)] = ...
    find_pair(rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
              rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K20ac (row 5)
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
