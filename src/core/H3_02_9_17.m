function H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%% H3_02_9_17 — Targeted layout for H3 peptide KSTGGKAPR (K9 / K14 PTMs)
%
% HIGH-LEVEL PURPOSE
%   Quantify the histone H3 peptide "KSTGGKAPR" (positions 9–17 region;
%   lysines at K9 and K14 are the PTM targets) across a panel of single and
%   combinatorial modifications:
%     - unmod, K9me1, K9me2, K9me3,
%     - K9ac, K14ac,
%     - K9me1K14ac, K9me2K14ac, K9me3K14ac,
%     - K9acK14ac
%   The routine anchors RT on the unmodified peptide, relocates expected RT
%   windows for each PTM (MS1-only or MS2-assisted depending on 'special'),
%   extracts retention times (RTs) and intensities, writes outputs, plots,
%   and optionally exports simple PSM labels for visualization tools.
%
% INPUTS
%   MS1_index : [nMS1×4]  — (MS1_scan, RT[min], cumMS1peakStartIdx, baseline)
%   MS1_peaks : [ΣMS1peaks×2] — concatenated (m/z, intensity) for all MS1 scans
%   MS2_index : [nMS2×8]  — (MS1_scan, MS1_RT, MS2_scan, m/z, z, FragType, MS2peakStartIdx, baseline)
%   MS2_peaks : [ΣMS2peaks×2] — concatenated (m/z, intensity) for all MS2 scans
%   ptol      : ppm tolerance (some helpers normalize 100→10 ppm for legacy behavior)
%   cur_outpath : output folder
%   special   : struct with run options:
%               • nDAmode: 0/1 = MS1-only; 2 = MS2-assisted
%               • ndebug : RT reference checking path switch
%               • nhmass : toggle alternative MS2 key-ion set
%               • raw_path: raw path used by check_ref(...)
%
% OUTPUTS (files)
%   <out>.mat      — quantitative results (per-PTM, per-charge)
%   <out>.plabel   — spectrum labels (if nDAmode==1 in this codebase convention)
%   Figures        — diagnostic plots via draw_layout(...)
%
% DEPENDENCIES (provided elsewhere in the codebase)
%   GetMods, calculate_pepmz, check_ref
%   get_rts, get_rts2, get_rts22
%   get_histone0, get_histone1, get_histone2, get_histone11
%   find_pair_new
%   output_histone, draw_layout, GetPSM
%
% NOTES
%   - 'unitdiff' (~1.0032 Da) is the 13C–12C mass difference used for
%     isotope pattern checks and m/z ladders.
%   - RT logic uses the unmodified form as anchor; modified states are
%     searched in heuristic windows relative to that anchor and to each
%     other (especially me2/me3 and K9/K14 acetyl pairs).
%

% --- 0) Idempotent skip if already processed
out_filename = 'H3_02_9_17';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end

% --- 1) Initialize histone peptide/PTM library (masses, charges, RT priors)
His = init_histone(cur_outpath,out_filename);

% --- 2) Compute per-PTM RT/intensity using MS1 (and optionally MS2)
unitdiff = 1.0032;       % 13C–12C difference
Mods     = GetMods();    % PTM mass table
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

% --- 3) Persist main quantitative outputs
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% --- 4) Plot diagnostic layout
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% --- 5) Optional PSM export (viewer-friendly labels)
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end


% ========================================================================
%                               NESTED HELPERS
% ========================================================================

function His = init_histone(cur_outpath,out_filename)
%% init_histone — Define KSTGGKAPR PTM panel, charge states, theoretical m/z, RT priors
%
% WHAT WE SET HERE
%   - Sequence: 'KSTGGKAPR' (targets at K9 and K14 along the H3 tail context)
%   - PTM panel covering single K9 methyl/acetyl, single K14 acetyl, and
%     combinatorials (K9me{1,2,3} + K14ac; K9acK14ac).
%   - Charges: 1, 2, 3 (columns). We compute theoretical m/z for each.
%   - RT priors (minutes) for each PTM state, later calibrated by observed
%     unmodified RT.
%
His.pep_seq = 'KSTGGKAPR';
His.mod_short = { ...
    'unmod'; ...
    'K9me1'; ...
    'K9me2'; ...
    'K9me3'; ...
    'K9ac'; ...
    'K14ac'; ...
    'K9me1K14ac'; ...
    'K9me2K14ac'; ...
    'K9me3K14ac'; ...
    'K9acK14ac' };

% 'mod_type' uses the internal "pos,mod;" encoding.
% Positions are 1-based over the peptide string 'K S T G G K A P R'
% Here: pos1=K9, pos6=K14 within the peptide context (for this peptide).
His.mod_type = { ...
    '0,pr;1,pr;6,pr;';  ... % unmodified (N-term + both Ks propionylated derivatization)
    '0,pr;1,me1;6,pr;'; ...
    '0,pr;1,me2;6,pr;'; ...
    '0,pr;1,me3;6,pr;'; ...
    '0,pr;1,ac;6,pr;';  ...
    '0,pr;1,pr;6,ac;';  ...
    '0,pr;1,me1;6,ac;'; ...
    '0,pr;1,me2;6,ac;'; ...
    '0,pr;1,me3;6,ac;'; ...
    '0,pr;1,ac;6,ac;'   };

% Consider 1+, 2+, 3+ for each PTM (columns).
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);

% Theoretical m/z per PTM×charge using Mods mass table and amino-acid masses.
His.pep_mz = calculate_pepmz(His);

% Empirical RT priors for the 10 PTM states (minutes).
His.rt_ref = [23.20; 25.91; 16.59; 16.46; 21.90; 21.91; 24.68; 14.95; 14.86; 20.23];
His.display = ones(length(His.mod_type),1);

% Keep some metadata for downstream output naming
His.outpath  = cur_outpath;
His.outfile  = out_filename;

% Align charge ordering: set a "main" charge column (use row1 col2 as reference)
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


function [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
        MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special)
%% calculate_layout — Anchor unmodified RT, relocate PTM windows, extract RT/intensities
%
% STRATEGY
%   1) Anchor the unmodified RT (with optional check_ref and MS2 aid).
%   2) Extract unmodified signal via get_histone0 and calibrate all priors.
%   3) Relocate modified forms using MS1-only (relocate) or MS2-assisted
%      (relocate2).
%   4) Fetch final RT/intensity for each PTM using get_histone1/2/11.
%
[npep,ncharge]    = size(His.pep_mz);
num_MS1           = size(MS1_index,1);
pep_rts           = zeros([npep,ncharge]);
pep_intens        = zeros([npep,ncharge]);
mono_isointens    = zeros([num_MS1,npep]);

% --- 1) Unmodified RT anchoring with heuristics
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only path with optional reference check
        [His.rt_ref(1),special.ndebug] = check_ref( ...
            special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if 2==special.ndebug
            % Wide search for unmod (hno=1) and me1 (hno=2) to triangulate anchor
            hno = 1; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2);
            hno = 2; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok<ASGLU>

            % Heuristic: if multiple large MS1 maxima exist, pick an early-pair
            % compromise; otherwise use me1 proximity or fallback to unmod apex.
            [tmp_sum,ix] = sort(inten_sum1,'descend'); %#ok<NASGU>
            tmp_rts = rts1(ix);
            if length(tmp_sum)>=3 && tmp_sum(3)>=tmp_sum(1)/30
                new_rts = sort(tmp_rts(1:3));
                if new_rts(2)-new_rts(1)>new_rts(3)-new_rts(2)
                    His.rt_ref(1) = new_rts(2);
                else
                    His.rt_ref(1) = new_rts(1);
                end
            elseif length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<7
                His.rt_ref(1) = min([tmp_rts(2),tmp_rts(1)]);
            else
                if ~isempty(top1_rt2)
                    if isempty(top1_rt1)
                        His.rt_ref(1) = top1_rt2-3;
                    else
                        p = find(rts1>top1_rt2-16 & rts1<top1_rt2);
                        if ~isempty(p)
                            [~,pp] = max(inten_sum1(p)); %#ok<ASGLU>
                            His.rt_ref(1) = rts1(p(pp));
                        end
                    end
                else
                    if ~isempty(top1_rt1)
                        His.rt_ref(1) = top1_rt1;
                    end
                end
            end
        end
    else
        % MS2-assisted unmodified anchoring
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0; t2 = MS1_index(num_MS1,2);
        else
            delta = 5; t1 = His.rt_ref(1)-delta; t2 = His.rt_ref(1)+delta;
        end
        hno = 1;
        [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                                   ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok<ASGLU>
        if ~isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% --- 2) Extract unmodified quantitative trace and calibrate all RT priors
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

if cur_rts(1)>0
    His.rt_ref(1)        = cur_rts(1);
    delta                = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end)    = His.rt_ref(2:end) + delta;
    pep_rts(hno,1:ncharge)            = cur_rts;
    pep_intens(hno,1:ncharge)         = cur_intens;
    mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
end

% --- 3) Relocation: define final expected windows for PTMs
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His); % external debug path
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His); % MS1-only
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass); % MS2-assisted
    end
end

% --- 4) Fetch final RT/intensity for each PTM

% K9 methyls (single K9me1/2/3)
for hno = 2:4
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)        = cur_rts;
        pep_intens(hno,1:ncharge)     = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

% K9ac / K14ac (some implementations treat these with a specialized extractor)
hno = 5; % starting at K9ac; K14ac is #6
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    % get_histone2 returns a 2×ncharge block for (K9ac, K14ac)
    pep_rts(hno:hno+1,1:ncharge)           = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge)        = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1)    = cur_mono_isointens(:,1:2);
end

% Combinatorials with K14ac: K9me1/2 + K14ac and K9acK14ac (skip me3K14ac special below)
for hno = [7 8 10]
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)        = cur_rts;
        pep_intens(hno,1:ncharge)     = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

% K9me3K14ac (some datasets need a slightly different XIC logic)
hno = 9;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end


function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%% relocate — MS1-only relocation of expected RTs for K9/K14 PTMs
%
% Heuristic windows relative to the unmodified anchor:
%   - K9me1 after unmod (anchor+0.5 .. anchor+16)
%   - K9me2/me3 before unmod (6 .. anchor-3), paired by find_pair_new
%   - K9ac between me2 and unmod; K14ac is adjusted relative to K9ac
%   - Combinatorials constrained between their single-PTM counterparts
%
delta  = 0.5;
nsplit = 1;

% -- K9me1 (after unmod)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
His.rt_ref(hno) = iffEmptyZero(rts2, top1_rt2);

% -- K9me2 (before unmod)
hno = 3;
t1 = 6; t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% -- K9me3 (also before unmod); pick pair with me2 (order may swap)
hno = 4;
t1 = 6; t2 = His.rt_ref(1)-3;
[rts4,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok<ASGLU>

[His.rt_ref(4),His.rt_ref(3)] = find_pair_new(top1_rt4, rts3, top1_rt3, inten_sum3, 0);

% -- K9ac between me2 and unmod
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
old_t = His.rt_ref(hno);
His.rt_ref(hno) = iffEmptyZero(rts5, top1_rt5);

% -- K14ac tracks relative shift vs. K9ac (if K9ac moved)
hno = 6;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1)
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end

% -- K9me1K14ac (after unmod but before K9me1)
hno = 7;
t1 = His.rt_ref(1)+delta;
if 0==His.rt_ref(2)
    t2 = His.rt_ref(1)+9;
else
    t2 = His.rt_ref(2)-delta;
end
[rts7,top1_rt7] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
His.rt_ref(hno) = iffEmptyZero(rts7, top1_rt7);

% -- K9me2K14ac (before unmod, near me2)
hno = 8;
t1 = 6;
if 0==His.rt_ref(3)
    t2 = His.rt_ref(1)-4;
else
    t2 = His.rt_ref(3)-delta;
end
[rts8,top1_rt8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok<ASGLU>

% -- K9me3K14ac (before unmod, near me3/me2 pair)
hno = 9;
t1 = 6;
if 0==His.rt_ref(3)
    t2 = His.rt_ref(1)-4;
else
    t2 = His.rt_ref(3)-delta;
end
[rts9,top1_rt9,inten_sum9] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
[His.rt_ref(8),His.rt_ref(9)] = find_pair_new(top1_rt8, rts9, top1_rt9, inten_sum9, 1);

% -- K9acK14ac (between me2 and unmod, but after K9ac if present)
hno = 10;
if 0==His.rt_ref(3)
    t1 = His.rt_ref(1)-18;
else
    t1 = His.rt_ref(3)+delta;
end
if 0==His.rt_ref(5)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(5)-delta;
end
[rts10,top1_rt10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
His.rt_ref(hno) = iffEmptyZero(rts10, top1_rt10);


function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%% relocate2 — MS2-assisted relocation for K9/K14 PTMs
%
% Same windows as 'relocate' but confirmation/selection uses MS2-aided
% routines (get_rts2 / get_rts22) and pairing via find_pair_new.
%
delta  = 0.5;
nsplit = 1;

% -- K9me1
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
His.rt_ref(hno) = iffEmptyZero(rts2, top1_rt2);

% -- K9me2
hno = 3;
t1 = 6; t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% -- K9me3
hno = 4;
t1 = 6; t2 = His.rt_ref(1)-3;
[rts4,top1_rt4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok<ASGLU>
[His.rt_ref(4),His.rt_ref(3)] = find_pair_new(top1_rt4, rts3, top1_rt3, inten_sum3, 0);

% -- K14ac (between me2 and unmod)
hno = 6;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
old_t = His.rt_ref(hno);
His.rt_ref(hno) = iffEmptyZero(rts6, top1_rt6);

% -- K9ac tracks K14ac if it moved
hno = 5;
if 0==His.rt_ref(hno+1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno+1)
    d = His.rt_ref(hno+1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end

% -- K9me1K14ac
hno = 7;
t1 = His.rt_ref(1)+delta;
if 0==His.rt_ref(2)
    t2 = His.rt_ref(1)+9;
else
    t2 = His.rt_ref(2)-delta;
end
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
His.rt_ref(hno) = iffEmptyZero(rts7, top1_rt7);

% -- K9me2K14ac
hno = 8;
t1 = 6;
if 0==His.rt_ref(3)
    t2 = His.rt_ref(1)-4;
else
    t2 = His.rt_ref(3)-delta;
end
[rts8,top1_rt8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok<ASGLU>

% -- K9me3K14ac
hno = 9;
t1 = 6;
if 0==His.rt_ref(3)
    t2 = His.rt_ref(1)-4;
else
    t2 = His.rt_ref(3)-delta;
end
[rts9,top1_rt9,inten_sum9] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
[His.rt_ref(8),His.rt_ref(9)] = find_pair_new(top1_rt8, rts9, top1_rt9, inten_sum9, 1);

% -- K9acK14ac
hno = 10;
if 0==His.rt_ref(3)
    t1 = His.rt_ref(1)-18;
else
    t1 = His.rt_ref(3)+delta;
end
if 0==His.rt_ref(5)
    t2 = His.rt_ref(1)-delta;
else
    t2 = His.rt_ref(5)-delta;
end
[rts10,top1_rt10] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
His.rt_ref(hno) = iffEmptyZero(rts10, top1_rt10);


% --- tiny utility for readability
function v = iffEmptyZero(rts, top1)
%% Return chosen RT or 0 when search failed
if isempty(rts), v = 0; else, v = top1; end
