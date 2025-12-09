function H3_02a_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%% H3_02a_9_17 — Targeted layout for H3 peptide KSTGGKAPR with S10 phosphorylation
%
% HIGH-LEVEL PURPOSE
%   Quantify the histone H3 peptide "KSTGGKAPR" (positions 9–17 region)
%   but focusing on phosphorylation at S10 (S10ph), optionally combined
%   with K9 methyl/acetyl states and K14 acetylation. This extends the
%   non-phospho panel in H3_02_9_17 to the phospho-centric panel.
%
% KEY DIFFERENCE VS H3_02_9_17
%   - Every PTM state here contains S10ph (Ser at peptide position 2),
%     except that "S10ph" alone acts as the unmodified anchor for this
%     panel (i.e., anchor == phospho baseline).
%   - MS2 behavior changes: for phosphate-containing precursors, consider
%     neutral loss in CID/HCD (98 Da), which may shift the most intense
%     fragments and improves MS2-assisted relocation; helpers (get_rts2)
%     should be robust to this via diagnostic-ion logic in your codebase.
%
% INPUTS
%   MS1_index, MS1_peaks : MS1 containers (see global conventions)
%   MS2_index, MS2_peaks : MS2 containers (see global conventions)
%   ptol      : ppm tolerance (100 is internally treated as 10 in helpers)
%   cur_outpath : output folder
%   special   : struct with run options (nDAmode, ndebug, nhmass, raw_path)
%
% OUTPUTS (files)
%   <out>.mat, <out>.plabel (optional), figures via draw_layout(...)
%
% DEPENDENCIES
%   GetMods, calculate_pepmz, check_ref, get_rts/get_rts2/get_rts22,
%   get_histone0/1/2/11, find_pair_new, output_histone, draw_layout, GetPSM
%

% --- 0) Idempotent skip if results already exist
out_filename = 'H3_02a_9_17';
% fprintf(1,'%s..',out_filename); % kept commented as in original
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end

% --- 1) Initialize phospho-centric histone panel (masses, charges, RT)
His = init_histone(cur_outpath,out_filename);

% --- 2) Compute per-PTM RT/intensity (MS1-only or MS2-assisted)
unitdiff = 1.0032;      % 13C–12C mass difference
Mods     = GetMods();
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

% --- 3) Save main quantitative outputs
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% --- 4) Plot diagnostic layout
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% --- 5) Optional PSM export
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end


% =========================================================================
%                           NESTED HELPERS
% =========================================================================
function His = init_histone(cur_outpath,out_filename)
%% init_histone — Define KSTGGKAPR PTM panel centered on S10ph
%
% WHAT IS NEW HERE
%   - Baseline “anchor” is S10ph (not fully unmodified peptide).
%   - Combine S10ph with K9 me1/me2/me3 or ac, and with K14ac.
%   - Charge set replicated for 1+, 2+, 3+; theoretical m/z computed.
%
His.pep_seq = 'KSTGGKAPR';
His.mod_short = {'S10ph';
    'K9me1S10ph';
    'K9me2S10ph';
    'K9me3S10ph';
    'K9acS10ph';
    'S10phK14ac';
    'K9me1S10phK14ac';
    'K9me2S10phK14ac';
    'K9me3S10phK14ac';
    'K9acS10phK14ac'};

% 'mod_type' uses the internal encoding "pos,mod;" where positions are
% indexed over the peptide string "K S T G G K A P R" => S10 is position 2,
% K9 is position 1, K14 is position 6 in this peptide context.
His.mod_type = {'0,pr;1,pr;2,ph;6,pr;';
    '0,pr;1,me1;2,ph;6,pr;';
    '0,pr;1,me2;2,ph;6,pr;';
    '0,pr;1,me3;2,ph;6,pr;';
    '0,pr;1,ac;2,ph;6,pr;';
    '0,pr;1,pr;2,ph;6,ac;';
    '0,pr;1,me1;2,ph;6,ac;';
    '0,pr;1,me2;2,ph;6,ac;';
    '0,pr;1,me3;2,ph;6,ac;';
    '0,pr;1,ac;2,ph;6,ac;'};

% Charges considered
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);

% Theoretical m/z per PTM×charge (computed with Mods mass table)
% (commented table kept in original source is omitted to avoid confusion)
His.pep_mz = calculate_pepmz(His);

% RT priors (minutes). These are initial anchors (empirical/prior);
% they will be re-centered based on the observed “anchor” (S10ph) RT.
His.rt_ref = [23.20
    25.91
    16.59
    16.46
    21.90
    21.91
    24.68
    14.95
    14.86
    20.23];
His.display = ones(length(His.mod_type),1);

His.outpath  = cur_outpath;
His.outfile  = out_filename;

% Align charge ordering to a main column (use row1 col2 as reference)
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
%% calculate_layout — Anchor on S10ph, relocate phospho-combinatorials
%
% STRATEGY
%   1) Use S10ph as the panel anchor (His.rt_ref(1)).
%   2) Extract S10ph with get_histone0 and re-center all downstream RT priors.
%   3) Relocate modified forms (me/acetyl + phospho) via MS1-only or MS2-assisted.
%   4) Extract final RT/intensity for each PTM state.
%
[npep,ncharge]    = size(His.pep_mz);
num_MS1           = size(MS1_index,1);
pep_rts           = zeros([npep,ncharge]);
pep_intens        = zeros([npep,ncharge]);
mono_isointens    = zeros([num_MS1,npep]);

% --- 1) Anchor on S10ph (serves as "unmod" within this phospho panel)
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only anchoring with optional reference check against S10ph spec
        [His.rt_ref(1),special.ndebug] = check_ref( ...
            special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if 2==special.ndebug
            % Wide search for anchor (here hno=1 → S10ph)
            hno = 1; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2);
            % Compare to K9me1S10ph to disambiguate if needed
            hno = 2; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok<ASGLU>
            % Selection heuristic: prefer earlier of two strong maxima if close
            [tmp_sum,ix] = sort(inten_sum1,'descend'); %#ok<NASGU>
            tmp_rts = rts1(ix);
            if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<7
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
        % MS2-assisted anchoring (benefits from phosphate NL behavior)
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

% --- 2) Extract anchor quantitative trace and calibrate RT priors
hno = 1; % S10ph
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

% --- 3) Relocation of phospho-combinatorials
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end
end

% --- 4) Final extraction per PTM state

% K9me1/2/3 + S10ph
for hno=2:4
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)        = cur_rts;
        pep_intens(hno,1:ncharge)     = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

% (K9acS10ph, S10phK14ac) via the specialized extractor
hno = 5;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge)           = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge)        = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1)    = cur_mono_isointens(:,1:2);
end

% K9me1/2 S10ph K14ac and K9acS10phK14ac
for hno=[7 8 10]
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)        = cur_rts;
        pep_intens(hno,1:ncharge)     = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

% K9me3S10phK14ac (may have slightly different XIC behavior → get_histone11)
hno = 9;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end


function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%% relocate — MS1-only relocation rules (phospho anchor)
%
% WINDOWS (heuristic; centered around S10ph anchor)
%   - K9me1S10ph after anchor:   [anchor+0.5, anchor+16]
%   - K9me2/3S10ph before anchor: [6, anchor-3], pairing me2 vs me3
%   - K9acS10ph between me2S10ph and anchor; S10phK14ac follows K9acS10ph
%   - Combinatorials constrained between their single-PTM bounds
%
delta  = 0.5;
nsplit = 1;

% -- K9me1S10ph
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts2), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt2; end

% -- K9me2S10ph
hno = 3;
t1 = 6; t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% -- K9me3S10ph
hno = 4;
t1 = 6; t2 = His.rt_ref(1)-3;
[rts4,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok<ASGLU>
[His.rt_ref(4),His.rt_ref(3)] = find_pair_new(top1_rt4, rts3, top1_rt3, inten_sum3, 0);

% -- K9acS10ph (between me2S10ph and anchor)
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
old_t = His.rt_ref(hno);
if 1==isempty(rts5), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt5; end

% -- S10phK14ac tracks K9acS10ph if K9acS10ph moved
hno = 6;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1)
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end

% -- K9me1S10phK14ac
hno = 7;
t1 = His.rt_ref(1)+delta;
if 0==His.rt_ref(2), t2 = His.rt_ref(1)+9; else, t2 = His.rt_ref(2)-delta; end
[rts7,top1_rt7] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts7), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt7; end

% -- K9me2S10phK14ac
hno = 8;
t1 = 6;
if 0==His.rt_ref(3), t2 = His.rt_ref(1)-4; else, t2 = His.rt_ref(3)-delta; end
[rts8,top1_rt8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok<ASGLU>

% -- K9me3S10phK14ac
hno = 9;
t1 = 6;
if 0==His.rt_ref(3), t2 = His.rt_ref(1)-4; else, t2 = His.rt_ref(3)-delta; end
[rts9,top1_rt9,inten_sum9] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
[His.rt_ref(8),His.rt_ref(9)] = find_pair_new(top1_rt8, rts9, top1_rt9, inten_sum9, 1);

% -- K9acS10phK14ac
hno = 10;
if 0==His.rt_ref(3), t1 = His.rt_ref(1)-18; else, t1 = His.rt_ref(3)+delta; end
if 0==His.rt_ref(5), t2 = His.rt_ref(1)-delta; else, t2 = His.rt_ref(5)-delta; end
[rts10,top1_rt10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts10), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt10; end


function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%% relocate2 — MS2-assisted relocation (benefits from phosphate NL behavior)
%
delta  = 0.5;
nsplit = 1;

% -- K9me1S10ph
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts2), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt2; end

% -- K9me2S10ph
hno = 3;
t1 = 6; t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% -- K9me3S10ph
hno = 4;
t1 = 6; t2 = His.rt_ref(1)-3;
[rts4,top1_rt4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok<ASGLU>
[His.rt_ref(4),His.rt_ref(3)] = find_pair_new(top1_rt4, rts3, top1_rt3, inten_sum3, 0);

% -- S10phK14ac (between me2S10ph and anchor)
hno = 6;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
old_t = His.rt_ref(hno);
if 1==isempty(rts6), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt6; end

% -- K9acS10ph tracks S10phK14ac if that moved
hno = 5;
if 0==His.rt_ref(hno+1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno+1)
    d = His.rt_ref(hno+1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end

% -- K9me1S10phK14ac
hno = 7;
t1 = His.rt_ref(1)+delta;
if 0==His.rt_ref(2), t2 = His.rt_ref(1)+9; else, t2 = His.rt_ref(2)-delta; end
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts7), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt7; end

% -- K9me2S10phK14ac
hno = 8;
t1 = 6;
if 0==His.rt_ref(3), t2 = His.rt_ref(1)-4; else, t2 = His.rt_ref(3)-delta; end
[rts8,top1_rt8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok<ASGLU>

% -- K9me3S10phK14ac
hno = 9;
t1 = 6;
if 0==His.rt_ref(3), t2 = His.rt_ref(1)-4; else, t2 = His.rt_ref(3)-delta; end
[rts9,top1_rt9,inten_sum9] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
[His.rt_ref(8),His.rt_ref(9)] = find_pair_new(top1_rt8, rts9, top1_rt9, inten_sum9, 1);

% -- K9acS10phK14ac
hno = 10;
if 0==His.rt_ref(3), t1 = His.rt_ref(1)-18; else, t1 = His.rt_ref(3)+delta; end
if 0==His.rt_ref(5), t2 = His.rt_ref(1)-delta; else, t2 = His.rt_ref(5)-delta; end
[rts10,top1_rt10] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts10), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt10; end
