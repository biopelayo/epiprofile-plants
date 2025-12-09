function H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%% H3_03_18_26 — Targeted panel for H3 KQLATKAAR (18–26)
%
% HIGH-LEVEL PURPOSE
%   Quantify the histone H3 peptide "KQLATKAAR" (positions 18–26) and
%   its PTM states at K18, T22, and K23 (me1/ac/pr combinations).
%   Workflow:
%     1) Initialize panel: sequences, PTMs, charge states, RT priors.
%     2) Anchor the unmodified form; re-center all RTs accordingly.
%     3) Refine RTs via relocation (MS1-only or MS2-assisted).
%     4) Extract per-PTM retention times (RTs) and intensities (XIC areas).
%     5) Save outputs, draw diagnostic layout, optionally export PSMs.
%
% KEY PEPTIDE POSITIONS (1-based within "KQLATKAAR"):
%   K18→index 1,  T22→index 5,  K23→index 6. "0" denotes peptide N-terminus.
%
% INPUTS
%   MS1_index, MS1_peaks : MS1 containers (rt grid + centroided peaks)
%   MS2_index, MS2_peaks : MS2 containers (scans + centroided peaks)
%   ptol                 : ppm tolerance (100 is internally handled as 10)
%   cur_outpath          : output folder for this RAW
%   special              : struct with fields nDAmode, ndebug, nhmass, raw_path
%
% SIDE-EFFECTS
%   Writes <out>.mat, figures via draw_layout(), and PSMs if enabled.

%% --- 0) Idempotent guard: skip if .mat already exists
out_filename = 'H3_03_18_26';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end

%% --- 1) Initialize panel definition
His = init_histone(cur_outpath,out_filename);

%% --- 2) Core computation: anchor → relocate → extract
unitdiff = 1.0032;      % ~13C–12C mass difference
Mods     = GetMods();   % global PTM mass table (used by paired extractors)
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

%% --- 3) Persist per-PTM quantifications
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

%% --- 4) Plot diagnostic layout (XICs / markers)
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

%% --- 5) Optional: export PSMs if DA/DIA mode is enabled
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end


% =========================================================================
%                               HELPERS
% =========================================================================
function His = init_histone(cur_outpath,out_filename)
%% init_histone — Define PTM panel for KQLATKAAR (18–26)
%
% PANEL CONTENT
%   - Unmodified backbone (propionyl tags per framework conventions).
%   - K18: me1 / ac; K23: me1 / ac; double me1 at K18+K23.
%   - T22: ac / pr, and a mixed K18ac+T22pr.
%   - Paired extractors are used where co-elution is plausible (K18/K23).
%
His.pep_seq = 'KQLATKAAR';
His.mod_short = {'unmod';
    'K23me1';
    'K18me1';
    'K18me1K23me1';
    'K18ac';
    'K23ac';
    'K18acK23ac';
    'K18acT22pr';
    'T22ac';
    'T22pr'};

% mod_type uses "pos,mod;" pairs:
%   0 = peptide N-terminus (framework propionylation tag)
%   1 = K18, 5 = T22, 6 = K23
His.mod_type = {'0,pr;1,pr;6,pr;';
    '0,pr;1,pr;6,me1;';
    '0,pr;1,me1;6,pr;';
    '0,pr;1,me1;6,me1;';
    '0,pr;1,ac;6,pr;';
    '0,pr;1,pr;6,ac;';
    '0,pr;1,ac;6,ac;';
    '0,pr;1,ac;5,pr;6,pr;';
    '0,pr;1,pr;5,ac;6,pr;';
    '0,pr;1,pr;5,pr;6,pr;'};

% Consider 1+, 2+, 3+ charge states for all entries
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);

% Theoretical m/z for each entry×charge (computed from Mods table)
% (keeping codebase function to avoid hard-coded numbers)
His.pep_mz = calculate_pepmz(His);

% RT priors (minutes) — used as initial guesses before re-centering
His.rt_ref = [37.41
    39.76
    39.77
    40.96
    36.13
    36.14
    34.64
    40.04
    40.05
    41.27];

% Display toggles in the layout; hide T22-focused entries by default
His.display = ones(length(His.mod_type),1);
His.display([8 9 10]) = 0;

His.outpath = cur_outpath;
His.outfile = out_filename;

% Normalize the charge-column order to keep a consistent "main" charge
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
%% calculate_layout — Anchor, relocate, and extract the K18/T22/K23 panel
%
% STRATEGY
%   1) Anchor using "unmod" (entry #1) with check_ref() + get_rts()/get_rts2().
%   2) Re-center all RT priors based on the observed unmod RT.
%   3) Relocate with MS1-only or MS2-assisted rules (see relocate*()).
%   4) Extract per-PTM RT/intensity (XICs); use paired extractors when
%      co-elution is expected (e.g., K18 vs K23 variants).
%
[npep,ncharge] = size(His.pep_mz);
num_MS1        = size(MS1_index,1);
pep_rts        = zeros([npep,ncharge]);
pep_intens     = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% --- Anchor on "unmod" (entry #1). If ndebug toggles, broaden search.
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        [His.rt_ref(1),special.ndebug] = check_ref( ...
            special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if 2==special.ndebug
            % Fallback: triangulate with K23me1 (hno=2) and K18ac (hno=5)
            hno = 1; t1=0; t2=MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);
            hno = 2; [rts2,top1_rt2]   = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2); %#ok<ASGLU>
            hno = 5; [rts5,top1_rt5]   = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok<ASGLU>
            if 0==isempty(top1_rt2) && 0==isempty(top1_rt5) && top1_rt2>top1_rt5
                if 1==isempty(top1_rt1)
                    % Use a weighted point between K18ac and K23me1 as anchor
                    His.rt_ref(1) = top1_rt5+(top1_rt2-top1_rt5)*0.4;
                else
                    % Snap to the strongest unmod peak between those bounds
                    p = find(rts1>top1_rt5+(top1_rt2-top1_rt5)*0.2 & rts1<top1_rt2);
                    if 0==isempty(p)
                        [~,pp] = max(inten_sum1(p));
                        His.rt_ref(1) = rts1(p(pp));
                    end
                end
            else
                if 0==isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt1;
                end
            end
        end
    else
        % MS2-assisted anchoring
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1=0; t2=MS1_index(num_MS1,2);
        else
            delta=5; t1=His.rt_ref(1)-delta; t2=His.rt_ref(1)+delta;
        end
        hno = 1;
        [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                                   ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok<ASGLU>
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% --- Extract unmod trace; re-center all priors to observed unmod RT
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);
if cur_rts(1)>0
    His.rt_ref(1)     = cur_rts(1);
    delta             = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end)+delta;
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% --- Relocation rules (MS1-only vs MS2-assisted)
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end
end

% --- K23me1 / K18me1 (pair or separated)
% If their expected RTs are >0.4 min apart, extract individually;
% otherwise, use a 2-target extractor to resolve partial co-elution.
if His.rt_ref(3)-His.rt_ref(2)>0.4
    for hno=2:3
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone10(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge)        = cur_rts;
            pep_intens(hno,1:ncharge)     = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end
    end
else
    hno = 2;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
        MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge)        = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge)     = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end
end

% --- K18me1K23me1 (double me1)
hno = 4;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% --- K18ac / K23ac (paired, often partially co-eluting)
hno = 5;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge)        = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge)     = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end

% --- K18acK23ac (di-acetyl)
hno = 7;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% --- T22 branch: T22ac vs K18acT22pr (pair or separated)
if His.rt_ref(9)-His.rt_ref(8)>0.4
    for hno=8:9
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone10(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge)        = cur_rts;
            pep_intens(hno,1:ncharge)     = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end
    end
else
    hno = 8;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2( ...
        MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge)        = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge)     = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end
end

% --- T22pr (unpaired)
hno = 10;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end


function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%% relocate — MS1-only relocation rules around the unmod anchor
%
% WINDOW HEURISTICS (minutes)
%   K23me1:      [anchor+0.1, anchor+14], with potential 2-peak tie handling
%   K18me1:      tracks K23me1 shift when co-eluting
%   K18me1K23me1:[anchor+Δ, anchor+~22] (adaptive to K18me1)
%   K18ac:       [anchor-14, anchor-0.1], K23ac tracks K18ac
%   K18acK23ac:  [anchor-20, min(K18ac,anchor)-0.1]
%   T22ac:       [anchor+0.1, anchor+14], T22pr tracks T22ac
%   T22pr:       [max(anchor, T22ac)+0.1, anchor+20]
%
delta  = 0.1;
nsplit = 1;

% -- K23me1 (and infer K18me1 shift if two close maxima are found)
hno = 2;
t1 = His.rt_ref(1)+delta; t2 = His.rt_ref(1)+14;
[rts2,top1_rt2,inten_sum2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2);
[tmp_sum,ix] = sort(inten_sum2,'descend');
tmp_rts = rts2(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
    His.rt_ref(hno)   = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts2), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt2; end
    % K18me1 follows the K23me1 shift when applicable
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1)
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end
end

% -- K18me1K23me1 (double me1)
hno = 4;
if 0==His.rt_ref(3), t1 = His.rt_ref(1)+delta; else, t1 = His.rt_ref(3)+delta; end
if 0==His.rt_ref(3), t2 = His.rt_ref(1)+22;    else, t2 = His.rt_ref(3)+(His.rt_ref(3)-His.rt_ref(1))+2; end
[rts4,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts4), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt4; end

% -- K18ac
hno = 5;
t1 = His.rt_ref(1)-14; t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
old_t = His.rt_ref(hno);
if 1==isempty(rts5), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt5; end

% -- K23ac tracks K18ac movement
hno = 6;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1)
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end

% -- K18acK23ac
hno = 7;
t1 = His.rt_ref(1)-20;
if 0==His.rt_ref(5), t2 = His.rt_ref(1)-delta; else, t2 = His.rt_ref(5)-delta; end
[rts7,top1_rt7] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts7), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt7; end

% -- T22ac (and infer K18acT22pr shift if two close maxima)
hno = 8;
t1 = His.rt_ref(1)+delta; t2 = His.rt_ref(1)+14;
[rts8,top1_rt8,inten_sum8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
[tmp_sum,ix] = sort(inten_sum8,'descend');
tmp_rts = rts8(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
    His.rt_ref(hno)   = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts8), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt8; end
    % K18acT22pr follows T22ac shift
    hno = 9;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1)
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end
end

% -- T22pr
hno = 10;
if 0==His.rt_ref(8), t1 = His.rt_ref(1)+delta; else, t1 = His.rt_ref(8)+delta; end
t2 = His.rt_ref(1)+20;
[rts10,top1_rt10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts10), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt10; end


function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%% relocate2 — MS2-assisted relocation rules (diagnostic fragments)
%
delta  = 0.1;
nsplit = 1;

% -- K23me1
hno = 2;
t1 = His.rt_ref(1)+delta; t2 = His.rt_ref(1)+14;
[rts2,top1_rt2,inten_sum2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,0,t1,t2,nhmass);
if 1==isempty(rts2), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt2; end

% -- K18me1
hno = 3;
t1 = His.rt_ref(1)+delta; t2 = His.rt_ref(1)+14;
[rts3,top1_rt3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,0,t1,t2,nhmass);
if 1==isempty(rts3), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt3; end

% If both present and equal/near-equal, resolve ordering or tie by intensity
if 0==isempty(rts2) && 0==isempty(rts3)
    if top1_rt2==top1_rt3
        hno = 2;
        [tmp_sum,ix] = sort(inten_sum2,'descend');
        tmp_rts = rts2(ix);
        if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
            His.rt_ref(hno)   = min([tmp_rts(2),tmp_rts(1)]);
            His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
        else
            old_t = His.rt_ref(hno);
            if 1==isempty(rts2), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt2; end
            % K18me1 follows K23me1 shift
            hno = 3;
            if 0==His.rt_ref(hno-1)
                His.rt_ref(hno) = 0;
            elseif old_t~= His.rt_ref(hno-1)
                d = His.rt_ref(hno-1) - old_t;
                His.rt_ref(hno) = His.rt_ref(hno) + d;
            end
        end
    elseif top1_rt2>top1_rt3 && abs(top1_rt2-top1_rt3)<2
        % Small RT gap: enforce ordering (K18me1 first, then K23me1)
        hno = 2;
        His.rt_ref(hno)   = top1_rt3;
        His.rt_ref(hno+1) = top1_rt2;
    end
end

% -- K18me1K23me1
hno = 4;
if 0==His.rt_ref(3), t1 = His.rt_ref(1)+delta; else, t1 = His.rt_ref(3)+delta; end
if 0==His.rt_ref(3), t2 = His.rt_ref(1)+22;    else, t2 = His.rt_ref(3)+(His.rt_ref(3)-His.rt_ref(1))+2; end
[rts4,top1_rt4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts4), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt4; end

% -- K23ac (MS2-based); K18ac will track K23ac motion below
hno = 6;
t1 = His.rt_ref(1)-14; t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
old_t = His.rt_ref(hno);
if 1==isempty(rts6), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt6; end

% -- K18ac tracks K23ac shift
hno = 5;
if 0==His.rt_ref(hno+1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno+1)
    d = His.rt_ref(hno+1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end

% -- K18acK23ac
hno = 7;
t1 = His.rt_ref(1)-20;
if 0==His.rt_ref(5), t2 = His.rt_ref(1)-delta; else, t2 = His.rt_ref(5)-delta; end
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts7), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt7; end

% -- T22ac (MS2-based) with possible two-peak resolution
hno = 8;
t1 = His.rt_ref(1)+delta; t2 = His.rt_ref(1)+14;
[rts8,top1_rt8,inten_sum8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
[tmp_sum,ix] = sort(inten_sum8,'descend');
tmp_rts = rts8(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
    His.rt_ref(hno)   = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if 1==isempty(rts8), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt8; end
    % K18acT22pr follows T22ac shift
    hno = 9;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1)
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end
end

% -- T22pr
hno = 10;
if 0==His.rt_ref(8), t1 = His.rt_ref(1)+delta; else, t1 = His.rt_ref(8)+delta; end
t2 = His.rt_ref(1)+20;
[rts10,top1_rt10] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts10), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt10; end
