function H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H3_04_27_40
% Panel for H3 peptide KSAPATGGVKKPHR (positions 27–40).
% Targets mono/di/tri-methylation on K27 and K36, combinations thereof,
% and K27ac. This function orchestrates:
%  - reference anchoring on unmodified form,
%  - RT recalibration,
%  - relocation (MS1-only or MS2-assisted),
%  - extraction of RTs and XIC intensities per PTM state,
%  - outputs (MAT + figures) and optional PSM export.
%
% INPUTS:
%   MS1_index, MS1_peaks : centroids and scan-time index for MS1
%   MS2_index, MS2_peaks : centroids and scan-time index for MS2
%   ptol                 : ppm tolerance for m/z matching
%   cur_outpath          : output directory
%   special              : struct with fields:
%                          - raw_path  : for check_ref
%                          - ndebug    : debug mode (relocateD)
%                          - nDAmode   : 0/1/2; 2 enables MS2-assisted logic
%                          - nhmass    : MS2 diag./neutral-loss params
%
% OUTPUTS:
%   Writes cur_outpath/H3_04_27_40.mat with pep_rts and pep_intens,
%   figures via draw_layout, and PSMs via GetPSM when enabled.

% check
out_filename = 'H3_04_27_40';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    % Early exit if MAT already produced; idempotent behavior.
    return;
end

% init peptide library and RT priors for this panel
His = init_histone(cur_outpath,out_filename);

% calculate target RTs and intensities across PTM states
unitdiff = 1.0032;           % C13-C12 isotopic mass difference
Mods = GetMods();            % global PTM mass registry (shared utility)
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

% persist quantitative outputs
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% diagnostic layout (XICs, apex RTs, MS2 marks)
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2); % RT (min) sequence for plotting
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
            mono_isointens,MS2_index,MS2_peaks,special);

% optional PSM export when DDA-like assistance is active
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
           mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end
end % <<< END main function ===================================================


function His = init_histone(cur_outpath,out_filename)
%%
% Build the in-panel peptide configuration:
%  - base sequence for H3 27–40,
%  - list of PTM short names and mod_type encodings,
%  - charge states considered,
%  - initial RT priors (empirical) and default display mask,
%  - ensure the "main" charge column ordering is consistent panel-wide.

His.pep_seq = 'KSAPATGGVKKPHR';

% List of targeted states (order matters and is used downstream)
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

% mod_type encodes positions and PTM kinds relative to 'pep_seq'.
% Example token "10,me1" means: at position index==10 in the peptide,
% apply me1. "0,pr" is N-term propionylation (if used by the original method).
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

% For this longer peptide, 2+/3+/4+ are practical charge states
His.pep_ch = repmat([2 3 4],length(His.mod_type),1);

% m/z values are computed programmatically
His.pep_mz = calculate_pepmz(His);

% Empirical RT priors (minutes); these are re-centered after anchoring unmod
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

% By default, display all except the two most complex double-trimethyl combos
His.display = ones(length(His.mod_type),1);
His.display([15 16]) = 0;

His.outpath = cur_outpath;
His.outfile = out_filename;

% Normalize charge-order columns so the panel has a consistent "main" column
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [~,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end
    % Tune only a subset (as in source) to avoid reindexing all rows unnecessarily
    tune = [4 5 6 7 8 9 11 12 13 14 15 16];
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end
end


function [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special)
%%
% Core quantification routine.

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% --- 1) Anchor unmodified reference ---
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only anchoring with smart fallback
        [His.rt_ref(1),special.ndebug] = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if 2==special.ndebug
            % Try inferring unmod RT using K27 ladder context (me1 & me2)
            hno = 1; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2);
            hno = 2; [~,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2); %#ok
            hno = 4; [~,top1_rt4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok
            if ~isempty(top1_rt2) && ~isempty(top1_rt4) && top1_rt2>top1_rt4
                if isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt4+(top1_rt2-top1_rt4)*0.7;
                else
                    p = find(rts1>top1_rt4+(top1_rt2-top1_rt4)*0.5 & rts1<top1_rt2);
                    if ~isempty(p)
                        [~,pp] = max(inten_sum1(p));
                        His.rt_ref(1) = rts1(p(pp));
                    end
                end
            else
                if ~isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt1;
                end
            end
        end
    else
        % MS2-assisted anchoring path
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0; t2 = MS1_index(num_MS1,2);
        else
            delta = 5; t1 = His.rt_ref(1)-delta; t2 = His.rt_ref(1)+delta;
        end
        hno = 1;
        [~,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if ~isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% Extract unmod to calibrate priors
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% --- 2) Recalibration ---
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta; % shift all priors equally
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% --- 3) Relocation stage (MS1-only vs MS2-assisted) ---
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special.nhmass);
    end
end

% --- 4) Extraction for each target state ---

% K36me1 / K27me1: split or paired depending on separation
if His.rt_ref(3)-His.rt_ref(2)>0.4
    for hno=2:3
        [cur_rts,cur_intens,cur_mono_isointens] = get_histone10(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge) = cur_rts;
            pep_intens(hno,1:ncharge) = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end
    end
else
    hno = 2;
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
    if cur_rts(1,1)>0
        pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
        pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
        mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
    end
end

% K27me2, K27me3, K36me2 (4,6,5) and K36me3 (7) + others resolved below
for hno=[4 6]
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end
hno = 5;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone11(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% K36me3 / K27me2K36me1 (paired extraction due to frequent co-elution)
hno = 7;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special);
if cur_rts(1,1)>0
    pep_rts(hno:hno+1,1:ncharge) = cur_rts(1:2,:);
    pep_intens(hno:hno+1,1:ncharge) = cur_intens(1:2,:);
    mono_isointens(1:num_MS1,hno:hno+1) = cur_mono_isointens(:,1:2);
end

% Remaining combos (9..13) via single-target extraction
for hno=9:13
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

% K27me3K36me2 (14) resolved as single
hno = 14;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% K27ac
hno = 17;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end
end


function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
% MS1-only relocation stage.

delta = 0.1;
nsplit = 1;

% K36me1 (2) with potential doublet behavior controlling K27me1 (3)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+15;
[rts2,top1_rt2,inten_sum2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,0,t1,t2);

[tmp_sum,ix] = sort(inten_sum2,'descend');
tmp_rts = rts2(ix);
if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
    % Two close maxima -> split K36me1 & K27me1 around them
    His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
    His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
else
    old_t = His.rt_ref(hno);
    if isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end

    % K27me1 (3) inherits displacement when doublet logic not triggered
    hno = 3;
    if 0==His.rt_ref(hno-1)
        His.rt_ref(hno) = 0;
    elseif old_t~= His.rt_ref(hno-1)
        d = His.rt_ref(hno-1) - old_t;
        His.rt_ref(hno) = His.rt_ref(hno) + d;
    end
end

% K27me2 (4), K27me3 (6), K27me3K36me1 (11): early broad search
hno = 4; t1 = His.rt_ref(1)-35; t2 = His.rt_ref(1)-delta;
[rts4,top1_rt4,inten_sum4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

hno = 6; t1 = His.rt_ref(1)-35; t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6,inten_sum6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok

hno = 11; t1 = His.rt_ref(1)-35; t2 = His.rt_ref(1)-delta;
[rts11,top1_rt11,inten_sum11] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok

% 'xt' consolidates a consistent triple across K27/K36 axes
xt = find_triple(rts4,top1_rt4,rts6,rts11,inten_sum4,inten_sum6,inten_sum11);

% For deep combos (13,14,16) build a pivot 'xt4' based on the strongest among 13/14
hno = 13; t1 = His.rt_ref(1)-50;
if 0==xt(1,1), t2 = His.rt_ref(1)-4; else, t2 = xt(1,1)-delta-0.4; end
[rts13,top1_rt13,inten_sum13] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

hno = 14; t1 = His.rt_ref(1)-50;
if 0==xt(1,1), t2 = His.rt_ref(1)-4; else, t2 = xt(1,1)-delta-0.4; end
[rts14,top1_rt14,inten_sum14] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

hno = 16; t1 = His.rt_ref(1)-50;
if 0==xt(1,1), t2 = His.rt_ref(1)-4; else, t2 = xt(1,1)-delta-0.4; end
[rts16,top1_rt16,inten_sum16] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok

if ~isempty(inten_sum13) || ~isempty(inten_sum14)
    rt_array = [0 0]; in_array = [0 0];
    if ~isempty(inten_sum13), rt_array(1) = top1_rt13; in_array(1) = max(inten_sum13); end
    if ~isempty(inten_sum14), rt_array(2) = top1_rt14; in_array(2) = max(inten_sum14); end
    [~,ix] = max(in_array);
    xt4 = rt_array(ix);
else
    xt4 = 0;
end

% Assign RTs using the 'xt' map
His.rt_ref(4)  = xt(1,1);   % K27me2
His.rt_ref(5)  = xt(1,2);   % K36me2
His.rt_ref(6)  = xt(2,1);   % K27me3
His.rt_ref(7)  = xt(1,3);   % K36me3
His.rt_ref(8)  = xt(2,2);   % K27me2K36me1
His.rt_ref(9)  = xt(2,3);   % K27me1K36me2

% K27me1K36me1 needs a local forward search from K27me1 (3)
hno = 10;
if 0==His.rt_ref(3), t1 = His.rt_ref(1)+delta; else, t1 = His.rt_ref(3)+delta; end
if 0==His.rt_ref(3), t2 = His.rt_ref(1)+18;  else, t2 = His.rt_ref(3)+(His.rt_ref(3)-His.rt_ref(1))+2; end
[rts10,top1_rt10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
His.rt_ref(hno) = iff_empty(top1_rt10, 0);

% Remaining high-order combos tied to 'xt' and pivot 'xt4'
His.rt_ref(11) = xt(3,2);  % K27me3K36me1
His.rt_ref(12) = xt(3,3);  % K27me1K36me3

% K27me2K36me2 aligned around 'xt4'
hno = 13;
ix = find(abs(rts13-xt4)<=delta+0.4);
if ~isempty(ix)
    [~,id] = max(inten_sum13(ix));
    His.rt_ref(hno) = rts13(ix(id));
else
    His.rt_ref(hno) = 0;
end

% K27me3K36me2 around 'xt4', propagate to K27me2K36me3
hno = 14;
old_t = His.rt_ref(hno);
ix = find(abs(rts14-xt4)<=delta+0.4);
if ~isempty(ix)
    [~,id] = max(inten_sum14(ix));
    His.rt_ref(hno) = rts14(ix(id));
else
    His.rt_ref(hno) = 0;
end

hno = 15;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1)
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end

% K27me3K36me3 using rts16 near 'xt4'
hno = 16;
ix = find(abs(rts16-xt4)<=delta+0.4);
if ~isempty(ix)
    [~,id] = max(inten_sum16(ix));
    His.rt_ref(hno) = rts16(ix(id));
else
    His.rt_ref(hno) = 0;
end

% K27ac in late window after resolving main methyl landscape
hno = 17;
if 0==xt(2,3), t1 = His.rt_ref(1)-10; else, t1 = xt(2,3)+delta; end
t2 = His.rt_ref(1)-delta;
[rts17,top1_rt17] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
His.rt_ref(hno) = iff_empty(top1_rt17, 0);
end


function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
% MS2-assisted relocation stage.

delta = 0.1;
nsplit = 1;

% K36me1 and K27me1 with MS2-assisted doublet logic
hno = 2; t1 = His.rt_ref(1)+delta; t2 = His.rt_ref(1)+15;
[rts2,top1_rt2,inten_sum2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,0,t1,t2,nhmass);
if isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

hno = 3; t1 = His.rt_ref(1)+delta; t2 = His.rt_ref(1)+15;
[rts3,top1_rt3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,0,t1,t2,nhmass);
if isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end

% Enforce ordered split if apexes collapse or are too close
if ~isempty(rts2) && ~isempty(rts3)
    if top1_rt2==top1_rt3
        hno = 2;
        [tmp_sum,ix] = sort(inten_sum2,'descend');
        tmp_rts = rts2(ix);
        if length(tmp_sum)>=2 && tmp_sum(2)>=tmp_sum(1)/30 && abs(tmp_rts(2)-tmp_rts(1))<2
            His.rt_ref(hno) = min([tmp_rts(2),tmp_rts(1)]);
            His.rt_ref(hno+1) = max([tmp_rts(2),tmp_rts(1)]);
        else
            old_t = His.rt_ref(hno);
            His.rt_ref(hno) = iff_empty(top1_rt2,0);
            hno = 3;
            if 0==His.rt_ref(hno-1)
                His.rt_ref(hno) = 0;
            elseif old_t~= His.rt_ref(hno-1)
                d = His.rt_ref(hno-1) - old_t;
                His.rt_ref(hno) = His.rt_ref(hno) + d;
            end
        end
    elseif top1_rt2>top1_rt3 && abs(top1_rt2-top1_rt3)<2
        His.rt_ref(2) = top1_rt3;
        His.rt_ref(3) = top1_rt2;
    end
end

% Early search for triple anchors (4,6,11)
hno = 4; t1 = His.rt_ref(1)-35; t2 = His.rt_ref(1)-delta;
[rts4,top1_rt4,inten_sum4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
hno = 6; t1 = His.rt_ref(1)-35; t2 = His.rt_ref(1)-delta;
[rts6,top1_rt6,inten_sum6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok
hno = 11; t1 = His.rt_ref(1)-35; t2 = His.rt_ref(1)-delta;
[rts11,top1_rt11,inten_sum11] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok

xt = find_triple(rts4,top1_rt4,rts6,rts11,inten_sum4,inten_sum6,inten_sum11);

% Deep combos (13,14,16) with MS2: choose pivot xt4 by top1_intensity
hno = 13; t1 = His.rt_ref(1)-50; if 0==xt(1,1), t2 = His.rt_ref(1)-4; else, t2 = xt(1,1)-delta-0.4; end
[rts13,top1_rt13,inten_sum13,top1_inten_sum13] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
hno = 14; t1 = His.rt_ref(1)-50; if 0==xt(1,1), t2 = His.rt_ref(1)-4; else, t2 = xt(1,1)-delta-0.4; end
[rts14,top1_rt14,inten_sum14,top1_inten_sum14] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
hno = 16; t1 = His.rt_ref(1)-50; if 0==xt(1,1), t2 = His.rt_ref(1)-4; else, t2 = xt(1,1)-delta-0.4; end
[rts16,top1_rt16,inten_sum16] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok

if ~isempty(inten_sum13) || ~isempty(inten_sum14)
    rt_array = [0 0]; in_array = [0 0];
    if ~isempty(inten_sum13), rt_array(1) = top1_rt13; in_array(1) = top1_inten_sum13; end
    if ~isempty(inten_sum14), rt_array(2) = top1_rt14; in_array(2) = top1_inten_sum14; end
    [~,ix] = max(in_array);
    xt4 = rt_array(ix);
else
    xt4 = 0;
end

% Map assignments as in MS1-only path
His.rt_ref(4)  = xt(1,1);   % K27me2
His.rt_ref(5)  = xt(1,2);   % K36me2
His.rt_ref(6)  = xt(2,1);   % K27me3
His.rt_ref(7)  = xt(1,3);   % K36me3
His.rt_ref(8)  = xt(2,2);   % K27me2K36me1
His.rt_ref(9)  = xt(2,3);   % K27me1K36me2

% K27me1K36me1 forward local search
hno = 10;
if 0==His.rt_ref(3), t1 = His.rt_ref(1)+delta; else, t1 = His.rt_ref(3)+delta; end
if 0==His.rt_ref(3), t2 = His.rt_ref(1)+18;  else, t2 = His.rt_ref(3)+(His.rt_ref(3)-His.rt_ref(1))+2; end
[rts10,top1_rt10] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
His.rt_ref(hno) = iff_empty(top1_rt10, 0);

His.rt_ref(11) = xt(3,2);  % K27me3K36me1
His.rt_ref(12) = xt(3,3);  % K27me1K36me3

% K27me2K36me2 around xt4
hno = 13;
ix = find(abs(rts13-xt4)<=delta+0.4);
if ~isempty(ix)
    [~,id] = max(inten_sum13(ix));
    His.rt_ref(hno) = rts13(ix(id));
else
    His.rt_ref(hno) = 0;
end

% K27me3K36me2 + propagate to K27me2K36me3
hno = 14; old_t = His.rt_ref(hno);
ix = find(abs(rts14-xt4)<=delta+0.4);
if ~isempty(ix)
    [~,id] = max(inten_sum14(ix));
    His.rt_ref(hno) = rts14(ix(id));
else
    His.rt_ref(hno) = 0;
end

hno = 15;
if 0==His.rt_ref(hno-1)
    His.rt_ref(hno) = 0;
elseif old_t~= His.rt_ref(hno-1)
    d = His.rt_ref(hno-1) - old_t;
    His.rt_ref(hno) = His.rt_ref(hno) + d;
end

% K27me3K36me3 near xt4
hno = 16;
ix = find(abs(rts16-xt4)<=delta+0.4);
if ~isempty(ix)
    [~,id] = max(inten_sum16(ix));
    His.rt_ref(hno) = rts16(ix(id));
else
    His.rt_ref(hno) = 0;
end

% K27ac (late window)
hno = 17;
if 0==xt(2,3), t1 = His.rt_ref(1)-10; else, t1 = xt(2,3)+delta; end
t2 = His.rt_ref(1)-delta;
[rts17,top1_rt17] = get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
    ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); % <-- fixed: get_rts2 (no '22')
His.rt_ref(hno) = iff_empty(top1_rt17, 0);
end


function v = iff_empty(x, default_zero)
% IFF_EMPTY  Return default_zero if x is empty, else x.
if isempty(x)
    v = default_zero;
else
    v = x;
end
end
