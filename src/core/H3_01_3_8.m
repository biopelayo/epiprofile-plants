function H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%% H3_01_3_8 — Targeted layout for H3 N-term peptide TKQTAR with K4 PTMs (unmod, me1, me2, me3, ac)
%
% HIGH-LEVEL PURPOSE (for non-proteomics readers)
%   This routine quantifies a specific histone H3 N-terminal tryptic peptide
%   "TKQTAR" (positions 1–6; target residue K4) across five modification
%   states: unmodified, K4me1, K4me2, K4me3, K4ac. It:
%     1) Initializes the peptide/charge/mass library for these PTMs.
%     2) Estimates peptide retention times (RTs) and intensities using MS1,
%        optionally assisted by MS2 (depending on 'special.nDAmode').
%     3) Aligns/calibrates RTs using the unmodified form as reference.
%     4) Writes quantitative outputs and diagnostic plots.
%     5) (Optional) Exports PSM-style labels for downstream spectrum viewers.
%
% INPUTS
%   MS1_index : [nMS1×4] table — columns: (scan_no, RT[min], cum_peak_start_idx, baseline)
%   MS1_peaks : [∑MS1peaks×2] — concatenated (m/z, intensity) for all MS1 scans
%   MS2_index : [nMS2×8] table — columns: (MS1_scan, MS1_RT, MS2_scan, m/z, z, FragType, MS2_peak_start_idx, baseline)
%   MS2_peaks : [∑MS2peaks×2] — concatenated (m/z, intensity) for all MS2 scans
%   ptol      : ppm tolerance (NOTE: if set to 100, the code normalizes it to 10 ppm internally)
%   cur_outpath : output folder for results (MAT files, plots, plabel)
%   special   : struct with runtime switches/metadata:
%                 • special.nDAmode : 0/1/2  (0/1 = MS1-only; 2 = MS2-assisted RT search)
%                 • special.ndebug  : flag for reference RT checking logic (see below)
%                 • special.nhmass  : if true, use alternative key-ion map in MS2 match
%                 • special.raw_path: raw file path used by check_ref(...)
%
% DEPENDENCIES (called by this routine)
%   get_rts, get_rts2, get_rts22, GetProfiles, GetTopBottom(11), GetLocal,
%   GetMods, calculate_pepmz, get_histone0, get_histone1, find_pair,
%   draw_layout, output_histone, GetPSM
%   (These are part of the EpiProfile/EpiProfile_PLANTS codebase.)
%
% OUTPUTS
%   Files are written under 'cur_outpath':
%     • <out_filename>.mat        — per-peptide quantitative summaries
%     • <out_filename>.plabel     — simple spectrum-label file for viewers
%     • figures / layout plots    — via draw_layout(...)
%   In-memory variables are local; persistency is through output files.
%
% NOTES
%   • 'unitdiff' is the 13C–12C mass difference per carbon (~1.0032 Da),
%     used to form isotopic ladders and precursor checks.
%   • The function uses the RT of unmodified TKQTAR as an anchor, then
%     relocates the expected windows for modified forms.
%   • When 'special.nDAmode==2', RT picking prefers MS2-supported windows
%     via get_rts2/get_rts22 (i.e., MS1–MS2 consensus).
%

% --- 0) Idempotency / skip if already processed for this out_filename
out_filename = 'H3_01_3_8';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    % If the MAT already exists, silently return (design choice for pipelines).
    return;
end

% --- 1) Initialize peptide/PTM library for H3 TKQTAR K4*
His = init_histone();

% --- 2) Compute per-peptide RTs/intensities using MS1, optionally with MS2 aid
unitdiff = 1.0032; % 13C–12C mass difference (Da), used for isotope spacing
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special);

% --- 3) Write primary quant result structures
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% --- 4) Plot: XICs, RT markers, etc., for diagnostics
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2); % RT vector per MS1 scan
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% --- 5) Optional: export simple PSM labels for viewer tools
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end


% ========================================================================
%                         NESTED HELPERS
% ========================================================================

function His = init_histone()
%% init_histone — Define the peptide sequence, PTMs, charge states and RT priors
%
% WHAT IT SETS
%   • His.pep_seq     : the bare sequence "TKQTAR"
%   • His.mod_short   : human-friendly short labels
%   • His.mod_type    : internal "position,mod;" list (0,pr = peptide N-term propionylation)
%   • His.pep_ch      : charge states considered per PTM (columns are charges)
%   • His.pep_mz      : theoretical m/z per PTM×charge (from calculate_pepmz)
%   • His.rt_ref      : initial RT priors (minutes) used to constrain searches
%   • His.display     : display flags (1 = plotted)
%   • Charge ordering : ensure all PTMs share a consistent "main" charge column
%
His.pep_seq = 'TKQTAR';
His.mod_short = { ...
    'unmod'; ...
    'K4me1'; ...
    'K4me2'; ...
    'K4me3'; ...
    'K4ac'  };
His.mod_type = { ...
    '0,pr;2,pr;';  ... % unmodified peptide N-term + Lys propionylation
    '0,pr;2,me1;'; ...
    '0,pr;2,me2;'; ...
    '0,pr;2,me3;'; ...
    '0,pr;2,ac;'   };

% Two charge states per PTM (col1, col2). Exact values are filled below.
His.pep_ch = repmat([1 2],length(His.mod_type),1);

% Theoretical m/z for each PTM×charge, using Mods and amino-acid masses
His.pep_mz = calculate_pepmz(His);

% Empirical RT priors (minutes) for the 5 forms; unmod in entry #1
His.rt_ref = [18.34; 21.58; 10.82; 10.80; 16.64];
His.display = ones(length(His.mod_type),1);

% Ensure a consistent "main" charge column across PTMs (align to first row col2)
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz); %#ok<ASGLU>
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end
    tune = [3 4]; % indices of PTMs to permute columns for consistent charge ordering
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end


function [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
        MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%% calculate_layout — Anchor unmodified RT, relocate PTM windows, extract RT/intensities
%
% STRATEGY
%   1) Anchor the unmodified form (#1) RT (His.rt_ref(1)):
%      • Optionally re-check via check_ref(...) (uses special.raw_path).
%      • If 'special.nDAmode==2', refine using MS2-assisted get_rts2.
%   2) Extract unmodified quantitative trace via get_histone0 (MS1 XIC).
%      If successful, RT-calibrate all other PTM priors by the observed delta.
%   3) Relocate expected windows for modified forms:
%      • relocate(...)  — MS1-only picking via get_rts
%      • relocate2(...) — MS2-assisted picking via get_rts2/get_rts22
%   4) For each PTM (2..5), call get_histone1(...) to fetch final RT/intensity.
%
[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);

pep_rts        = zeros([npep,ncharge]);     % per PTM × charge
pep_intens     = zeros([npep,ncharge]);     % per PTM × charge
mono_isointens = zeros([num_MS1,npep]);     % per scan × PTM (monoisotopic MS1 XIC)

% --- 1) Unmodified RT anchoring
His.rt_unmod_orig = His.rt_ref(1);  % keep original prior
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only mode: optionally re-check the prior RT against raw path metadata
        [His.rt_ref(1),special.ndebug] = check_ref( ...
            special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if 2==special.ndebug
            % Fallback: scan-wide search using get_rts for unmod and me1 to triangulate the anchor
            hno = 1; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts1,top1_rt1,inten_sum1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok<ASGLU>
            hno = 2; t1 = 0; t2 = MS1_index(num_MS1,2);
            [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,1,t1,t2); %#ok<ASGLU>
            if ~isempty(top1_rt2)
                if isempty(top1_rt1)
                    His.rt_ref(1) = top1_rt2 - 3;  % heuristic offset
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
    else
        % MS2-assisted anchoring of unmodified
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

% --- 2) Extract unmodified final RT/intensity and calibrate priors
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0( ...
    MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

if cur_rts(1) > 0
    % Observed RT for unmod anchor; shift other priors by delta
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1) - His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    % Store quantitative outputs
    pep_rts(hno,1:ncharge)               = cur_rts;
    pep_intens(hno,1:ncharge)            = cur_intens;
    mono_isointens(1:num_MS1,hno)        = cur_mono_isointens;
end

% --- 3) Relocate expected windows for the modified forms depending on mode
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);        % debug relocation (external)
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);     % MS1-only
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                        ptol,unitdiff,His,nhmass);                  % MS2-assisted
    end
end

% --- 4) Fetch final RT/intensity for each modified PTM (K4me1, me2, me3, ac)
for hno = 2:5
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone1( ...
        MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1) > 0
        pep_rts(hno,1:ncharge)        = cur_rts;
        pep_intens(hno,1:ncharge)     = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end


function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%% relocate — MS1-only relocation of expected RTs for K4 PTMs
%
% Uses get_rts over constrained RT windows around the unmodified anchor to
% propose final RT references for: K4me1, K4me2, K4me3, K4ac.
%
delta  = 0.1;  % small guard band
nsplit = 1;

% -- K4me1 expected after unmod (same charge ordering)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok<ASGLU>
His.rt_ref(hno) = iffEmptyUseZero(rts2, top1_rt2);

% -- K4me2 expected before unmod (earlier RT)
hno = 3;
t1 = 4;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok<ASGLU>

% -- K4me3 also before unmod; pair (me2,me3) by co-elution logic
hno = 4;
t1 = 4;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok<ASGLU>

[His.rt_ref(3),His.rt_ref(4)] = find_pair( ...
    rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
    rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% -- K4ac typically between me2 and unmod
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2); %#ok<ASGLU>
His.rt_ref(hno) = iffEmptyUseZero(rts5, top1_rt5);


function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%% relocate2 — MS2-assisted relocation for K4 PTMs
%
% Same intent as 'relocate', but MS2-confirmed windows are preferred:
%   • get_rts2 for K4me1 and K4ac
%   • get_rts22 for K4me2 / K4me3 (slightly different similarity voting)
%
delta  = 0.1;
nsplit = 1;

% -- K4me1 (after unmod)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+16;
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok<ASGLU>
His.rt_ref(hno) = iffEmptyUseZero(rts2, top1_rt2);

% -- K4me2 / K4me3 (before unmod), vote as a pair
hno = 3;
t1 = 4;
t2 = His.rt_ref(1)-3;
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = ...
    get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
              ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok<ASGLU>

hno = 4;
t1 = 4;
t2 = His.rt_ref(1)-3;
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = ...
    get_rts22(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
              ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok<ASGLU>

[His.rt_ref(3),His.rt_ref(4)] = find_pair( ...
    rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
    rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% -- K4ac (between me2 and unmod)
hno = 5;
t1 = (His.rt_ref(1)+His.rt_ref(3))/2;
t2 = His.rt_ref(1)-delta;
[rts5,top1_rt5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                           ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass); %#ok<ASGLU>
His.rt_ref(hno) = iffEmptyUseZero(rts5, top1_rt5);


% --- tiny utility for readability
function val = iffEmptyUseZero(rts, top1)
%% Return chosen RT or 0 if search failed
if isempty(rts)
    val = 0;
else
    val = top1;
end

