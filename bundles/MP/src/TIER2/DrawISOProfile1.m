function DrawISOProfile1(raw_path,raw_names,ptol,special)
%% DrawISOProfile1 — Main runner for histone_normal mode (MP bundle)
%
%  TIER: T2 (modified from upstream EpiProfile 2.0 basic)
%  ORIGIN: Forked from AT bundle for EpiProfile_PLANTS MP bundle.
%          MP-specific: H3_14_27_40 uses KSAPSTGGVKKPHR; H4_02v_20_23 adds KVFR.
%          In upstream EpiProfile 2.0, this function orchestrates the
%          per-run quantification of all histone peptide regions.
%          This version adds plant-specific TIER3 modules (sequence
%          variants derived from MSA) and robustness guards.
%
%  WHAT IT DOES:
%    1. Creates the output folder tree (histone_layouts/NN_rawname/detail/)
%    2. For each RAW file:
%       a) Loads MS1/MS2 scan and peak data
%       b) Detects acquisition mode (DDA vs DIA)
%       c) Calls every histone peptide module (H3, H4, H2A)
%          — TIER1 modules: upstream panels with full PTM quantification
%          — TIER3 modules: plant sequence-variant panels
%       d) Generates per-region Snapshots (H3, H4)
%    3. Aggregates results across runs:
%       — OutputTogether:    cohort-level histone_ratios.xls
%       — OutputSinglePTMs:  single-PTM marginal summaries
%       — OutputFigures:     QC figures (bar, boxplot, PCA, heatmaps)
%
%  INPUTS:
%    raw_path   — folder containing MS1/ and MS2/ subdirectories
%    raw_names  — cell array of run basenames (no extension)
%    ptol       — mass tolerance in ppm (typically 10 for FT)
%    special    — struct with fields: raw_path, nsource, nsubtype,
%                 norganism, soutput, nfigure, ndebug, nhmass
%
%  OUTPUTS:
%    Files written under raw_path/histone_layouts/:
%      — per-run .mat files with His/auc structs
%      — per-run XIC PDF plots in detail/
%      — histone_ratios.xls (cohort export)
%      — histone_ratios_single_PTMs.xls
%      — QC figures (if nfigure==1)

%% 1. Create output folder tree
layout_path = fullfile(raw_path,'histone_layouts');
if 0==exist(layout_path,'dir') && 0==mkdir(layout_path)
    fprintf(1,'can not create: %s\n',layout_path);
    return;
end;

for i=1:length(raw_names)
    if i<10
        prefix = ['0',num2str(i)];
    else
        prefix = num2str(i);
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',raw_names{i}]);
    if 0==exist(cur_outpath,'dir') && 0==mkdir(cur_outpath)
        fprintf(1,'can not create: %s\n',cur_outpath);
        return;
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',raw_names{i}],'detail');
    if 0==exist(cur_outpath,'dir') && 0==mkdir(cur_outpath)
        fprintf(1,'can not create: %s\n',cur_outpath);
        return;
    end;
end;

%% 2. Per-run quantification
for i=1:length(raw_names)
    fprintf(1,'\n%s\n',raw_names{i});
    cur_rawname = raw_names{i};

    % Load MS1 data
    MS1_scanfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1scans.mat']);
    MS1_peakfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1peaks.mat']);
    load(MS1_scanfile);% MS1_index
    load(MS1_peakfile);% MS1_peaks

    % Load MS2 data
    MS2_scanfile = fullfile(raw_path,'MS2',[cur_rawname,'_MS2scans.mat']);
    MS2_peakfile = fullfile(raw_path,'MS2',[cur_rawname,'_MS2peaks.mat']);
    load(MS2_scanfile);% MS2_index
    load(MS2_peakfile);% MS2_peaks

    % Detect DDA vs DIA
    nlen = length(unique(MS2_index(:,4)));%#ok
    c0 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+0,4) && MS2_index(2,4)==MS2_index(2+nlen+0,4);
    c1 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+1,4) && MS2_index(2,4)==MS2_index(2+nlen+1,4);
    c2 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+2,4) && MS2_index(2,4)==MS2_index(2+nlen+2,4);
    c3 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+3,4) && MS2_index(2,4)==MS2_index(2+nlen+3,4);
    c4 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+4,4) && MS2_index(2,4)==MS2_index(2+nlen+4,4);
    if nlen<270 && (c0 || c1 || c2 || c3 || c4)% (1100-300)/3=267
        special.nDAmode = 2;% DIA
    else
        special.nDAmode = 1;% DDA
    end;

    % Build output path for this run
    if i<10
        prefix = ['0',num2str(i)];
    else
        prefix = num2str(i);
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');

    %% ---- TIER1: Upstream H3 modules (PTM quantification) ----
    H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_02a_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_02b_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_04a_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_04v3a_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_05_41_49(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_06_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_06a_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_09u_64_135(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);

    %% ---- TIER2: Plant-adapted H3 modules ----
    H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);

    %% ---- TIER3: Plant sequence-variant modules (from MSA) ----
    H3_11_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_12_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_13_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_14_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_16_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_17_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3_18_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);

    %% ---- TIER1: Upstream H4 modules ----
    H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_02v_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_02a_18_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_02b_20_35(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_03_24_35(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_04_40_45(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_05_68_78(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_06_79_92(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4_07u_24_102(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);

    %% ---- TIER2: Plant H2A unmod diagnostic peptides ----
    HH2A_01u_1_7(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);

    %% ---- TIER2: Plant H2B unmod diagnostic peptides ----
    HH2B_01u_104_145(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);

    %% ---- TIER3: Plant H3 Snapshot (aggregator) ----
    H3_Snapshot(cur_outpath);
    H4_Snapshot(cur_outpath);

    fprintf(1,'done.\n');
end;

%% 3. Cohort-level aggregation
OutputTogether(layout_path,raw_names);
OutputSinglePTMs(layout_path,raw_names);

if 1==special.nfigure
    OutputFigures(layout_path,raw_names);
end;
