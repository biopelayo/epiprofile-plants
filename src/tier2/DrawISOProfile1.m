function DrawISOProfile1(raw_path,raw_names,ptol,special)
%%
% ========================================================================
% DrawISOProfile1
% ------------------------------------------------------------------------
% PLAIN-ENGLISH PURPOSE (non-proteomics audience):
%   Orchestrate a full per-sample histone peptide profiling run:
%   - Prepare/ensure an output directory structure under "histone_layouts/"
%   - For each raw sample name:
%       * Load precomputed MS1/MS2 scan indices and peak lists (MAT files)
%       * Heuristically decide if data is DIA-like vs DDA-like
%       * Dispatch to a series of peptide-segment extractor functions
%         (e.g., H3_01_3_8, H3_02_9_17, H4_01_4_17, etc.)
%       * Optionally recalculate within-window ratios across variant splits
%         (e.g., H3_02/H3_02a/H3_02b or H3_04/H3_04a)
%       * Produce per-sample snapshots/benchmarks
%   - Finally, aggregate all outputs:
%       * OutputTogether: merged table of RT/Area/Ratio across samples
%       * OutputSinglePTMs: single-PTM ratio table (for normal/SILAC sources)
%       * OutputFigures: optional multi-sample visual summaries
%
% WHAT THIS FUNCTION ASSUMES EXISTS:
%   For each sample name 'X' in raw_names, four MAT files already exist:
%     MS1_index/peaks → "<raw_path>/MS1/X_MS1scans.mat" / "X_MS1peaks.mat"
%     MS2_index/peaks → "<raw_path>/MS2/X_MS2scans.mat" / "X_MS2peaks.mat"
%   Each MAT contains the expected variables (e.g., MS1_index, MS1_peaks, ...).
%
% KEY INPUTS:
%   raw_path   : Base input folder with subfolders 'MS1' and 'MS2' holding MAT files
%   raw_names  : Cell array of sample base names (strings without extensions)
%   ptol       : Mass tolerance setting used by downstream extractors (unit depends on pipeline)
%   special    : Struct carrying run-wide switches, e.g.:
%                  .nsource   → 1: histone_normal, 2: histone_SILAC,
%                                3: histone_C13, 4: histone_N15, 5: histone_13CD3
%                  .nsubtype  → subtype for N15 (0..5, see below)
%                  .soutput   → character array/string controlling which H3/H4 panels to run
%                  .nfigure   → if true(1) and enough samples, draw summary figures
%                  .nDAmode   → will be SET HERE per sample (1=DDA, 2=DIA) via heuristic
%
% OUTPUTS:
%   No direct variables returned; side effects:
%     - Creates "histone_layouts/<##>_<sample>/detail/" directories
%     - Writes many .mat/.pdf/.xls outputs via the called H3_*/H4_* functions
%     - Calls OutputTogether / OutputSinglePTMs / OutputFigures to aggregate
%
% NOTES:
%   * The "H3_*/H4_*" function names encode histone/proteolytic window
%     (e.g., H3_02_9_17 stands for H3, window residues 9–17).
%   * The "recal_ratio_*" helpers re-normalize the third column (Ratio)
%     across split windows by recomputing area fractions.
%   * This function never modifies raw data; it only reads and orchestrates.
% ========================================================================

% ------------------------------------------------------------------------
% Ensure base layout folder exists: "<raw_path>/histone_layouts"
%   If it doesn’t exist, try to create it; abort on failure.
layout_path = fullfile(raw_path,'histone_layouts');
if 0==exist(layout_path,'dir') && 0==mkdir(layout_path)
    fprintf(1,'can not create: %s\n',layout_path);
    return;
end;

% ------------------------------------------------------------------------
% Pre-create per-sample output folders:
%   For each raw name, make:
%     "<layout_path>/<##>_<name>/" and ".../detail/"
%   The <##> is a zero-padded 2-digit running index (01, 02, ...).
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

% ------------------------------------------------------------------------
% MAIN PER-SAMPLE LOOP:
%   For each sample:
%     - Load required MS1/MS2 *index* and *peak* MAT files
%     - Decide DIA vs DDA by a simple heuristic (sets special.nDAmode)
%     - Recompute per-window panels depending on 'nsource' and 'soutput'
for i=1:length(raw_names)
    fprintf(1,'\n%s\n',raw_names{i});
    cur_rawname = raw_names{i};

    % --- Resolve MAT file paths for this sample and load them --------------
    MS1_scanfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1scans.mat']);
    MS1_peakfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1peaks.mat']);
    MS2_scanfile = fullfile(raw_path,'MS2',[cur_rawname,'_MS2scans.mat']);
    MS2_peakfile = fullfile(raw_path,'MS2',[cur_rawname,'_MS2peaks.mat']);
    load(MS1_scanfile);% expects variable: MS1_index
    load(MS1_peakfile);% expects variable: MS1_peaks
    load(MS2_scanfile);% expects variable: MS2_index
    load(MS2_peakfile);% expects variable: MS2_peaks

    % --- Heuristic DIA vs DDA decision ------------------------------------
    % nlen = number of unique precursor m/z values observed in MS2_index(:,4).
    % The code then checks a few "offset alignments" c0..c4 to see if patterns
    % repeat every nlen scans (indicative of systematically stepped isolation
    % windows as seen in DIA). Threshold nlen<270 is used (empirical).
    nlen = length(unique(MS2_index(:,4)));%#ok
    c0 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+0,4) && MS2_index(2,4)==MS2_index(2+nlen+0,4);
    c1 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+1,4) && MS2_index(2,4)==MS2_index(2+nlen+1,4);
    c2 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+2,4) && MS2_index(2,4)==MS2_index(2+nlen+2,4);
    c3 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+3,4) && MS2_index(2,4)==MS2_index(2+nlen+3,4);
    c4 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+4,4) && MS2_index(2,4)==MS2_index(2+nlen+4,4);
    if nlen<270 && (c0 || c1 || c2 || c3 || c4)% (1100-300)/3=267
        special.nDAmode = 2;% Interpret as DIA-like
    else
        special.nDAmode = 1;% Interpret as DDA-like
    end;

    % --- Build per-sample output detail folder -----------------------------
    if i<10
        prefix = ['0',num2str(i)];
    else
        prefix = num2str(i);
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');

    % ======================================================================
    % DISPATCH BY SOURCE TYPE:
    %   special.nsource:
    %     1 → histone_normal
    %     2 → histone_SILAC
    %     3 → histone_C13
    %     4 → histone_N15 (subtypes apply)
    %     5 → histone_13CD3
    %   Within 1/2 we further branch by special.soutput(1):
    %     '1' → "H3H4 basic"
    %     '2' → "H3H4 basic + H3S10ph" (+ a ratio recalibration)
    %     else → "H3H4 all" (+ multiple ratio recalibrations)
    % ======================================================================
    if 1==special.nsource || 2==special.nsource
        % histone_normal or histone_SILAC
        if '1'==special.soutput(1)
            % --- CASE '1': H3H4 basic panels only -------------------------
            H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_06_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_11_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_12_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_13_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_14_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_16_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_17_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_18_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
            
            H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_04_40_45(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');

        elseif '2'==special.soutput(1)
            % --- CASE '2': H3H4 basic + H3S10ph --------------------------
            H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02a_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_06_54_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');

            % Re-normalize combined ratios for H3_02 windows (02 + 02a)
            recal_ratio_H3_0201(cur_outpath);
            
            H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_04_40_45(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');

        else
            % --- CASE 'else': H3H4 all (comprehensive panels) -------------
            H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02a_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02b_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);

            H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04a_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04v3a_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_05_41_49(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_06_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_06a_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_09u_64_135(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');

            % Re-normalize combined ratios for H3_02 (02+02a+02b) and H3_04 (04+04a and 04v3+04v3a)
            recal_ratio_H3_0202(cur_outpath);
            recal_ratio_H3_0402(cur_outpath);
            
            H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02a_18_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02b_20_35(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02c_20_36(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_03_24_35(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_04_40_45(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_05_68_78(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_06_79_92(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_07u_24_102(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
        end;

        % --- Post-panels: build quick summaries for H3 and H4 --------------
        H3_Snapshot(cur_outpath);
        H4_Snapshot(cur_outpath);
        
        % --- Optional: subset of histone H2A panels based on soutput(2) ----
        if '1'==special.soutput(2)
            HH2A_01u_1_7(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        end;
        
        % --- If we concluded this was DDA data, produce benchmark outputs ---
        if 1==special.nDAmode
            GetBenchmark(cur_outpath);
        end;

    elseif 3==special.nsource
        % ==================================================================
        % histone_C13 source: specific subset of panels (H3, H4, H2A)
        % ==================================================================
        H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        HH2A_02m1_4_11(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        HH2A_04oV_1_19(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        HH2A_04oZ_1_19(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        HH2A_05m1_12_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        fprintf(1,'\n');

    elseif 4==special.nsource
        % ==================================================================
        % histone_N15 source with subtypes:
        %   nsubtype: 0=N14 light Mods, 1=N15 light Mods, 2=N14 heavy Mods,
        %             3=N15 heavy Mods, 4=0+1, 5=0+3 (combinations)
        % ==================================================================
        if 0==special.nsubtype || 4==special.nsubtype || 5==special.nsubtype
            % Unlabeled or combined light/heavy mixes
            H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_02m1_4_11(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_05m3_12_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
        else
            % Explicit N15 variants (H3N/H4N/HH2AN panels)
            H3N_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3N_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3N_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3N_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3N_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4N_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4N_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2AN_02m1_4_11(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2AN_05m3_12_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
        end;

    elseif 5==special.nsource
        % ==================================================================
        % histone_13CD3 source: subset tailored to this labeling strategy
        % ==================================================================
        H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_06_54_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        fprintf(1,'\n');
        
        H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H4_06_79_92(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        fprintf(1,'\n');
    end;
end;

% ------------------------------------------------------------------------
% FINAL AGGREGATIONS ACROSS ALL SAMPLES:
%   - OutputTogether: single spreadsheet with (by default) columns for
%     Ratios, Areas, RTs per peptide row across samples.
%   - OutputSinglePTMs: per-PTM ratio tables (only for sources 1 or 2).
%   - OutputFigures: optional figure summaries if enough samples and enabled.
OutputTogether(layout_path,raw_names);
if 1==special.nsource || 2==special.nsource
    OutputSinglePTMs(layout_path,raw_names);
end;

if (1==special.nsource || 2==special.nsource) && 1==special.nfigure && length(raw_names)>2
    OutputFigures(layout_path,raw_names);
end;

% ========================================================================
% LOCAL HELPERS — Re-normalize ratio columns within split windows
% ------------------------------------------------------------------------
% These functions re-open the .mat results for specific windows (e.g., H3_02
% and its variants H3_02a/H3_02b; also H3_04/H3_04a and H3_04v3/H3_04v3a),
% compute the sum of Areas (auc(:,2)) across the set, and then recompute
% auc(:,3) = auc(:,2) / total_sum so that column 3 holds a within-window
% fraction (i.e., normalized ratio) after combining variants.
% ========================================================================

function recal_ratio_H3_0201(cur_outpath)
%%
% Combine H3_02 + H3_02a: set Ratio = Area / (Area_02 + Area_02a)

% H3_02
out_filename = 'H3_02_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s1 = sum(auc(:,2));%#ok

out_filename = 'H3_02a_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s2 = sum(auc(:,2));

s = s1+s2;

out_filename = 'H3_02_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

out_filename = 'H3_02a_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;%#ok
save(mat_file,'His','auc');

function recal_ratio_H3_0202(cur_outpath)
%%
% Combine H3_02 + H3_02a + H3_02b: set Ratio = Area / (Area_02 + Area_02a + Area_02b)

% H3_02
out_filename = 'H3_02_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s1 = sum(auc(:,2));%#ok

out_filename = 'H3_02a_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s2 = sum(auc(:,2));

out_filename = 'H3_02b_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s3 = sum(auc(:,2));

s = s1+s2+s3;

out_filename = 'H3_02_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

out_filename = 'H3_02a_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

out_filename = 'H3_02b_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;%#ok
save(mat_file,'His','auc');

function recal_ratio_H3_0402(cur_outpath)
%%
% Combine H3_04 + H3_04a  (and)  H3_04v3 + H3_04v3a
% For each pair, set Ratio = Area / (Area_sum_pair)

% H3_04
out_filename = 'H3_04_27_40';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s1 = sum(auc(:,2));%#ok

out_filename = 'H3_04a_27_40';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s2 = sum(auc(:,2));

s = s1+s2;

out_filename = 'H3_04_27_40';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

out_filename = 'H3_04a_27_40';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

% H3_04v3
out_filename = 'H3_04v3_27_40';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s1 = sum(auc(:,2));

out_filename = 'H3_04v3a_27_40';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s2 = sum(auc(:,2));

s = s1+s2;

out_filename = 'H3_04v3_27_40';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

out_filename = 'H3_04v3a_27_40';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;%#ok
save(mat_file,'His','auc');
