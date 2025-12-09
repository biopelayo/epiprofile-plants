function H3_05_41_49(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%% H3_05_41_49 — Targeted panel for H3 peptide YRPGTVALR (residues 41–49)
%
% Purpose
% -------
% This routine performs targeted extraction and quantification of the short H3
% peptide "YRPGTVALR" (aa 41–49 canonical numbering) across three PTM states:
%   1) unmodified (baseline),
%   2) Y41pr (historical shorthand: "pr" stands for propionylation in EpiProfile),
%   3) Y41ph (tyrosine phosphorylation).
%
% The panel:
%   - Initializes a histone/PTM layout (precursor m/z per charge, reference RTs).
%   - Re-locates expected RTs using MS1 (and optionally MS2, if DA mode is active).
%   - Extracts monoisotopic chromatograms and summarizes intensities/RTs per PTM.
%   - Saves results to disk and draws layout figures.
%   - Optionally exports PSMs when diagnostic-assisted mode is requested.
%
% Inputs
% ------
% MS1_index : [nMS1 x 2] double
%   Scan index (col 1) and retention time in minutes (col 2) for MS1 spectra.
%
% MS1_peaks : struct/cell (EpiProfile internal)
%   Centroided MS1 peaks indexed by scan. Used by low-level extractors.
%
% MS2_index : [nMS2 x 2] double
%   Scan index (col 1) and retention time in minutes (col 2) for MS2 spectra.
%
% MS2_peaks : struct/cell (EpiProfile internal)
%   Centroided MS2 peaks for diagnostic fragment checks (used in DA mode).
%
% ptol : double
%   Mass tolerance (typically ppm) forwarded to the extraction helpers.
%
% cur_outpath : char
%   Output directory for .mat artifacts and figures.
%
% special : struct with control flags and optional parameters
%   Fields commonly consumed here:
%     - ndebug    : 1 → use relocateD (diagnostic/simple relocation).
%     - nDAmode   : implementation-dependent; this panel treats
%                   (nDAmode ~= 2) as “MS1-only mode” and (nDAmode == 2) as
%                   “MS2-assisted DA mode” for the relocation stage.
%                   Separately, if (nDAmode == 1) we trigger PSM export below.
%                   (Yes, legacy cores used mixed semantics; we preserve them.)
%     - nhmass    : neutral-loss/precursor mass hint for get_rts2 (DA mode).
%     - raw_path  : raw vendor path used by check_ref when adjusting unmodified RT.
%
% Side-effects (disk outputs)
% ---------------------------
% - <cur_outpath>/H3_05_41_49.mat  (tabular quant for this panel)
% - Layout figures (XIC plots/heatmaps depending on your draw_layout build)
% - Optional PSM exports if (special.nDAmode == 1)
%
% Dependencies (not defined here; part of EpiProfile core)
% --------------------------------------------------------
%   calculate_pepmz, get_histone0, get_histone1, get_rts, get_rts2,
%   relocateD, output_histone, draw_layout, GetPSM, check_ref
%
% Notes
% -----
% - We do not change business logic; only add protective comments and minor
%   clarifications. The “pr” code in mod_type follows the EpiProfile convention
%   for propionylation (chemical derivatization). If your pipeline uses
%   different derivatization, update mod_type accordingly.
% - Charge grid intentionally includes 1+/2+/3+ because this short peptide can
%   produce significant singly charged signal in some TOF/low-mass settings.

    % ------------------------ %
    % 0) Guard: already done? %
    % ------------------------ %
    out_filename = 'H3_05_41_49';
    fprintf(1,'%s..',out_filename);
    out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
    if 0~=exist(out_file0,'file')
        % If results already exist, exit silently to avoid recomputation.
        return;
    end

    % ------------------------- %
    % 1) Initialize PTM layout %
    % ------------------------- %
    His = init_histone();

    % ------------------------------------------- %
    % 2) Compute: relocate + extract intensities  %
    % ------------------------------------------- %
    unitdiff = 1.0032; % ^13C isotopic spacing
    [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
        MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special);

    % -------------------------- %
    % 3) Persist tabular result %
    % -------------------------- %
    output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

    % ------------------------ %
    % 4) Draw summary figures  %
    % ------------------------ %
    num_MS1 = size(MS1_index,1);
    isorts  = MS1_index(1:num_MS1,2); % RT axis for plotting
    draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
                isorts,mono_isointens,MS2_index,MS2_peaks,special);

    % ------------------------------------------- %
    % 5) Optional: export PSMs for DA mode (n=1) %
    % ------------------------------------------- %
    if 1==special.nDAmode
        GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
               mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
    end
end


% ====================================================================== %
%                              SUBFUNCTIONS                              %
% ====================================================================== %

function His = init_histone()
%% Assemble the "His" struct with peptide/PTM definition and charge grid
%  Fields:
%    - pep_seq    : base sequence (YRPGTVALR)
%    - mod_short  : human-readable PTM names (3 rows)
%    - mod_type   : compact encoding used by EpiProfile extractors
%    - pep_ch     : charge grid (replicated [1 2 3] for each PTM row)
%    - pep_mz     : precursor m/z matrix [nPTM x nCharge], calculated
%    - rt_ref     : seed/reference RTs used to initialize search windows
%    - display    : visibility mask for plotting (e.g., hide Y41pr by default)
%
%  About mod_type:
%    - "idx,tag;" pairs target residue indices in EpiProfile’s internal mapping
%      (not redefined here). "pr" = propionylation, "ph" = phosphorylation.
%    - The first row (unmod) only contains peptide-level propionylation (if any).
%    - The second/third rows add residue-level PTMs (Y41pr / Y41ph).

    His.pep_seq = 'YRPGTVALR';

    His.mod_short = { ...
        'unmod'; ...
        'Y41pr'; ...
        'Y41ph' ...
    };

    % NOTE: Keep this aligned with your core’s indexing. We do not alter it here.
    His.mod_type = { ...
        '0,pr;';        ... % unmodified peptide (baseline derivatization only)
        '0,pr;1,pr;';   ... % Y41 propionylation (legacy shorthand)
        '0,pr;1,ph;'    ... % Y41 phosphorylation
    };

    % Allow 1+/2+/3+ (short peptide often strong at 1+ in some instruments)
    His.pep_ch = repmat([1 2 3], length(His.mod_type), 1);

    % If you ever need hard-coded m/z, leave them in the commented block above.
    % We keep it dynamic here:
    His.pep_mz = calculate_pepmz(His);

    % Initial RT references (minutes). These are seed anchors only.
    His.rt_ref = [ ...
        31.20;  ... unmod
        38.18;  ... Y41pr (often later due to derivatization)
        34.00   ... Y41ph (can shift earlier/later depending on system)
    ];

    % Control layout visibility (hide Y41pr in plots; still quantified)
    His.display = ones(length(His.mod_type),1);
    His.display(2) = 0;

    % Reorder charge columns such that the “main” charge (middle, i.e. 2+) is first
    main_ch = His.pep_ch(1,2);
    if main_ch ~= His.pep_ch(1,1)
        [npep,ncharge] = size(His.pep_mz);
        new_ch = [main_ch, setdiff(His.pep_ch(1,:), main_ch)];
        x = zeros([1,ncharge]);
        for ino=1:ncharge
            x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
        end
        tune = 1:npep;
        His.pep_mz(tune,:) = His.pep_mz(tune,x);
        His.pep_ch(tune,:) = His.pep_ch(tune,x);
    end
end


function [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%% Relocate the expected RTs and extract intensities/RTs for each PTM row
%
% Outputs
% -------
% pep_rts         : [nPTM x nCharge] RT per PTM/charge (0 if no peak found)
% pep_intens      : [nPTM x nCharge] XIC areas or intensity proxies
% mono_isointens  : [nMS1 x nPTM]    monoisotopic XIC over RT (for plotting)
%
% Strategy
% --------
% 1) Try to refine the unmodified anchor (His.rt_ref(1)) via check_ref:
%    - In MS1-only mode, this may adjust both the RT and ndebug flag.
%    - In DA mode (here defined as nDAmode==2 for relocation stage), optionally
%      run get_rts2 around the current RT to refine the anchor further.
% 2) Extract unmodified signal via get_histone0 and “calibrate” the RT grid:
%    - If unmodified is found, delta = (found_RT - original_unmod_RT) is added
%      to all other RT references so that downstream windows track the shift.
% 3) Relocate PTM windows and extract PTMs Y41pr/Y41ph using:
%    - relocate  (MS1-only)   or
%    - relocate2 (MS2-assisted).
%
% Notes
% -----
% - get_histone0/get_histone1 are core EpiProfile helpers that, given His and an
%   index (hno), perform XIC integration/peak picking for the unmodified/PTM row.

    [npep,ncharge] = size(His.pep_mz);
    num_MS1 = size(MS1_index,1);

    pep_rts        = zeros([npep,ncharge]);
    pep_intens     = zeros([npep,ncharge]);
    mono_isointens = zeros([num_MS1,npep]);

    % ---- 1) Unmodified anchor check/refinement ----
    His.rt_unmod_orig = His.rt_ref(1);
    if 1~=special.ndebug
        if 2~=special.nDAmode
            % MS1-only: check/adjust unmodified reference (may also update ndebug)
            [His.rt_ref(1), special.ndebug] = ...
                check_ref(special.raw_path, [His.pep_seq,His.mod_type{1}], His.rt_ref(1), special.ndebug);
        else
            % DA mode (relocation via MS2 later): still adjust unmodified anchor
            nhmass = special.nhmass;
            His.rt_ref(1) = check_ref(special.raw_path, [His.pep_seq,His.mod_type{1}], His.rt_ref(1), special.ndebug);

            % If check_ref yields no change, widen the search around unmodified RT
            if His.rt_unmod_orig==His.rt_ref(1)
                t1 = 0;
                t2 = MS1_index(num_MS1,2);
            else
                delta = 5; % +/- 5 min around refined anchor
                t1 = His.rt_ref(1)-delta;
                t2 = His.rt_ref(1)+delta;
            end
            hno = 1; % unmodified row
            [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                                       ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok<ASGLU>
            if ~isempty(top1_rt1)
                His.rt_ref(1) = top1_rt1;
            end
        end
    end

    % ---- 2) Extract unmodified and calibrate RT grid ----
    hno = 1; % unmodified
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

    if cur_rts(1)>0
        % Calibrate the reference RTs so all downstream windows track the shift
        His.rt_ref(1) = cur_rts(1);
        delta = cur_rts(1) - His.rt_unmod_orig;
        His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

        pep_rts(hno,1:ncharge)            = cur_rts;
        pep_intens(hno,1:ncharge)         = cur_intens;
        mono_isointens(1:num_MS1,hno)     = cur_mono_isointens;
    end

    % ---- 3) Relocation for modified states (MS1-only vs DA) ----
    if 1==special.ndebug
        His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        if 2~=special.nDAmode
            His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
        else
            % nhmass was defined above in the DA preparation branch
            His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
        end
    end

    % ---- 4) Extract modified rows (Y41pr, Y41ph) ----
    for hno = 2:3
        [cur_rts,cur_intens,cur_mono_isointens] = ...
            get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
        if cur_rts(1)>0
            pep_rts(hno,1:ncharge)        = cur_rts;
            pep_intens(hno,1:ncharge)     = cur_intens;
            mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
        end
    end
end


function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%% MS1-only relocation of expected RT windows for Y41pr / Y41ph
%
% Heuristic windows
% -----------------
% - We anchor both modified states to the refined unmodified RT (His.rt_ref(1))
%   and search downstream. This peptide is short and often exhibits clean RT
%   alignment; tight deltas suffice (0.1 min).
%
% Notes
% -----
% - nsplit is left at 1 (no special multi-split regime).
% - If a row returns empty from get_rts, we set its RT to 0.

    delta  = 0.1;
    nsplit = 1;

    % --- Y41pr ---
    hno = 2;
    t1 = His.rt_ref(1) + delta;
    t2 = His.rt_ref(1) + 25;
    [rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

    if isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end

    % --- Y41ph ---
    hno = 3;
    t1 = His.rt_ref(1) + delta;
    t2 = His.rt_ref(1) + 15;
    [rts3,top1_rt3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

    if isempty(rts3)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt3;
    end
end


function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%% MS2-assisted relocation (DA mode) for Y41pr / Y41ph
%
% Same windows as MS1-only but relying on get_rts2() to inject diagnostic
% evidence from MS2 (neutral losses/fragment ions) when available. This
% improves specificity under co-elution or low S/N.

    delta  = 0.1;
    nsplit = 1;

    % --- Y41pr ---
    hno = 2;
    t1 = His.rt_ref(1) + delta;
    t2 = His.rt_ref(1) + 25;
    [rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                               ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

    if isempty(rts2)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt2;
    end

    % --- Y41ph ---
    hno = 3;
    t1 = His.rt_ref(1) + delta;
    t2 = His.rt_ref(1) + 15;
    [rts3,top1_rt3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                               ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

    if isempty(rts3)
        His.rt_ref(hno) = 0;
    else
        His.rt_ref(hno) = top1_rt3;
    end
end
