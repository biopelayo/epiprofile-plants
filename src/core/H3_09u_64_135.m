function H3_09u_64_135(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H3_09u_64_135
% -------------------------------------------------------------------------
% Purpose:
%   Utility panel quantifying three short H3 peptides in one run:
%     1) 64-69  KLPFQR   (mid-late elution anchor)
%     2) 129-134 RIRGER  (early elution)
%     3) 130-135 IRGERA  (early-mid, close to RIRGER)
%
%   Unlike PTM-centric panels, each "row" here is a different peptide
%   sequence. This is useful for RT anchoring, digestion QC, and assessing
%   chromatographic stability.
%
% Inputs:
%   MS1_index, MS1_peaks : centroided MS1 (core format)
%   MS2_index, MS2_peaks : centroided MS2 (only used to refine RT seeds
%                          when nDAmode == 2)
%   ptol                 : mass tolerance (ppm/instrument units)
%   cur_outpath          : output directory
%   special              : struct
%       .ndebug   -> 1: relocateD (diagnostic)
%       .nDAmode  -> historical semantics:
%                     2 => refine seeds using MS2 (get_rts2)
%                     1 => export PSMs after quantification (GetPSM2)
%       .nhmass   -> neutral-loss/precursor mass hint for get_rts2
%       .raw_path -> vendor raw path (check_ref)
%
% Outputs (to disk):
%   <cur_outpath>/H3_09u_64_135.mat  (pep_rts, pep_intens, mono_isointens)
%   Figures via draw_layout; PSM tables via GetPSM2 if nDAmode == 1
% -------------------------------------------------------------------------

% ------------------------------- check -----------------------------------
out_filename = 'H3_09u_64_135';
% fprintf(1,'%s..',out_filename); % progress print kept commented like original
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    % Idempotency: if a prior run exists, skip computation.
    return;
end

% -------------------------------- init -----------------------------------
His = init_histone();

% ------------------------------ calculate --------------------------------
unitdiff = 1.0032; % ~13C spacing for isotopic handling
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special);

% -------------------------------- output ---------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% --------------------------------- draw ----------------------------------
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2); % RT axis (minutes)
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
            mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------- Get PSM ---------------------------------
if 1==special.nDAmode
    % Historical choice for this utility panel: use GetPSM2.
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
            mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % main

% =========================================================================
function His = init_histone()
%%
% Build the "His" struct for a multi-peptide utility panel.
% Here, pep_seq is a sentinel ('unmod'), and the actual peptide sequences
% are encoded in 'mod_short'. 'mod_type' captures derivatization rules per
% core conventions (do not change unless your core maps indices differently).

His.pep_seq = 'unmod'; % sentinel; actual sequences live in mod_short

% Each entry in mod_short is a distinct peptide sequence to be quantified.
His.mod_short = {'KLPFQR';   % 64-69
                 'RIRGER';   % 129-134
                 'IRGERA'};  % 130-135

% Compact PTM encodings (core-dependent):
% - '0,pr' => peptide termini derivatization
% - '1,pr' => extra terminal handling for first peptide (historical)
His.mod_type = {'0,pr;1,pr;';  % KLPFQR
                '0,pr;';       % RIRGER
                '0,pr;'};      % IRGERA

% Allow 1+/2+/3+ for each row
His.pep_ch = repmat([1 2 3],length(His.mod_type),1);

% Legacy hard-coded m/z retained as comment for cross-checking only:
%{
His.pep_mz = [900.5302  450.7687  300.8482;  % KLPFQR
              842.4955  421.7514  281.5034;  % RIRGER
              757.4315  379.2194  253.1487]; % IRGERA
%}

% Compute m/z values portably from composition + mod_type
His.pep_mz = calculate_pepmz(His);

% Seed RTs (min) per peptide row
His.rt_ref = [39.63;  % KLPFQR
              17.36;  % RIRGER
              22.00]; % IRGERA

% Utility panel: hidden in plots by default (still quantified)
His.display = zeros(length(His.mod_type),1);

% ----------------------- charge column reordering ------------------------
% Panel-wide convention: make 2+ the first column for consistency.
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

end % init_histone

% =========================================================================
function [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% High-level orchestration for a multi-peptide utility panel:
%   - For each row (peptide), refine its RT seed via check_ref, optionally
%     probing with get_rts2 (DA mode) around that seed.
%   - Extract each peptide via get_histone0 (MS1-driven extraction).
%   - Build matrices: pep_rts, pep_intens, mono_isointens.

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% ---------------------- refine RT seeds (all rows) -----------------------
His.rt_unmod_orig = His.rt_ref(1); % retained for consistency (unused below)

if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only seed refinement using check_ref
        for hno=1:3
            [His.rt_ref(hno),special.ndebug] = check_ref( ...
                special.raw_path,[His.mod_short{hno},His.mod_type{hno}], ...
                His.rt_ref(hno),special.ndebug);
        end
    else
        % DA-assisted seed refinement: after check_ref, probe with get_rts2
        nhmass = special.nhmass;
        for hno=1:3
            rt_unmod_orig = His.rt_ref(hno);
            His.rt_ref(hno) = check_ref( ...
                special.raw_path,[His.mod_short{hno},His.mod_type{hno}], ...
                His.rt_ref(hno),special.ndebug);

            % Determine probing window: full-range if unchanged; otherwise ±5 min
            if rt_unmod_orig==His.rt_ref(hno)
                t1 = 0;
                t2 = MS1_index(num_MS1,2);
            else
                delta = 5;
                t1 = His.rt_ref(hno)-delta;
                t2 = His.rt_ref(hno)+delta;
            end

            % Probe with DA (MS2-assisted)
            [~,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                                    ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
            if 0==isempty(top1_rt1)
                His.rt_ref(hno) = top1_rt1;
            end
        end
        % Set ndebug to 1 afterwards (historical behavior preserved)
        special.ndebug = 1;
    end
end

% ----------------------------- extraction --------------------------------
% Extract each peptide row using get_histone0 (MS1-driven).
% Rows: 1) KLPFQR  2) RIRGER  3) IRGERA
for hno=1:3
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone0( ...
        MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

end % calculate_layout
