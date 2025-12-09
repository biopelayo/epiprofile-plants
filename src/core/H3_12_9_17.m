function H3_12_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% === Purpose (English) =======================================================
% Targeted quantification and visualization for four H3 peptides (see
% init_histone.mod_short) across charge states defined in His.pep_ch.
% Workflow:
%   1) Skip if results already exist (idempotent execution).
%   2) Initialize the histone-peptide panel (sequences/mods/charges/mz/RTs).
%   3) Compute per-peptide retention times and intensities using MS1/MS2.
%   4) Persist results to disk and draw diagnostic figures.
%   5) Optionally, extract PSM-level evidence in data-assisted mode.

% check
out_filename = 'H3_12_9_17';                                 % fixed output stem (used for files)
% fprintf(1,'%s..',out_filename);                            % optional console logging (kept muted)
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);     % expected .mat result path
if 0~=exist(out_file0,'file')                                 % if the result file already exists…
    return;                                                   % …exit early to avoid recomputation
end;

% init
His = init_histone();                                         % build the peptide panel specification:
                                                              %   - mod_short/mod_type (derivatization map)
                                                              %   - pep_ch (charge states)
                                                              %   - pep_mz (m/z per charge, computed)
                                                              %   - rt_ref (reference RT anchors)
                                                              %   - display flags

% calculate
unitdiff = 1.0032;                                            % isotopic mass diff (~C13-C12), used to
                                                              % detect monoisotopic/isotopic spacing
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);              % main quantification (RTs, intensities,
                                                              % monoisotopic traces)

% output
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts); % persist His + quantification to .mat

% draw
num_MS1 = size(MS1_index,1);                                  % number of MS1 rows/scans
isorts = MS1_index(1:num_MS1,2);                              % RT column used for sorting/x-axis
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special); % render diagnostic plots

% Get PSM
if 1==special.nDAmode                                          % if data-assisted mode enabled (1)
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
             isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff); % extract PSM support
end;

% ============================================================================
% Local function: histone-peptide panel and static metadata
% ============================================================================
function His = init_histone()
%%
% === Purpose ================================================================
% Define peptide labels, modification encodings, charge states, theoretical
% m/z values, reference RTs and display flags for the H3 panel used here.

His.pep_seq = 'unmod';                                        % panel descriptor (unmodified backbone)
His.mod_short = {'KSTGGKAPR';                                 % four short peptide labels
    'KSTGGKGPR';
    'KSHGGKAPR';
    'ISTGGKAPR'};
His.mod_type = {'0,pr;1,pr;6,pr;';                            % modification map per peptide:
    '0,pr;1,pr;6,pr;';                                        % 'pr' typically = propionylation (see README)
    '0,pr;1,pr;6,pr;';
    '0,pr;6,pr;'};                                            % positions indexed in peptide coordinates

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);          % consider charge states z = 1, 2, 3
%{
His.pep_mz = [900.5302	450.7687	300.8482                   % OPTIONAL hard-coded m/z table
    842.4955	421.7514	281.5034
    757.4315	379.2194	253.1487];
%}
His.pep_mz = calculate_pepmz(His);                             % compute m/z per charge from seq+mods
His.rt_ref = [39.63                                           % reference retention times (minutes)
    17.36
    22
    22];
His.display = ones(length(His.mod_type),1);                    % 1 = include peptide in outputs/plots

% main ch
main_ch = His.pep_ch(1,2);                                    % choose charge 2 as the "main" charge
if main_ch~=His.pep_ch(1,1)                                   % if main charge is not leftmost column…
    [npep,ncharge] = size(His.pep_mz);                        % …we reorder columns so main_ch is first
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];      % new column order: main_ch then others
    x = zeros([1,ncharge]);                                   % preallocate index mapping
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));          % map from old order → new order
    end;
    tune = 1:npep;                                            % apply reorder to all peptides
    His.pep_mz(tune,:) = His.pep_mz(tune,x);                  % reorder m/z columns
    His.pep_ch(tune,:) = His.pep_ch(tune,x);                  % reorder charge columns
end;

% ============================================================================
% Local function: compute RTs/intensities and monoisotopic traces
% ============================================================================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%
% === Purpose ================================================================
% For each peptide and charge state:
%   - Determine representative RTs (pep_rts).
%   - Integrate intensities (pep_intens).
%   - Build monoisotopic intensity traces across RT (mono_isointens).
% Also: recalibrate RT anchors if not in debug mode.

[npep,ncharge] = size(His.pep_mz);                             % dimensions based on m/z table
num_MS1 = size(MS1_index,1);                                   % number of MS1 scans
pep_rts = zeros([npep,ncharge]);                               % preallocate RTs (peptide × charge)
pep_intens = zeros([npep,ncharge]);                            % preallocate intensities
mono_isointens = zeros([num_MS1,npep]);                        % monoisotopic trace per peptide (over RT)

% calibrate the rt_ref
His.rt_unmod_orig = His.rt_ref(1);                             % store original RT for first peptide
if 1~=special.ndebug                                           % allow RT anchor recalibration unless debug
    if 2~=special.nDAmode                                      % standard mode (nDAmode ~= 2)
        for hno=1:4                                            % loop over 4 peptides
            [His.rt_ref(hno),special.ndebug] = ...
                check_ref(special.raw_path, ...
                          [His.mod_short{hno},His.mod_type{hno}], ...
                          His.rt_ref(hno), ...
                          special.ndebug);                     % confirm/adjust reference RT anchors
        end;
    else                                                       % DA mode (nDAmode == 2): flexible refinement
        nhmass = special.nhmass;                               % DA helper mass for get_rts2 (if applicable)
        for hno=1:4
            rt_unmod_orig = His.rt_ref(hno);                   % keep current anchor
            His.rt_ref(hno) = check_ref(special.raw_path, ...
                                         [His.mod_short{hno},His.mod_type{hno}], ...
                                         His.rt_ref(hno), ...
                                         special.ndebug);      % try confirming anchor
            if rt_unmod_orig==His.rt_ref(hno)                  % if unchanged (no clear evidence)…
                t1 = 0;                                        % …search full RT range: start at 0
                t2 = MS1_index(num_MS1,2);                     % …end at last RT from MS1 index
            else                                               % if updated anchor found…
                delta = 5;                                     % …focus search in ±5 min around anchor
                t1 = His.rt_ref(hno)-delta;
                t2 = His.rt_ref(hno)+delta;
            end;
            [rts1,top1_rt1] = ...
                get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                         ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok % refine top RT candidate
            if 0==isempty(top1_rt1)                            % if a best RT was identified…
                His.rt_ref(hno) = top1_rt1;                    % …update anchor to that RT
            end;
        end;
        special.ndebug = 1;                                    % lock recalibrated anchors (no further changes)
    end;
end;

% 64-69KLPFQR
% 129-134RIRGER
% 130-135IRGERA
% (Notes above: historical fragment indexing for H3; left as original comments.)

for hno=1:4                                                    % iterate over all 4 peptides
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special); % extract per-peptide data
    if cur_rts(1)>0                                            % guard: only keep if a valid RT exists
        pep_rts(hno,1:ncharge) = cur_rts;                      % store RTs for all charges
        pep_intens(hno,1:ncharge) = cur_intens;                % store intensities for all charges
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;    % store monoisotopic RT trace
    end;
end;
