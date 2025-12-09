function H4_05_68_78(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ========================================================================
% H4_05_68_78 — Targeted panel for histone H4 peptide "DAVTYTEHAKR" (aa 68–78)
%               States: unmodified (derivatization baseline) and "Y72pr".
% ------------------------------------------------------------------------
% Design notes:
% - "unmod": peptide N-term and internal Lys are propionylated (usual baseline
%   after propionic anhydride derivatization to block primary amines).
% - "Y72pr": additional propionyl group on Tyr (peptide index 5). This can
%   appear under certain derivatization conditions and typically increases
%   hydrophobicity -> later retention than "unmod".
%
% Inputs:
%   MS1_index, MS1_peaks : centroided MS1 index and peaks (project schema)
%   MS2_index, MS2_peaks : centroided MS2 (optional but recommended for DA)
%   ptol                 : precursor tolerance in ppm
%   cur_outpath          : folder to write .mat and figures
%   special              : struct with fields:
%                          - ndebug   : 1 -> diagnostic relocateD path
%                          - nDAmode  : 0 MS1-only; 2 MS2-assisted; 1 also PSM export
%                          - nhmass   : neutral mass helper for get_rts2 (if used)
%                          - raw_path : raw file path used by check_ref
%
% Side-effects (files):
%   <cur_outpath>/H4_05_68_78.mat with His, pep_rts, pep_intens, mono_isointens
%   XIC/layout figures via draw_layout; optional PSM via GetPSM when nDAmode==1.
% ========================================================================

% ------------------------------
% 1) Skip if results already exist (idempotent behavior)
% ------------------------------
out_filename = 'H4_05_68_78';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end

% ------------------------------
% 2) Initialize peptide/PTMs/charges/mz/RT seeds
% ------------------------------
His = init_histone();

% ------------------------------
% 3) Main quantification: anchor, relocate, extract
% ------------------------------
unitdiff = 1.0032; % ~C13 spacing used by monoisotopic helpers
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);

% ------------------------------
% 4) Persist quantification to disk
% ------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% ------------------------------
% 5) Plot layout/XIC for visual QC
% ------------------------------
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2); % RT vector (minutes)
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------
% 6) Optional PSM export (if requested)
% ------------------------------
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % H4_05_68_78


% ========================================================================
% init_histone — define sequence, PTM rows, charges, theoretical m/z, RT
% ========================================================================
function His = init_histone()
%%
% Peptide sequence (H4 aa 68–78)
His.pep_seq = 'DAVTYTEHAKR';

% Two rows (states):
% - 'unmod' : derivatization baseline (N-term + Lys propionylated)
% - 'Y72pr' : additional propionyl on Tyr (peptide index 5)
His.mod_short = {'unmod';
    'Y72pr'};

% PTM encoding (peptide-relative indices):
%   0,pr   -> N-terminal propionylation
%   10,pr  -> Lys propionylation (the K in this peptide)
%   5,pr   -> Tyr propionylation (non-canonical but observed under some conditions)
His.mod_type = {'0,pr;10,pr;';
    '0,pr;5,pr;10,pr;'};

% Charge grid considered for extraction (1+,2+,3+,4+). Utilities will reorder
% so that the main charge (2+) becomes the first column for alignment.
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);

%{
% Legacy hard-coded m/z values (kept for debugging reference); now we compute
% dynamically from sequence + PTM specification using calculate_pepmz(His).
His.pep_mz = [1402.6961  701.8517  468.2369  351.4295
              1458.7223  729.8648  486.9123  365.4360];
%}

% Compute theoretical m/z for each (row x charge) from sequence + PTMs
His.pep_mz = calculate_pepmz(His);

% Seed RTs (minutes). Expectation: Y72pr retains longer than unmod.
His.rt_ref = [31.08
    37.08];

% Display flags (1=include in composite plots)
His.display = zeros(length(His.mod_type),1);

% Ensure the "main" charge sits in the first column (cosmetic consistency)
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


% ========================================================================
% calculate_layout — anchor baseline, drift-correct, relocate Y72pr, extract
% ========================================================================
function [pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special)
%%
[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);

% Preallocate outputs
pep_rts        = zeros([npep,ncharge]);
pep_intens     = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% ---- 1) Anchor the "unmod" row using check_ref (MS1-only or DA)
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only anchoring: refine RT seed by scanning raw path metadata
        [His.rt_ref(1),special.ndebug] = ...
            check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                      His.rt_ref(1),special.ndebug);
    else
        % MS2-assisted anchoring with a local window around check_ref
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                                  His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0; t2 = MS1_index(num_MS1,2);  % full gradient fallback
        else
            delta = 5; t1 = His.rt_ref(1)-delta; t2 = His.rt_ref(1)+delta;
        end
        hno = 1; % unmod
        [rts1,top1_rt1] = ...
            get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,hno,1,t1,t2,nhmass); %#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end
    end
end

% ---- 2) Extract unmod to compute delta and correct remaining seeds
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

if cur_rts(1)>0
    % Update anchor and propagate drift to other seeds
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    % Write back extracted metrics for the anchored row
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% ---- 3) Relocate Y72pr search window
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His); % diagnostic path
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His); % MS1-only window
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                        ptol,unitdiff,His,nhmass); % DA-assisted window
    end
end

% ---- 4) Extract Y72pr (row 2) per charge
hno = 2;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if cur_rts(1)>0
    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

end % calculate_layout


% ========================================================================
% relocate — MS1-only relocation for Y72pr (later than unmod)
% ========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
delta  = 0.1; % small exclusion to avoid overlapping anchor apex
nsplit = 1;

% Y72pr: expected to elute AFTER unmod (more hydrophobic)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts2,top1_rt2] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;     % not found
else
    His.rt_ref(hno) = top1_rt2; % adopt the top-peak RT
end

end % relocate


% ========================================================================
% relocate2 — MS2-assisted relocation for Y72pr (same window, get_rts2)
% ========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
delta  = 0.1;
nsplit = 1;

hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts2,top1_rt2] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

end % relocate2
