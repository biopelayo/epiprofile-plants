function H4_06_79_92(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ========================================================================
% H4_06_79_92 — Targeted panel for histone H4 peptide "KTVTAMDVVYALKR" (aa 79–92)
%               States: unmod, M84ox, Y88pr, M84oxY88pr
% ------------------------------------------------------------------------
% RATIONALE
% - Derivatization baseline ("unmod"): peptide N-term and Lys residues are
%   propionylated to block primary amines and simplify charge states.
% - Oxidized methionine (M84ox): +15.9949 Da; typically increases polarity,
%   often elutes earlier than the unmodified species.
% - Propionylated tyrosine (Y88pr): +56.0262 Da on the phenolic ring (can
%   occur under certain conditions), typically increases hydrophobicity and
%   elutes later than unmod.
% - Double mod (M84oxY88pr): combined effect; RT often between the two.
%
% INPUTS
%   MS1_index, MS1_peaks : centroided MS1 index and peaks (project schema)
%   MS2_index, MS2_peaks : centroided MS2 (for DA anchoring/PSM)
%   ptol                 : precursor tolerance in ppm
%   cur_outpath          : output folder for .mat artifacts and figures
%   special              : struct with fields:
%                          - ndebug   : 1 -> diagnostic relocateD path
%                          - nDAmode  : 0 MS1-only; 2 MS2-assisted; 1 also exports PSM
%                          - nhmass   : neutral mass helper for get_rts2 (if needed)
%                          - raw_path : raw file path for check_ref
%
% OUTPUT FILES (side effects)
%   <cur_outpath>/H4_06_79_92.mat with His/pep_rts/pep_intens/mono_isointens
%   layout/XIC figures via draw_layout
%   optional PSM dump when nDAmode==1 via GetPSM
% ========================================================================

% ------------------------------
% 1) Idempotency: skip if this panel is already processed
% ------------------------------
out_filename = 'H4_06_79_92';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end

% ------------------------------
% 2) Define peptide, PTM grid, charges, mz, initial RT seeds
% ------------------------------
His = init_histone();

% ------------------------------
% 3) Anchor/relocate/extract across states
% ------------------------------
unitdiff = 1.0032; % ~C13 spacing used by monoisotopic helpers
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);

% ------------------------------
% 4) Persist (matrix + metadata)
% ------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% ------------------------------
% 5) Plot composite layout/XIC for visual QC
% ------------------------------
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2); % RT vector (minutes)
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------
% 6) Optional: export PSM if requested
% ------------------------------
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % H4_06_79_92


% ========================================================================
% init_histone — define sequence, PTM states, charge grid, theoretical m/z,
%                 RT seeds and display flags
% ========================================================================
function His = init_histone()
%%
% Peptide sequence (H4 aa 79–92)
His.pep_seq = 'KTVTAMDVVYALKR';

% Rows = states to be quantified
His.mod_short = {'unmod';
    'M84ox';
    'Y88pr';
    'M84oxY88pr'};

% PTM encoding (peptide-relative indices; 0=N-term)
% - unmod        : 0,pr; 1,pr; 13,pr   (N-term, Lys1, Lys13 propionylated)
% - M84ox        : ...; 6,ox; ...      (Met at peptide pos 6 is oxidized)
% - Y88pr        : ...; 10,pr; ...     (Tyr at peptide pos 10 has propionyl)
% - M84oxY88pr   : both 6,ox and 10,pr
His.mod_type = {'0,pr;1,pr;13,pr;';
    '0,pr;1,pr;6,ox;13,pr;';
    '0,pr;1,pr;10,pr;13,pr;';
    '0,pr;1,pr;6,ox;10,pr;13,pr;'};

% Charge grid to consider; utilities will reorder so that the main charge
% (2+) appears first for alignment consistency across matrices.
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);

%{
% Legacy hard-coded m/z values (kept as historical reference only).
His.pep_mz = [1762.9772  881.9922  588.3306  441.4997
              1778.9721  889.9897  593.6622  445.4985
              1819.0034  910.0053  607.0060  455.5063
              1834.9983  918.0028  612.3376  459.5050];
%}

% Compute theoretical m/z for each row x charge based on seq + PTMs
His.pep_mz = calculate_pepmz(His);

% Seed RTs (minutes): [unmod; M84ox; Y88pr; M84oxY88pr]
His.rt_ref = [48.5
    46.8
    49.7
    47.7];

% Visibility flags for plotting (1=include). Here we keep unmod visible.
His.display = zeros(length(His.mod_type),1);
His.display(1) = 1;

% Ensure the "main" charge (2+) sits in the first column in pep_mz/pep_ch
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
% calculate_layout — anchor unmod, propagate drift, relocate windows for
%                    M84ox/Y88pr/double, extract each state
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

% ---- 1) Anchor the "unmod" row using check_ref (MS1 or MS2-assisted)
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only anchoring: refine RT seed by scanning the raw metadata
        [His.rt_ref(1),special.ndebug] = ...
            check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                      His.rt_ref(1),special.ndebug);
    else
        % DA anchoring with a local window around the check_ref suggestion
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

% ---- 2) Extract "unmod" and propagate any drift to the rest of seeds
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% ---- 3) Relocate search windows for modified states
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His); % diagnostic path
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His); % MS1-only
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                        ptol,unitdiff,His,nhmass); % MS2-assisted
    end
end

% ---- 4) Extract each modified row (M84ox, Y88pr, double)
for hno=2:4
    [cur_rts,cur_intens,cur_mono_isointens] = ...
        get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge)        = cur_rts;
        pep_intens(hno,1:ncharge)     = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

end % calculate_layout


% ========================================================================
% relocate — MS1-only relocation rules per state (windows relative to unmod)
% ========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
delta  = 0.1; % small exclusion to avoid overlapping the unmod apex
nsplit = 1;

% M84ox — typically earlier than unmod (more polar)
hno = 2;
t1 = His.rt_ref(1)-15;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% Y88pr — typically later than unmod (more hydrophobic)
hno = 3;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+20;
[rts3,top1_rt3] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end

% M84oxY88pr — around unmod (combined effects); search a tight window
hno = 4;
t1 = His.rt_ref(1)-10;
t2 = His.rt_ref(1)+5;
[rts4,top1_rt4] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt4;
end

end % relocate


% ========================================================================
% relocate2 — MS2-assisted relocation (same windows, uses get_rts2)
% ========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
delta  = 0.1;
nsplit = 1;

% M84ox
hno = 2;
t1 = His.rt_ref(1)-15;
t2 = His.rt_ref(1)-delta;
[rts2,top1_rt2] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

% Y88pr
hno = 3;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+20;
[rts3,top1_rt3] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts3)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt3;
end

% M84oxY88pr
hno = 4;
t1 = His.rt_ref(1)-10;
t2 = His.rt_ref(1)+5;
[rts4,top1_rt4] = ...
    get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
             ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

if 1==isempty(rts4)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt4;
end

end % relocate2
