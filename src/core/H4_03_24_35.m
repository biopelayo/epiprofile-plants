function H4_03_24_35(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ========================================================================
% H4_03_24_35 — Targeted panel for histone H4 peptide "DNIQGITKPAIR" (aa 24–35)
%               Monitoring: unmodified (derivatized baseline) and T30pr.
% ------------------------------------------------------------------------
% Why this panel?
% - Short H4 peptide covering Thr30 and Lys31 region.
% - "T30pr" captures O-propionylation on Thr30 (method-dependent artifact),
%   useful as a QA/chemistry marker. "unmod" has N-term+Lys propionylation.
%
% Pipeline (consistent with the core across panels):
%   1) Panel definition (sequence, PTMs, charges, theoretical m/z, seed RTs)
%   2) Anchor the unmodified row via check_ref (optionally refined by MS2)
%   3) Extract unmodified (get_histone0) and correct global RT drift (delta)
%   4) Relocate window for T30pr (MS1 relocate / DA relocate2)
%   5) Extract T30pr (get_histone1), persist, draw, optional PSM export
%
% Inputs:
%   MS1_index, MS1_peaks       : centroided MS1 (index with RT in minutes)
%   MS2_index, MS2_peaks       : centroided MS2 (DA anchoring / PSM export)
%   ptol                       : mass tolerance (ppm)
%   cur_outpath                : output directory
%   special.ndebug             : 1 -> minimal diag (relocateD), else relocate
%   special.nDAmode            : 2 -> MS2-assisted relocation; 1 -> export PSM
%   special.nhmass             : neutral-loss / helper (if applicable)
%   special.raw_path           : raw path for check_ref anchor
%
% Outputs (to disk):
%   <cur_outpath>/H4_03_24_35.mat  with His, pep_rts, pep_intens, mono_isointens
%   Layout figures (draw_layout) and optional PSM files if nDAmode==1
% ========================================================================

% ------------------------------
% 1) Idempotent check
% ------------------------------
out_filename = 'H4_03_24_35';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return; % respect original short-circuit if results already exist
end

% ------------------------------
% 2) Define panel (sequence/PTMs/charges/mz/RT seeds)
% ------------------------------
His = init_histone();

% ------------------------------
% 3) Compute layout (anchor, relocate, extract)
% ------------------------------
unitdiff = 1.0032; % ~C13 spacing used by monoisotopic helpers
[pep_rts,pep_intens,mono_isointens] = ...
    calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                     ptol,unitdiff,His,special);

% ------------------------------
% 4) Persist results
% ------------------------------
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% ------------------------------
% 5) QC plots
% ------------------------------
num_MS1 = size(MS1_index,1);
isorts  = MS1_index(1:num_MS1,2); % RT vector
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
            isorts,mono_isointens,MS2_index,MS2_peaks,special);

% ------------------------------
% 6) Optional PSM export (nDAmode==1)
% ------------------------------
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens, ...
           isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

end % H4_03_24_35


% ========================================================================
% init_histone — sequence, PTM rows, charge grid, theoretical m/z, RT seeds
% ========================================================================
function His = init_histone()
%%
% Core H4 peptide (aa 24–35)
His.pep_seq = 'DNIQGITKPAIR';

% PTM catalog:
% - "unmod": derivatized baseline (N-term + Lys propionylation)
% - "T30pr": additional O-propionylation at Thr30 (index 7 in this peptide)
His.mod_short = {'unmod';
    'T30pr'};

% Compact encoding (peptide-relative indices):
%   0,pr  -> N-terminal propionylation
%   8,pr  -> Lys propionylation (present in both rows)
%   7,pr  -> Thr O-propionylation (T30pr row only)
His.mod_type = {'0,pr;8,pr;';
    '0,pr;7,pr;8,pr;'};

% Expected charge states; "main" charge will be arranged at first column
His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);

%{
% If needed for debugging, historical fixed m/z (now superseded by calculate_pepmz):
His.pep_mz = [1437.8060  719.4066  479.9402  360.2070
              1493.8322  747.4197  498.6156  374.2135];
%}

% Compute theoretical m/z from sequence + PTM specification
His.pep_mz = calculate_pepmz(His);

% Seed retention times (minutes)
His.rt_ref = [42.45
    45.9];

% Auxiliary panel (hidden in composite plots by default)
His.display = zeros(length(His.mod_type),1);

% Cosmetic: set the "main" charge (here 2+) as first column for alignment
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
% calculate_layout — anchor unmod, drift-correct, relocate, extract T30pr
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

% ---- Anchor: unmodified (derivatized baseline)
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        % MS1-only anchoring
        [His.rt_ref(1),special.ndebug] = ...
            check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                      His.rt_ref(1),special.ndebug);
    else
        % DA anchoring: refine in ±5 min around check_ref or scan full gradient
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}], ...
                                  His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0; t2 = MS1_index(num_MS1,2);
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

% Extract unmodified (gives monoisotopic trace as well)
hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = ...
    get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% If an anchor RT was observed, apply global RT drift correction (delta)
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    delta = cur_rts(1)-His.rt_unmod_orig;
    His.rt_ref(2:end) = His.rt_ref(2:end) + delta;

    pep_rts(hno,1:ncharge)        = cur_rts;
    pep_intens(hno,1:ncharge)     = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end

% ---- Relocate the T30pr window
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His); % minimal diag
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His); % MS1-only
    else
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks, ...
                        ptol,unitdiff,His,nhmass); % DA path
    end
end

% ---- Extract T30pr (row 2)
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
% relocate — MS1-only relocation for T30pr
% ========================================================================
function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
delta  = 0.1;
nsplit = 1;

% T30pr expected after unmodified (derivatized baseline)
hno = 2;
t1 = His.rt_ref(1)+delta;
t2 = His.rt_ref(1)+30;
[rts2,top1_rt2] = ...
    get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

if 1==isempty(rts2)
    His.rt_ref(hno) = 0;
else
    His.rt_ref(hno) = top1_rt2;
end

end % relocate


% ========================================================================
% relocate2 — MS2-assisted relocation (identical window, get_rts2 backend)
% ========================================================================
function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
delta  = 0.1;
nsplit = 1;

% T30pr with DA support
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
