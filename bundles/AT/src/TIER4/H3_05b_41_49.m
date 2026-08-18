function H3_05b_41_49(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% ============================ DRAFT — NOT REGISTERED ============================
% H3.1 sibling of H3_05_41_49 (which encodes the H3.3 form YRPGTVALR).
%
% AT H3.1 carries F41 (UniProt P59226, genes HTR1/2/3/9/13) -> peptide FRPGTVALR.
% Position 41 is the canonical plant H3.1/H3.3 discriminator (F in H3.1, Y in
% H3.3). F41 has NO PTM site (no Tyr-OH, no Lys), so this peptide is unmod-only
% and serves only as a normaliser / to capture the H3.1 pool at 41-49.
%
% This module is a scaffold. Before enabling it:
%   1. Confirm with a PSM check that FRPGTVALR is actually present in the data.
%   2. Calibrate His.rt_ref (the value below is a placeholder copied from the
%      Y41 unmod RT; F is slightly more hydrophobic than Y so expect a small
%      shift).
%   3. Register it in TIER2/init_histone0.m by adding, after the H3_05 block:
%         no = no + 1;
%         His.out_filename{no,1} = 'H3_05b_41_49';
%         His.pep_seq{no,1}      = 'FRPGTVALR';
%         His.mod_type{no,1}     = '0,pr;';
%         His.pep_ch(no,1)       = 2;
%         His.pep_mz(no,1)       = calculate_pepmz0(His,no,special);
%         new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
%         His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));
% ===============================================================================

% check
out_filename = 'H3_05b_41_49';
fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end;

% init
His = init_histone();

% calculate
unitdiff = 1.0032;
[pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special);

% output
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% draw
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS2_index,MS2_peaks,special);

% Get PSM
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

function His = init_histone()
%%

His.pep_seq = 'FRPGTVALR';% AT H3.1: F41 (vs H3.3 YRPGTVALR). UniProt P59226
His.mod_short = {'unmod'};
His.mod_type = {'0,pr;'};% N-term propionyl only; F41 carries no PTM site

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.rt_ref = 31.2;% PLACEHOLDER (copied from Y41 unmod) — CALIBRATE before use
His.display = ones(length(His.mod_type),1);

% main ch
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end;
    tune = 1:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end;

function [pep_rts,pep_intens,mono_isointens] = calculate_layout(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,special)
%%

[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% unmod (only entry — FRPGTVALR has no PTM site)
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        [His.rt_ref(1),special.ndebug] = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
    else
        nhmass = special.nhmass;
        His.rt_ref(1) = check_ref(special.raw_path,[His.pep_seq,His.mod_type{1}],His.rt_ref(1),special.ndebug);
        if His.rt_unmod_orig==His.rt_ref(1)
            t1 = 0;
            t2 = MS1_index(num_MS1,2);
        else
            delta = 5;
            t1 = His.rt_ref(1)-delta;
            t2 = His.rt_ref(1)+delta;
        end;
        hno = 1;% unmod
        [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,1,t1,t2,nhmass);%#ok
        if 0==isempty(top1_rt1)
            His.rt_ref(1) = top1_rt1;
        end;
    end;
end;

hno = 1;
[cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);

% calibrate the rt_ref
if cur_rts(1)>0
    His.rt_ref(1) = cur_rts(1);
    pep_rts(hno,1:ncharge) = cur_rts;
    pep_intens(hno,1:ncharge) = cur_intens;
    mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
end;
