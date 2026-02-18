function H4_02v_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%% H4_02v_20_23 — H4K20 quantification for Marchantia polymorpha
% TIER: T3 (new for MP bundle)
%
% Sequence: KVFR (MP-specific; human/AT = KVLR, L->F substitution at pos 21)
% PTMs: unmod, K20me1, K20me2, K20me3, K20ac
%
% This module quantifies both the MP-specific KVFR and the canonical KVLR
% to allow detection of either isoform in the sample.

% check
out_filename = 'H4_02v_20_23';
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

His.pep_seq = 'unmod';
His.mod_short = {'KVLR';
    'KVFR'};
His.mod_type = {'0,pr;1,pr;';
    '0,pr;1,pr;'};

His.pep_ch = repmat([1 2],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [29.52
    30];
His.display = ones(length(His.mod_type),1);

% main ch
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);%#ok
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
num_isowin = 4;
pep_rts = zeros(npep,ncharge);
pep_intens = zeros(npep,ncharge);
mono_isointens = zeros(npep,ncharge,num_isowin);

for ipep=1:npep
    for ich=1:ncharge
        cur_ch = His.pep_ch(ipep,ich);
        cur_mz = His.pep_mz(ipep,ich);

        [pair_MS1_index,pair_MS1_peaks] = find_pair(MS1_index,MS1_peaks,cur_mz,ptol);
        [pair_MS2_index,pair_MS2_peaks] = find_pair(MS2_index,MS2_peaks,cur_mz,ptol);

        iso_mzs = cur_mz + (0:num_isowin-1)*unitdiff/cur_ch;
        [cur_rt,cur_inten,cur_isointens] = get_histone0(pair_MS1_index,pair_MS1_peaks,pair_MS2_index,pair_MS2_peaks,iso_mzs,ptol,num_MS1,cur_ch,special);

        if ipep==1 && ich==1
            rt_unmod = cur_rt;
        end;

        pep_rts(ipep,ich) = cur_rt;
        pep_intens(ipep,ich) = cur_inten;
        mono_isointens(ipep,ich,:) = cur_isointens;
    end;
end;
