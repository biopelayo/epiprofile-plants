function HH2A_01u_1_7(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%% HH2A_01u_1_7 — Unmod-only H2A diagnostic peptides (AT)
%
%  TIER: T2 (EpiProfile-PLANTS, Arabidopsis thaliana)
%  HISTONE: H2A (variant-diagnostic short peptides)
%  REGION: Core domain (conserved across H2A canonical/Z/X/W subtypes)
%
%  PEPTIDES: 8 diagnostic peptides covering 13 AT H2A variants (HTA1-HTA13)
%    1. DNKKSR   — canonical H2A (HTA1, HTA2, HTA10, HTA13)
%    2. DNKKTR   — H2A.Z (HTA8, HTA9, HTA11)
%    3. DNKKNR   — H2A.X (HTA3, HTA5)
%    4. HLLLAIR  — canonical variant
%    5. HLQLAIR  — canonical/Z variant
%    6. HLCLAIR  — H2A.W (HTA6, HTA7, HTA12)
%    7. HIQLAVR  — H2A.Z variant
%    8. HVLLAVR  — divergent (HTA4)
%
%  NOTES:
%    - Unmod search only. K residues propionylated.
%    - rt_ref = 0 (to be calibrated on first run).

% check
out_filename = 'HH2A_01u_1_7';
%fprintf(1,'%s..',out_filename);
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
    GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end;

function His = init_histone()
%%

His.pep_seq = 'unmod';
His.mod_short = {'DNKKSR';
    'DNKKTR';
    'DNKKNR';
    'HLLLAIR';
    'HLQLAIR';
    'HLCLAIR';
    'HIQLAVR';
    'HVLLAVR'};
His.mod_type = {'0,pr;3,pr;4,pr;';
    '0,pr;3,pr;4,pr;';
    '0,pr;3,pr;4,pr;';
    '0,pr;';
    '0,pr;';
    '0,pr;';
    '0,pr;';
    '0,pr;'};

His.pep_ch = repmat([1 2 3 4],length(His.mod_type),1);
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [0;0;0;0;0;0;0;0];
His.display = zeros(length(His.mod_type),1);

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
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        for hno=1:npep
            [His.rt_ref(hno),special.ndebug] = check_ref(special.raw_path,[His.mod_short{hno},His.mod_type{hno}],His.rt_ref(hno),special.ndebug);
        end;
    else
        nhmass = special.nhmass;
        for hno=1:npep
            rt_unmod_orig = His.rt_ref(hno);
            His.rt_ref(hno) = check_ref(special.raw_path,[His.mod_short{hno},His.mod_type{hno}],His.rt_ref(hno),special.ndebug);
            if rt_unmod_orig==His.rt_ref(hno)
                t1 = 0;
                t2 = MS1_index(num_MS1,2);
            else
                delta = 5;
                t1 = His.rt_ref(hno)-delta;
                t2 = His.rt_ref(hno)+delta;
            end;
            [rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,1,t1,t2,nhmass);%#ok
            if 0==isempty(top1_rt1)
                His.rt_ref(hno) = top1_rt1;
            end;
        end;
        special.ndebug = 1;
    end;
end;

for hno=1:npep
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;
