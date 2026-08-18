function HH2A_08u_4_99(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%

% check
out_filename = 'HH2A_08u_4_99';
% fprintf(1,'%s..',out_filename(2:end));
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
His.mod_short = {'HLQLAIR';
    'HLQLAVR';
    'GGKKKSTKTSR';
    'SGKKKMSKLSR';
    'IHRHLKTR';
    'IHRHLKSR';
    'NDEELNKLLGR';
    'AGLQFPVGR';
    'VHRLLR';
    'IHPELLAKKR';
    'YIKKGHPKYR'};
His.mod_type = {'0,pr;';
    '0,pr;';
    '0,pr;3,pr;4,pr;5,pr;8,pr;';
    '0,pr;3,pr;4,pr;5,pr;8,pr;';
    '0,pr;6,pr;';
    '0,pr;6,pr;'
    '0,pr;7,pr;';
    '0,pr;';
    '0,pr;';
    '0,pr;8,pr;9,pr;';
    '0,pr;3,pr;4,pr;8,pr;'};

His.pep_ch = repmat([1 2 3],length(His.mod_type),1);
%{
His.pep_mz = [906.5520	453.7796	302.8555
    892.5363	446.7718	298.1836
    1457.8322	729.4197	486.6156
    1529.8719	765.4396	510.6288
    1172.7011	586.8542	391.5719
    1158.6854	579.8464	386.9000
    1412.7380	706.8726	471.5842
    1000.5574	500.7824	334.1907
    849.5417	425.2745	283.8521
    1372.8311	686.9192	458.2819
    1513.8525	757.4299	505.2890];
%}
His.pep_mz = calculate_pepmz(His);
His.rt_ref = [34.47
    32.31
    30.92
    30
    38.84
    22.07
    44.64
    39.95
    27.04
    30
    30];
His.display = zeros(length(His.mod_type),1);

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

% calibrate the rt_ref
His.rt_unmod_orig = His.rt_ref(1);
if 1~=special.ndebug
    if 2~=special.nDAmode
        for hno=1:11
            [His.rt_ref(hno),special.ndebug] = check_ref(special.raw_path,[His.mod_short{hno},His.mod_type{hno}],His.rt_ref(hno),special.ndebug);
        end;
    else
        nhmass = special.nhmass;
        for hno=1:11
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

% 82-88HLQLAIR
% 82-88HLQLAVR
% 4-14GGKKKSTKTSR
% 4-14SGKKKMSKLSR
% 32-39IHRHLKTR
% 32-39IHRHLKSR
% 89-99NDEELNKLLGR
% 21-29AGLQFPVGR
% 30-35VHRLLR
% 108-117IHPELLAKKR
% 30-39IKKGHPKYR
for hno=1:11
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end;
end;