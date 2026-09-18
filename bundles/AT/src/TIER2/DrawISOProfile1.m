function DrawISOProfile1(raw_path,raw_names,ptol,special)
%%

layout_path = fullfile(raw_path,'histone_layouts');
if 0==exist(layout_path,'dir') && 0==mkdir(layout_path)
    fprintf(1,'can not create: %s\n',layout_path);
    return;
end;

for i=1:length(raw_names)
    if i<10
        prefix = ['0',num2str(i)];
    else
        prefix = num2str(i);
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',raw_names{i}]);
    if 0==exist(cur_outpath,'dir') && 0==mkdir(cur_outpath)
        fprintf(1,'can not create: %s\n',cur_outpath);
        return;
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',raw_names{i}],'detail');
    if 0==exist(cur_outpath,'dir') && 0==mkdir(cur_outpath)
        fprintf(1,'can not create: %s\n',cur_outpath);
        return;
    end;
end;

for i=1:length(raw_names)
    fprintf(1,'\n%s\n',raw_names{i});
    cur_rawname = raw_names{i};
    MS1_scanfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1scans.mat']);
    MS1_peakfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1peaks.mat']);
    MS2_scanfile = fullfile(raw_path,'MS2',[cur_rawname,'_MS2scans.mat']);
    MS2_peakfile = fullfile(raw_path,'MS2',[cur_rawname,'_MS2peaks.mat']);
    load(MS1_scanfile);% MS1_index
    load(MS1_peakfile);% MS1_peaks
    load(MS2_scanfile);% MS2_index
    load(MS2_peakfile);% MS2_peaks
    nlen = length(unique(MS2_index(:,4)));%#ok
    c0 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+0,4) && MS2_index(2,4)==MS2_index(2+nlen+0,4);
    c1 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+1,4) && MS2_index(2,4)==MS2_index(2+nlen+1,4);
    c2 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+2,4) && MS2_index(2,4)==MS2_index(2+nlen+2,4);
    c3 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+3,4) && MS2_index(2,4)==MS2_index(2+nlen+3,4);
    c4 = size(MS2_index,1)>2+nlen+4 && MS2_index(1,4)==MS2_index(1+nlen+4,4) && MS2_index(2,4)==MS2_index(2+nlen+4,4);
    if nlen<270 && (c0 || c1 || c2 || c3 || c4)% (1100-300)/3=267
        special.nDAmode = 2;% DIA
    else
        special.nDAmode = 1;% DDA
    end;
    if i<10
        prefix = ['0',num2str(i)];
    else
        prefix = num2str(i);
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');
    if 1==special.nsource || 2==special.nsource
        % histone_normal or histone_SILAC
        if '1'==special.soutput(1)
            % 1, H3H4 basic
            H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_06_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_11_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_12_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_13_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_14_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_16_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_17_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_18_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
            
            H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_04_40_45(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
        elseif '2'==special.soutput(1)
            % 2, H3H4 basic + H3S10ph
            H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02a_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_06_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
            recal_ratio_H3_0201(cur_outpath);
            
            H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_04_40_45(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
        else
            % 3, H3H4 all
            H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02a_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02b_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04a_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %H3_04v3a_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% Tier 4: mammalian peptide, removed
            H3_05_41_49(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_06_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_06a_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_09u_64_135(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
            recal_ratio_H3_0202(cur_outpath);
            recal_ratio_H3_0402(cur_outpath);
            
            H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02a_18_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02b_20_35(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02c_20_36(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %H4_03_24_35(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% file never existed in AT ZIP
            H4_04_40_45(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_05_68_78(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_06_79_92(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_07u_24_102(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
        end;
        H3_Snapshot(cur_outpath);
        H4_Snapshot(cur_outpath);
        HH2A_AT_Snapshot(cur_outpath);
        HH2B_AT_Snapshot(cur_outpath);
        HH1_AT_Snapshot(cur_outpath);
        
        if '1'==special.soutput(2)
            % H2A AT
            HH2A_01u_1_7(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_AT01_can_20_31(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_AT02_can1_5_19(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_AT03_Z_33_41(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_AT04_X_15_25(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_AT05_Xa_41_50(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_AT06_W12_16_27(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_AT07_W12_28_39(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2A_AT08_W6_1_14(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            % [Retirado 2026-05-24] HH2A_AT09_can_36_42 (KGKYAER) - peptido fantasma.
            HH2A_06m1_72_77(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
            % H2B AT
            HH2B_AT01_shared_94_100(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2B_AT02_shared_108_117(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2B_AT03_shared_118_134(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2B_AT04_67_disc(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2B_AT05_14_disc(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2B_AT06_10_disc(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH2B_AT07_311_disc(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            % [Retirado 2026-05-24] HH2B_AT08_shared_80_86 (LARYNKK) - peptido fantasma.
            HH2BMo_03u_1_100(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
            % H1 AT
            HH1_AT01_H11_104_112(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH1_AT02_H11_80_95(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH1_AT03_H12_104_112(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            HH1_AT04_H12_80_103(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
        end;
        
        if 1==special.nDAmode
            GetBenchmark(cur_outpath);
        end;
    elseif 3==special.nsource
        % histone_C13
        H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        %HH2A_02m1_4_11(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% removed: human-only module
        %HH2A_04oV_1_19(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% removed: human-only module
        %HH2A_04oZ_1_19(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% removed: human-only module
        %HH2A_05m1_12_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% removed: human-only module
        fprintf(1,'\n');
    elseif 4==special.nsource
        % histone_N15 — disabled for AT bundle (human-only H2A modules removed)
        if 0==special.nsubtype || 4==special.nsubtype || 5==special.nsubtype
            H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %HH2A_02m1_4_11(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% removed: human-only
            %HH2A_05m3_12_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% removed: human-only
            fprintf(1,'\n');
        else
            %H3N_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);% N15 modules not in AT bundle
            %H3N_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %H3N_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %H3N_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %H3N_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %H4N_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %H4N_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %HH2AN_02m1_4_11(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            %HH2AN_05m3_12_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
            fprintf(1,'\n');
        end;
    elseif 5==special.nsource
        % histone_13CD3
        H3_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_06_53_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H3_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        fprintf(1,'\n');
        
        H4_01_4_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H4_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        H4_06_79_92(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
        fprintf(1,'\n');
    end;
end;

OutputTogether(layout_path,raw_names);
if 1==special.nsource || 2==special.nsource
    OutputSinglePTMs(layout_path,raw_names);
end;

if (1==special.nsource || 2==special.nsource) && 1==special.nfigure && length(raw_names)>2
    OutputFigures(layout_path,raw_names);
end;

function recal_ratio_H3_0201(cur_outpath)
%%

% H3_02
out_filename = 'H3_02_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s1 = sum(auc(:,2));%#ok

out_filename = 'H3_02a_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s2 = sum(auc(:,2));

s = s1+s2;

out_filename = 'H3_02_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

out_filename = 'H3_02a_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;%#ok
save(mat_file,'His','auc');

function recal_ratio_H3_0202(cur_outpath)
%%

% H3_02
out_filename = 'H3_02_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s1 = sum(auc(:,2));%#ok

out_filename = 'H3_02a_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s2 = sum(auc(:,2));

out_filename = 'H3_02b_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
s3 = sum(auc(:,2));

s = s1+s2+s3;

out_filename = 'H3_02_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

out_filename = 'H3_02a_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;
save(mat_file,'His','auc');

out_filename = 'H3_02b_9_17';
mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
load(mat_file);
auc(:,3) = auc(:,2)/s;%#ok
save(mat_file,'His','auc');

function recal_ratio_H3_0402(cur_outpath)
%%

% H3_04 + H3_04a (only if both exist)
mat_file_04 = fullfile(cur_outpath,'H3_04_27_40.mat');
mat_file_04a = fullfile(cur_outpath,'H3_04a_27_40.mat');
if 0~=exist(mat_file_04,'file') && 0~=exist(mat_file_04a,'file')
    load(mat_file_04);
    s1 = sum(auc(:,2));%#ok
    load(mat_file_04a);
    s2 = sum(auc(:,2));
    s = s1+s2;
    load(mat_file_04);
    auc(:,3) = auc(:,2)/s;
    save(mat_file_04,'His','auc');
    load(mat_file_04a);
    auc(:,3) = auc(:,2)/s;
    save(mat_file_04a,'His','auc');
end;

% H3_04v3 + H3_04v3a (only if H3_04v3a exists -- removed for AT, Tier 4)
out_filename = 'H3_04v3a_27_40';
mat_file_v3a = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(mat_file_v3a,'file')
    out_filename = 'H3_04v3_27_40';
    mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
    load(mat_file);
    s1 = sum(auc(:,2));
    load(mat_file_v3a);
    s2 = sum(auc(:,2));
    s = s1+s2;
    out_filename = 'H3_04v3_27_40';
    mat_file = fullfile(cur_outpath,[out_filename,'.mat']);
    load(mat_file);
    auc(:,3) = auc(:,2)/s;
    save(mat_file,'His','auc');
    load(mat_file_v3a);
    auc(:,3) = auc(:,2)/s;%#ok
    save(mat_file_v3a,'His','auc');
end;% end of H3_04v3a check