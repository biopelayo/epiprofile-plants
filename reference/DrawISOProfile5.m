function DrawISOProfile5(raw_path,raw_names,ptol,special)
%%

special.onlyme = 1;
layout_path = fullfile(raw_path,'histone_layouts','13CD3');
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
    H3M_01_3_8(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3M_02_9_17(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3M_03_18_26(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3M_04_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3M_04v3_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3M_06_54_63(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3M_07_73_83(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H3M_08_117_128(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4M_02_20_23(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    H4M_06_79_92(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special);
    fprintf(1,'\n');
end;

OutputTogether(layout_path,raw_names);