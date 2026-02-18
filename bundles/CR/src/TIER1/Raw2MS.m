function Raw2MS(raw_path,raw_names)

raw2ms1_file = fullfile(pwd,'RawToMS1.exe');
if 0==exist(raw2ms1_file,'file')
    fprintf(1,'please first put RawToMS1.exe under the source code folder!\n');
    return;
end;

xtract_file = fullfile(pwd,'xtract.exe');
if 0==exist(xtract_file,'file')
    fprintf(1,'please first put xtract.exe under the source code folder!\n');
    return;
end;

ms1_path = fullfile(raw_path,'MS1');
if 0==exist(ms1_path,'dir') && 0==mkdir(ms1_path)
    fprintf(1,'can not create: %s\n',ms1_path);
    return;
end;

ms2_path = fullfile(raw_path,'MS2');
if 0==exist(ms2_path,'dir') && 0==mkdir(ms2_path)
    fprintf(1,'can not create: %s\n',ms2_path);
    return;
end;

for i=1:length(raw_names)
    raw_file = fullfile(raw_path,[raw_names{i},'.RAW']);
    ms1_file = fullfile(ms1_path,[raw_names{i},'.MS1']);
    ms2_file = fullfile(ms2_path,[raw_names{i},'.ms2']);
    if 0==exist(ms1_file,'file')
        command = [raw2ms1_file,' ',['"',raw_file,'"']];
        if 0~=system(command)
            fprintf(1,'can not convert: %s\n',raw_file);
            return;
        end;
        movefile(fullfile(raw_path,[raw_names{i},'.MS1']),ms1_path);
    end;
    if 0==exist(ms2_file,'file')
        command = [xtract_file,' -a -i 1 -m 5 -ms2 -o ',['"',ms2_path,'"'],' ',['"',raw_file,'"']];
        [status,message] = system(command);%#ok
        if 0~=status
            fprintf(1,'can not convert: %s\n',raw_file);
            return;
        end;
    end;
end;