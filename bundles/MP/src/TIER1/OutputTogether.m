function OutputTogether(layout_path,raw_names)
%%

% npeps
npeps = 0;
i = 1;
cur_rawname = raw_names{i};
if i<10
    prefix = ['0',num2str(i)];
else
    prefix = num2str(i);
end;
cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1);
    npeps = npeps + size(His.pep_mz,1);
end;

% info
info = zeros([npeps,length(raw_names)*3]);
for i=1:length(raw_names)
    cur_rawname = raw_names{i};
    if i<10
        prefix = ['0',num2str(i)];
    else
        prefix = num2str(i);
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');
    m = 0;
    for j = 1:length(matfiles)
        matfile1 = fullfile(cur_outpath,matfiles(j).name);
        load(matfile1);
        nlen = size(His.pep_mz,1);
        info(m+1:m+nlen,(i-1)*3+1:(i-1)*3+3) = auc;
        m = m+nlen;
    end;
end;

% output
[raw_path,fname] = fileparts(layout_path);
if 1==strcmp(fname,'histone_layouts')
    inten_file = fullfile(raw_path,'histone_ratios.xls');
else
    inten_file = fullfile(fileparts(raw_path),['histone_ratios_',fname,'.xls']);
end;
fp = fopen(inten_file,'w');
if -1==fp
    disp(['can not open: ',inten_file]);
    return;
end;

%{
for i=1:length(raw_names)
    fprintf(fp,'\t%d,%s\t\t',i,raw_names{i});
end;
fprintf(fp,'\r\n');

fprintf(fp,'Peptide');
for i=1:length(raw_names)
    fprintf(fp,'\tRT(min)\tArea\tRatio');
end;
fprintf(fp,'\r\n');
%}
for i=1:length(raw_names)
    fprintf(fp,'\t%d,%s',i,raw_names{i});
end;
fprintf(fp,'\t');
for i=1:length(raw_names)
    fprintf(fp,'\t%d,%s',i,raw_names{i});
end;
fprintf(fp,'\t');
for i=1:length(raw_names)
    fprintf(fp,'\t%d,%s',i,raw_names{i});
end;
fprintf(fp,'\r\n');

fprintf(fp,'Peptide');
for i=1:length(raw_names)
    fprintf(fp,'\tRatio');
end;
fprintf(fp,'\t');
for i=1:length(raw_names)
    fprintf(fp,'\tArea');
end;
fprintf(fp,'\t');
for i=1:length(raw_names)
    fprintf(fp,'\tRT(min)');
end;
fprintf(fp,'\r\n');

i = 1;
cur_rawname = raw_names{i};
if i<10
    prefix = ['0',num2str(i)];
else
    prefix = num2str(i);
end;
cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');

m = 0;
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1);
    pepname = matfiles(j).name(1:end-4);
    if 1==strcmp(pepname(1:2),'HH')
        pepname = pepname(2:end);
    end;
    p = strfind(pepname,'_');
    version = pepname(p(1)+1:p(2)-1);
    
    if 1==strcmp(fname,'SILAC') || 1==strcmp(fname,'C13')
        if 1==strcmp(pepname(1:7),'H3_04v3')
            pep_seq = ['H33',His.pep_seq(3:end)];
        else
            pep_seq = His.pep_seq;
        end;
        fprintf(fp,'%s\r\n',pep_seq);
    else
        if length(version)>=4
            if 1==strcmp(pepname(1:7),'H3_04v3')
                version = version(4);
            else
                version = version(4:end);
            end;
            histone_pos = [pepname(1:p(1)-1),version,'_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
        else
            histone_pos = [pepname(1:p(1)-1),'_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
        end;
        fprintf(fp,'%s(%s)\r\n',His.pep_seq,histone_pos);
    end;
    
    nlen = size(His.pep_mz,1);
    
    for ino = m+1:m+nlen
        %{
        fprintf(fp,'%s',His.mod_short{ino-m});
        for jno=1:length(raw_names)
            fprintf(fp,'\t%.2f\t%e\t%f',info(ino,(jno-1)*3+1),info(ino,(jno-1)*3+2),info(ino,(jno-1)*3+3));
        end;
        fprintf(fp,'\r\n');
        %}
        if 1==strcmp(fname,'SILAC') || 1==strcmp(fname,'C13')
            fprintf(fp,'%s',His.mod_short{ino-m});
        else
            fprintf(fp,'%s %s',histone_pos,His.mod_short{ino-m});
        end;
        for jno=1:length(raw_names)
            fprintf(fp,'\t%f',info(ino,(jno-1)*3+3));
        end;
        fprintf(fp,'\t');
        for jno=1:length(raw_names)
            fprintf(fp,'\t%e',info(ino,(jno-1)*3+2));
        end;
        fprintf(fp,'\t');
        for jno=1:length(raw_names)
            fprintf(fp,'\t%.2f',info(ino,(jno-1)*3+1));
        end;
        fprintf(fp,'\r\n');
    end;
    m = m+nlen;
end;
fclose(fp);