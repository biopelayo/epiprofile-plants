function OutputSinglePTMs(layout_path,raw_names)
%%

%% nos, peptides, and info
i = 1;
cur_rawname = raw_names{i};
if i<10
    prefix = ['0',num2str(i)];
else
    prefix = num2str(i);
end;
cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');

% separate nos
npeps = 0;
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
% Fix A-CRIT-004 / C-CRIT-04: orden estable
[~,idx_mf] = sort({matfiles.name});
matfiles = matfiles(idx_mf);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1);
    cnp = size(His.pep_mz,1);
    pepname = matfiles(j).name(1:end-4);
    if 1==strcmp(pepname(1:2),'HH')
        pepname = pepname(2:end);
    end;
    npeps = npeps + cnp;% AT: count ALL histones (H3+H4+H2A+H2B+H1)
end;

% peptides
m = 0;
peptides = repmat({''},[npeps,1]);
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
% Fix A-CRIT-004 / C-CRIT-04: orden estable
[~,idx_mf] = sort({matfiles.name});
matfiles = matfiles(idx_mf);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1);
    pepname = matfiles(j).name(1:end-4);
    if 1==strcmp(pepname(1:2),'HH')
        pepname = pepname(2:end);
    end;
    p = strfind(pepname,'_');
    version = pepname(p(1)+1:p(2)-1);
    
    bvar = 0;
    if length(version)>=4
        if 1==strcmp(pepname(1:7),'H3_04v3')
            version = version(4);
        else
            version = version(4:end);
        end;
        histone_pos = [pepname(1:p(1)-1),version,'_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    elseif length(version)==3 && version(3)=='v'
        bvar = 1;
        histone_pos = ['_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    else
        histone_pos = [pepname(1:p(1)-1),'_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    end;
    
    nlen = size(His.pep_mz,1);
    for ino = m+1:m+nlen
        if 1==bvar
            x = strfind(His.mod_short{ino-m},'.');
            if 1==strcmp(pepname(1:3),'H2B')
                peptides{ino,1} = ['H2B',His.mod_short{ino-m}(1:x(1)-1),histone_pos];
            else
                peptides{ino,1} = [His.mod_short{ino-m}(1:x(1)-1),histone_pos];
            end;
        else
            peptides{ino,1} = [histone_pos,' ',His.mod_short{ino-m}];
        end;
    end;
    m = m+nlen;
    if m==npeps
        break;
    end;
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
        if m==npeps
            break;
        end;
    end;
end;

% collect RT, Area, Ratio, respectively
info0 = zeros([npeps,length(raw_names)*3]);
for jno=1:length(raw_names)
    info0(:,jno) = info(:,(jno-1)*3+1);
end;
for jno=1:length(raw_names)
    info0(:,length(raw_names)+jno) = info(:,(jno-1)*3+2);
end;
for jno=1:length(raw_names)
    info0(:,2*length(raw_names)+jno) = info(:,(jno-1)*3+3);
end;
info = info0;
info_ratio = info(1:npeps,2*length(raw_names)+1:3*length(raw_names));

%% output
targets = {'H3K4me1';
    'H3K4me2';
    'H3K4me3';
    'H3K4ac';
    'H3K9me1';
    'H3K9me2';
    'H3K9me3';
    'H3K9ac';
    'H3S10ph';
    'H3K14ac';
    'H3K18me1';
    'H3K18ac';
    'H3K23me1';
    'H3K23ac';
    'H31K27me1';
    'H31K27me2';
    'H31K27me3';
    'H31K27ac';
    'H31K36me1';
    'H31K36me2';
    'H31K36me3';
    'H33K27me1';
    'H33K27me2';
    'H33K27me3';
    'H33K27ac';
    'H33K36me1';
    'H33K36me2';
    'H33K36me3';
    'H3K56me1';
    'H3K56me2';
    'H3K56me3';
    'H3K56ac';
    'H3K79me1';
    'H3K79me2';
    'H3K79me3';
    'H3K79ac';
    'H3K122ac';
    'H4K5ac';
    'H4K8ac';
    'H4K12ac';
    'H4K16ac';
    'H4K20me1';
    'H4K20me2';
    'H4K20me3';
    'H4K20ac'};

% Fix 2026-06-16 (P.G.): añadir targets H2A/H2B/H1 (un target por modulo+marca).
% Recupera hits propios (H2BK4me3, H2A.W K8/K9...) que antes se mezclaban en H3/H4.
extra = {};
for p=1:npeps
    nm = peptides{p};
    sp = strfind(nm,' ');
    if 1==isempty(sp), continue; end;
    modp = nm(1:sp(1)-1);
    if length(modp)<2, continue; end;
    if 1==strcmp(modp(1:2),'H3') || 1==strcmp(modp(1:2),'H4'), continue; end;
    toks = regexp(nm(sp(1)+1:end),'[A-Z]\d+(me1|me2|me3|ac|ph|ub|cr|so)','match');
    for u=1:length(toks)
        key = [modp,' ',toks{u}];
        if 0==any(strcmp(extra,key)) && 0==any(strcmp(targets,key))
            extra{end+1,1} = key; %#ok
        end;
    end;
end;
targets = [targets; extra];

% Fix 2026-06-16 (P.G.): cada target se restringe a su MODULO (region). Antes
% X1=1:npeps para todo salvo H31/H33, de modo que una marca KNN compartida por otra
% histona (p.ej. H2BK4 -> H3K4) contaminaba el target. Targets 'modulo marca' (con
% espacio) son los H2A/H2B/H1 anadidos: region = el propio modulo.
sratios = zeros([length(targets),length(raw_names)]);
for t=1:length(targets)
    c_pep = targets{t};
    sp = strfind(c_pep,' ');
    if 0==isempty(sp)
        region = c_pep(1:sp(1)-1);   c_mod = c_pep(sp(1)+1:end);
    elseif 1==strcmp(c_pep(1:3),'H31')
        region = 'H3_27_40';   c_mod = c_pep(4:end);
    elseif 1==strcmp(c_pep(1:3),'H33')
        region = 'H33_27_40';  c_mod = c_pep(4:end);
    elseif 1==strcmp(c_pep(1:4),'H3K4')
        region = 'H3_3_8';     c_mod = c_pep(3:end);
    elseif 1==strcmp(c_pep(1:4),'H3K9') || 1==strcmp(c_pep(1:5),'H3S10') || 1==strcmp(c_pep(1:5),'H3K14')
        region = 'H3_9_17';    c_mod = c_pep(3:end);
    elseif 1==strcmp(c_pep(1:5),'H3K18') || 1==strcmp(c_pep(1:5),'H3K23')
        region = 'H3_18_26';   c_mod = c_pep(3:end);
    elseif 1==strcmp(c_pep(1:5),'H3K56')
        region = 'H3_53_63';   c_mod = c_pep(3:end);
    elseif 1==strcmp(c_pep(1:5),'H3K79')
        region = 'H3_73_83';   c_mod = c_pep(3:end);
    elseif 1==strcmp(c_pep(1:6),'H3K122')
        region = 'H3_117_128'; c_mod = c_pep(3:end);
    elseif 1==strcmp(c_pep(1:5),'H4K20')
        region = 'H4_20_23';   c_mod = c_pep(3:end);
    elseif 1==strcmp(c_pep(1:3),'H4K')
        region = 'H4_4_17';    c_mod = c_pep(3:end);
    else
        region = '';           c_mod = c_pep(3:end);
    end;

    X1 = [];
    no = 0;
    for p=1:npeps
        if 1==isempty(region) || 0==isempty(strfind(peptides{p},region))
            no = no + 1;
            X1(no) = p;%#ok
        end;
    end;

    X2 = [];
    no = 0;
    for x=1:length(X1)
        if 0==isempty(strfind(peptides{X1(x)},c_mod))
            no = no + 1;
            X2(no) = X1(x);%#ok
        end;
    end;

    if 0==no
        continue;
    end;
    for x=1:length(X2)
        sratios(t,:) = sratios(t,:) + info_ratio(X2(x),:);
    end;
end;

mat_file = fullfile(layout_path,'histone_ratios_single_PTMs.mat');
save(mat_file,'targets','sratios');

inten_file = fullfile(layout_path,'histone_ratios_single_PTMs.xls');
fp = fopen(inten_file,'w');
if -1==fp
    disp(['can not open: ',inten_file]);
    return;
end;

for i=1:length(raw_names)
    fprintf(fp,'\t%d,%s',i,raw_names{i});
end;
fprintf(fp,'\r\n');

for t=1:length(targets)
    fprintf(fp,'%s',targets{t});
    for jno=1:length(raw_names)
        fprintf(fp,'\t%f',sratios(t,jno));
    end;
    fprintf(fp,'\r\n');
end;

fclose(fp);