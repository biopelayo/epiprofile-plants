function OutputFigures(layout_path,raw_names)
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

% total nos
npeps = 0;
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1);
    npeps = npeps + size(His.pep_mz,1);
end;

% separate nos
c_npeps = zeros([2,1]);
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1);
    cnp = size(His.pep_mz,1);
    pepname = matfiles(j).name(1:end-4);
    if 1==strcmp(pepname(1:2),'HH')
        pepname = pepname(2:end);
    end;
    if 1==strcmp(pepname(1:2),'H3')
        c_npeps(1) = c_npeps(1) + cnp;
    elseif 1==strcmp(pepname(1:2),'H4')
        c_npeps(2) = c_npeps(2) + cnp;
    end;
end;

% peptides
m = 0;
peptides = repmat({''},[npeps,1]);
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
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

% info_rt = info(1:npeps,1:length(raw_names));%++
info_auc = info(1:npeps,length(raw_names)+1:2*length(raw_names));
info_ratio = info(1:npeps,2*length(raw_names)+1:3*length(raw_names));

% info_rt = info_rt - floor(max(max(info_rt,[],1))/2);%++
info_ratio = log2((info_ratio+1e-10)*2^5);%++
zscore_ratio = zscore(info_ratio,[],2);


%% figures
% raw_names
for i=1:length(raw_names)
    raw_names{i} = [num2str(i),',',raw_names{i}];
    raw_names{i}(strfind(raw_names{i},'_')) = '-';%++
end;

% fig_path
fig_path = fullfile(layout_path,'heatmap_clustering');
if 0==exist(fig_path,'dir') && 0==mkdir(fig_path)
    fprintf(1,'can not create: %s\n',fig_path);
    return;
end;

%% 01---------------------
% bar
cur_fig = '01_bar_peptide_number';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
set(gcf,'visible','off');

pepnos = zeros([length(raw_names),1]);
for jno=1:length(raw_names)
    pepnos(jno) = length(find(info_auc(1:npeps,jno)>0));
end;
bar(pepnos);
title('bar of peptide number');
xlabel('samples');
ylabel('# histone peptides');
xlim([0 ceil(length(raw_names)*4/3)]);
ylim([0 250]);
text(0,npeps,num2str(npeps));
for jno=1:length(raw_names)
    text(length(raw_names)+1,250-jno*8,raw_names{jno});
end;
% caption
text(0,-24,['There are ',num2str(npeps),' histone peptides theoretically. The number of identified peptides in each sample is shown.'],'FontSize',9);
%text(0,-24,'There are 205 histone peptides theoretically. The number of identified peptides in each sample is shown.','FontSize',9);

print('-dpdf',out_file);
close();

%% 02---------------------
% boxplot
cur_fig = '02_boxplot_peptide_intensity';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
set(gcf,'visible','off');

boxplot(log10(info_auc),'datalim',[2 12]);
title('boxplot of peptide intensity');
xlabel('samples');
ylabel('histone peptide intensity (log_1_0)');
xlim([0 ceil(length(raw_names)*4/3)]);
ylim([0 15]);
set(gca,'YTick',0:2:15);
for jno=1:length(raw_names)
    text(length(raw_names)+1,15-jno*0.5,raw_names{jno});
end;
% caption
text(0,-1.5,'Boxplot produces a distribution of peptide intensity in each sample. There is one box per sample.','FontSize',9);
text(0,-2,'On each box, the central mark is the median, the edges of the box are the 25th and 75th percentiles,','FontSize',9);
text(0,-2.5,'and the whiskers extend to the most extreme data points not considered outliers.','FontSize',9);

print('-dpdf',out_file);
close();

%% 03---------------------
% principal component analysis (pca)
cur_fig = '03_pca';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
set(gcf,'visible','off');

cur_ratio = zscore_ratio;
[wcoeff,score] = pca(cur_ratio');%#ok row is observation (sample)
plot(score(:,1),score(:,2),'+');
title('PCA of samples');
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
xlim([min(score(:,1))-2 max(score(:,1))+2+(max(score(:,1))-min(score(:,1)))/3]);
ylim([min(score(:,2))-2 max(score(:,2))+2]);
if max(score(:,2))+2<50
    step = 0.7;
else
    step = 4;
end;
for jno=1:length(raw_names)
    text(score(jno,1),score(jno,2),num2str(jno),'color','r');
    text(max(score(:,1))+2,max(score(:,2))+2-jno*step,raw_names{jno});
end;
% caption
text(min(score(:,1))-2,min(score(:,2))-2-3*step,'Ratios are projected onto the first two principal components by PCA. Distances or clusters of samples are shown.','FontSize',9);

print('-dpdf',out_file);
close();

%%
ratio_labels = {'0.1%','0.2%','0.4%','0.8%','1.6%','3.125%','6.25%','12.5%','25%','50%','100%'};

%% 11---------------------
% HeatMap
heatmap_histone(fig_path,info_ratio,raw_names,peptides,c_npeps,ratio_labels);

%% 12---------------------
% clustering
% clustering_histone(fig_path,info_ratio,raw_names,peptides,ratio_labels);

write_txt(fig_path);

function heatmap_histone(fig_path,info_ratio,raw_names,peptides,c_npeps,ratio_labels)
%%

% 10---------------------
% HeatMap of single PTMs
mat_file = fullfile(fileparts(fig_path),'histone_ratios_single_PTMs.mat');
if 0~=exist(mat_file,'file')
    load(mat_file);% targets, sratios
    sratios = log2((sratios+1e-10)*2^5);%#ok
    
    % heatmap
    cur_fig = '10_heatmap_ratio_single_PTMs';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    set(gcf,'visible','off');
    
    cur_ratio = sratios;
    hmo = HeatMap(cur_ratio,'RowLabels',targets,'ColumnLabels',raw_names,'DisplayRange',5);
    addTitle(hmo,'HeatMap of ratios for single PTMs','FontSize',10);
    plot(hmo);
    colorbar('Position',[0.91 0.226 0.01 0.698],'YTick',1:11,'YTickLabel',ratio_labels);
    
    print('-dpdf','-r300',out_file);
    close();
    delete(hmo);
    
    % clustergram
    cur_fig = '10_heatmap_zscore_single_PTMs';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    set(gcf,'visible','off');
    
    cur_ratio = zscore(sratios,[],2);
    cgo = clustergram(cur_ratio,'RowLabels',targets,'ColumnLabels',raw_names,'Cluster',1,'DisplayRatio',[1e-6 0.2]);
    addTitle(cgo,'HeatMap of zscores for single PTMs','FontSize',10);
    plot(cgo);
    
    print('-dpdf','-r300',out_file);
    close();
    delete(cgo);
end;

% 11---------------------
all_figs = {'11_heatmap_ratio_H3','11_heatmap_ratio_H4'};
index = [1;cumsum(c_npeps)+1];

for ino=1:length(c_npeps)
    cur_fig = all_figs{ino};
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    set(gcf,'visible','off');
    
    p = strfind(cur_fig,'_');
    cur_fig(p) = '-';
    IX = index(ino):index(ino+1)-1;
    cur_ratio = info_ratio(IX,1:length(raw_names));
    hmo = HeatMap(cur_ratio,'RowLabels',peptides(IX),'ColumnLabels',raw_names,'DisplayRange',5);
    addTitle(hmo,['HeatMap of ratios for ',cur_fig(p(3)+1:end)],'FontSize',10);
    plot(hmo);
    colorbar('Position',[0.91 0.226 0.01 0.698],'YTick',1:11,'YTickLabel',ratio_labels);
    
    print('-dpdf','-r300',out_file);
    close();
    delete(hmo);
end;

function clustering_histone(fig_path,info_ratio,raw_names,peptides,ratio_labels)%#ok
%%

zscore_ratio = zscore(info_ratio,[],2);

% 12---------------------
cur_fig = '12_clustering_peptides';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
set(gcf,'visible','off');

nK = 6;
[cidx,ctrs] = kmeans(zscore_ratio,nK,'rep',5);
for c = 1:nK
    subplot(2,nK,c);
    plot(1:length(raw_names),zscore_ratio(cidx==c,:)');
    axis tight;
end;
for c = 1:nK
    subplot(2,nK,nK+c);
    plot(1:length(raw_names),ctrs(c,:)');
    axis tight;
    axis off;
end;
suptitle('Clustering of peptides');

print('-dpdf','-r300',out_file);
close();

% 13---------------------
for ino=1:nK
    cur_ratio = info_ratio(cidx==ino,:);
    
    cur_fig = ['13_cluster_',num2str(ino)];
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    set(gcf,'visible','off');
    
    hmo = HeatMap(cur_ratio,'RowLabels',peptides(cidx==ino),'ColumnLabels',raw_names,'DisplayRange',5);
    addTitle(hmo,['cluster ',num2str(ino)],'FontSize',10);
    plot(hmo);
    colorbar('Position',[0.91 0.226 0.01 0.698],'YTick',1:11,'YTickLabel',ratio_labels);
    
    print('-dpdf','-r300',out_file);
    close();
    delete(hmo);
end;

function write_txt(fig_path)
%%

txt_file = fullfile(fig_path,'Figure Legends.txt');

fp = fopen(txt_file,'w');
if -1==fp
    fprintf(1,'can not open: %s\r\n',txt_file);
    return;
end;

fprintf(fp,'Figure 01: the bar plot represents the number of histone peptides quantified using EpiProfile. \r\n');
fprintf(fp,'The number on the top left of the graph represents all detectable peptides, while each bar \r\n');
fprintf(fp,'represents how many of those peptides were detected with a quantification different than zero \r\n');
fprintf(fp,'in each sample. Sample legend is listed on the top right of the plot. This graph should be used \r\n');
fprintf(fp,'as quality control; too many missing values are indicative of low sensitive analysis.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 02: the box plot displays the Log10 transformed total peptide intensity. The central red \r\n');
fprintf(fp,'line indicates the median abundance across all peptides, and the vertical width the dynamic range \r\n');
fprintf(fp,'of the quantification. Sample legend is listed on the top right of the plot. The graph should be \r\n');
fprintf(fp,'used as quality control of the experiments, i.e. all box plots should have a comparable central \r\n');
fprintf(fp,'value and width. If this is not the case the low abundance samples were probably injected in a \r\n');
fprintf(fp,'lower amount than recommended.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 03: principal component analysis (PCA) of all performed runs. The graph displays in two \r\n');
fprintf(fp,'dimensions of the n-dimensional dataset, aiming to simplify the distance between samples into a \r\n');
fprintf(fp,'spatial 2D graph. Replicates of the same sample should cluster close to each other, while different \r\n');
fprintf(fp,'condition should be further apart. If replicates do not cluster it might be appropriate to verify \r\n');
fprintf(fp,'the quality of the LC-MS runs, as statistics will be poor in case of non reproducible analyses.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 10: the heatmap represents the relative abundance of single histone PTMs, quantified using \r\n');
fprintf(fp,'the area detected for the respective peptides. For PTMs that belong into multiple combinations \r\n');
fprintf(fp,'(e.g. K9me3 and K9me3K14ac) the relative abundance was obtained by summing the relative abundance \r\n');
fprintf(fp,'of all peptides containing the given modification. Color coding represents the percentage of \r\n');
fprintf(fp,'occupancy for each given mark (only histone H3 and H4 marks are displayed).\r\n');

fclose(fp);