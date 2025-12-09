function OutputFigures(layout_path,raw_names)
%% OutputFigures — Create QC/summary figures from EpiProfile outputs (MAT-files).
%  This function expects that, for each raw file (sample), you already have a
%  folder "<layout_path>/<index>_<raw_name>/detail" containing multiple *.mat
%  files produced by prior steps (e.g., peptide windows). Each *.mat must
%  define at least:
%     - His: a struct with fields like 'pep_mz' (N_peptides x 1)
%     - auc: a (N_peptides x 3) matrix with [RT, Area, Ratio] columns
%     - mod_short: (optional) a cell array with short modification labels (one per peptide)
%
%  The function:
%    1) Counts the total number of peptides ("npeps") from the FIRST sample.
%    2) Computes per-histone-subclass peptide counts (c_npeps) for reporting.
%    3) Builds a peptide label list ("peptides") based on file names and mod_short.
%    4) Aggregates [RT, Area, Ratio] across all samples into a big matrix.
%    5) Extracts and transforms Area/Ratio matrices for plotting (log/scale).
%    6) Generates PDF figures: bar of identified peptides, boxplots, PCA, and
%       heatmaps (if auxiliary MAT exists). Also writes a "Figure Legends.txt".
%
%  IMPORTANT design detail:
%    - The total peptide list and counts are derived from the FIRST sample.
%      This assumes that the peptide definition/order is consistent across samples.

%% nos, peptides, and info
i = 1; % We start by inspecting ONLY the first raw file to discover peptide layout.
cur_rawname = raw_names{i};
if i<10
    prefix = ['0',num2str(i)]; % Left-pad index with a zero for 1–9 (e.g., "01", "02")
else
    prefix = num2str(i);
end;
cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail'); % "<layout_path>/<01_raw>/detail"

% total nos
npeps = 0; % Will hold the TOTAL number of peptides across all MATs for the FIRST sample
out_file1 = fullfile(cur_outpath,'*.mat'); % Enumerate *.mat peptide-window files
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1); % Must load 'His' with field 'pep_mz'
    % Each window contributes size(His.pep_mz,1) peptides (rows).
    npeps = npeps + size(His.pep_mz,1);
end;

% separate nos
% c_npeps partitions peptide counts into 9 buckets by histone class/subclass:
% [H3 group 1 (no<=3), H3 group 2 (no==4), H3 group 3 (else),
%  H4, H1, H2A group 1 (no<=3), H2A group 2 (no==4), H2A group 3 (else), H2B]
c_npeps = zeros([9,1]);
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1); % 'His' must exist again in each file
    cnp = size(His.pep_mz,1); % peptide count for this window
    pepname = matfiles(j).name(1:end-4); % filename without ".mat"
    if 1==strcmp(pepname(1:2),'HH')
        pepname = pepname(2:end); % normalize names that start with "HH" (drop leading H)
    end;
    p = strfind(pepname,'_');        % underscores split tokens within the filename
    version = pepname(p(1)+1:p(2)-1);% token between first and second underscore (e.g., "01", "02a", "04v3"...)
    curno = str2double(version(1:2));% numeric window id = first two chars (e.g., "01"->1)
    if 1==strcmp(pepname(1:2),'H3')
        % H3 is split into 3 groups based on curno
        if curno<=3
            c_npeps(1) = c_npeps(1) + cnp;
        elseif curno==4
            c_npeps(2) = c_npeps(2) + cnp;
        else
            c_npeps(3) = c_npeps(3) + cnp;
        end;
    elseif 1==strcmp(pepname(1:2),'H4')
        c_npeps(4) = c_npeps(4) + cnp;
    elseif 1==strcmp(pepname(1:2),'H1')
        c_npeps(5) = c_npeps(5) + cnp;
    elseif 1==strcmp(pepname(1:3),'H2A')
        % H2A is also split into 3 groups based on curno
        if curno<=3
            c_npeps(6) = c_npeps(6) + cnp;
        elseif curno==4
            c_npeps(7) = c_npeps(7) + cnp;
        else
            c_npeps(8) = c_npeps(8) + cnp;
        end;
    elseif 1==strcmp(pepname(1:3),'H2B')
        c_npeps(9) = c_npeps(9) + cnp;
    end;
end;

% peptides
% Build a cell array of human-readable peptide identifiers (N = npeps),
% based on filename tokens + the per-peptide short modification label.
m = 0; % running index into the big peptide list
peptides = repmat({''},[npeps,1]);
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1); % must define 'His' and often 'mod_short'
    pepname = matfiles(j).name(1:end-4);
    if 1==strcmp(pepname(1:2),'HH')
        pepname = pepname(2:end); % normalize 'HH' prefix to 'H'
    end;
    p = strfind(pepname,'_');             % underscores positions
    version = pepname(p(1)+1:p(2)-1);     % middle token with window/version info
    
    bvar = 0; % flag for a "variant" naming pattern (length(version)==3 && version(3)=='v')
    if length(version)>=4
        % Special case: files like "H3_04v3_*" compress version to a single char
        if 1==strcmp(pepname(1:7),'H3_04v3')
            version = version(4);
        else
            version = version(4:end); % otherwise drop "??v" and keep the suffix
        end;
        % histone_pos forms a canonical "H?_X_Y_Z" position descriptor
        histone_pos = [pepname(1:p(1)-1),version,'_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    elseif length(version)==3 && version(3)=='v'
        % If version is exactly 3 chars and ends with 'v', mark as variant
        bvar = 1;
        % Build a shorter pos tag: just "_Y_Z"
        histone_pos = ['_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    else
        % Standard case: "H?_X_Y_Z"
        histone_pos = [pepname(1:p(1)-1),'_',pepname(p(2)+1:p(3)-1),'_',pepname(p(3)+1:end)];
    end;
    
    nlen = size(His.pep_mz,1); % number of peptides in this MAT (window)
    for ino = m+1:m+nlen
        if 1==bvar
            % Variant style: label = (optionally "H2B") + mod_short (without its .suffix) + histone_pos
            x = strfind(His.mod_short{ino-m},'.'); % find '.' to strip tail from mod_short
            if 1==strcmp(pepname(1:3),'H2B')
                peptides{ino,1} = ['H2B',His.mod_short{ino-m}(1:x(1)-1),histone_pos];
            else
                peptides{ino,1} = [His.mod_short{ino-m}(1:x(1)-1),histone_pos];
            end;
        else
            % Default style: label = histone_pos + space + full mod_short
            peptides{ino,1} = [histone_pos,' ',His.mod_short{ino-m}];
        end;
    end;
    m = m+nlen; % advance global peptide index
end;

% info
% Build the big data matrix "info" of size (npeps x (#raw_names*3)),
% concatenating the 3 columns [RT, Area, Ratio] from each sample:
% sample1(RT,Area,Ratio) | sample2(RT,Area,Ratio) | ...
info = zeros([npeps,length(raw_names)*3]);
for i=1:length(raw_names)
    cur_rawname = raw_names{i};
    if i<10
        prefix = ['0',num2str(i)];
    else
        prefix = num2str(i);
    end;
    cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');
    m = 0; % running row pointer within "info"
    for j = 1:length(matfiles)
        matfile1 = fullfile(cur_outpath,matfiles(j).name);
        load(matfile1); % must define 'auc' with 3 columns: [RT, Area, Ratio]
        nlen = size(His.pep_mz,1);
        % Place the 3 columns of 'auc' in the corresponding sample block:
        info(m+1:m+nlen,(i-1)*3+1:(i-1)*3+3) = auc;
        m = m+nlen;
    end;
end;

% collect RT, Area, Ratio, respectively
% Reorder columns so that we get blocks of [RT-block | Area-block | Ratio-block]
info0 = zeros([npeps,length(raw_names)*3]);
for jno=1:length(raw_names)
    info0(:,jno) = info(:,(jno-1)*3+1); % RT columns packed first
end;
for jno=1:length(raw_names)
    info0(:,length(raw_names)+jno) = info(:,(jno-1)*3+2); % Area columns next
end;
for jno=1:length(raw_names)
    info0(:,2*length(raw_names)+jno) = info(:,(jno-1)*3+3); % Ratio columns last
end;
info = info0; % overwrite with reordered version

% info_rt = info(1:npeps,1:length(raw_names));%++  (kept commented in the original)
info_auc = info(1:npeps,length(raw_names)+1:2*length(raw_names));       % Area-only block
info_ratio = info(1:npeps,2*length(raw_names)+1:3*length(raw_names));   % Ratio-only block

% info_rt = info_rt - floor(max(max(info_rt,[],1))/2);%++ (kept commented in the original)
% Stabilize and scale ratios: shift by a tiny +1e-10 to avoid log2(0), multiply by 2^5,
% then take log2. This rescales the ratio dynamic range for visualization.
info_ratio = log2((info_ratio+1e-10)*2^5);%++
% Z-score per peptide (row-wise): center and scale each peptide across samples.
zscore_ratio = zscore(info_ratio,[],2);

%% figures
% raw_names
% Prettify sample names: prefix with "index," and replace underscores with hyphens
for i=1:length(raw_names)
    raw_names{i} = [num2str(i),',',raw_names{i}];
    raw_names{i}(strfind(raw_names{i},'_')) = '-';%++
end;

% fig_path
% Ensure output folder exists: "<layout_path>/heatmap_clustering"
fig_path = fullfile(layout_path,'heatmap_clustering');
if 0==exist(fig_path,'dir') && 0==mkdir(fig_path)
    fprintf(1,'can not create: %s\n',fig_path);
    return;
end;

%% 01---------------------
% bar — Number of identified peptides per sample
cur_fig = '01_bar_peptide_number';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
set(gcf,'visible','off');

% Count peptides with Area>0 per sample as "identified"
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
text(0,npeps,num2str(npeps)); % annotate total theoretical peptides
for jno=1:length(raw_names)
    text(length(raw_names)+1,250-jno*8,raw_names{jno}); % draw legend on the right margin
end;
% caption
text(0,-24,['There are ',num2str(npeps),' histone peptides theoretically. The number of identified peptides in each sample is shown.'],'FontSize',9);
%text(0,-24,'There are 205 histone peptides theoretically. The number of identified peptides in each sample is shown.','FontSize',9);

print('-dpdf',out_file);
close();

%% 02---------------------
% boxplot — Distribution of log10(Area) across peptides per sample
cur_fig = '02_boxplot_peptide_intensity';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
set(gcf,'visible','off');

boxplot(log10(info_auc),'datalim',[2 12]); % log10 area; limit whiskers for consistent display
title('boxplot of peptide intensity');
xlabel('samples');
ylabel('histone peptide intensity (log_1_0)');
xlim([0 ceil(length(raw_names)*4/3)]);
ylim([0 15]);
set(gca,'YTick',0:2:15);
for jno=1:length(raw_names)
    text(length(raw_names)+1,15-jno*0.5,raw_names{jno}); % legend on the right
end;
% caption
text(0,-1.5,'Boxplot produces a distribution of peptide intensity in each sample. There is one box per sample.','FontSize',9);
text(0,-2,'On each box, the central mark is the median, the edges of the box are the 25th and 75th percentiles,','FontSize',9);
text(0,-2.5,'and the whiskers extend to the most extreme data points not considered outliers.','FontSize',9);

print('-dpdf',out_file);
close();

%% 03---------------------
% PCA — Principal Component Analysis on the z-scored ratio matrix
cur_fig = '03_pca';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
set(gcf,'visible','off');

cur_ratio = zscore_ratio; % rows=peptides, cols=samples (we will transpose)
[wcoeff,score] = pca(cur_ratio');%#ok row is observation (sample)
% 'score' are sample coordinates in PC space; plot PC1 vs PC2
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
    text(score(jno,1),score(jno,2),num2str(jno),'color','r'); % numeric sample id near each point
    text(max(score(:,1))+2,max(score(:,2))+2-jno*step,raw_names{jno}); % legend block at top-right
end;
% caption
text(min(score(:,1))-2,min(score(:,2))-2-3*step,'Ratios are projected onto the first two principal components by PCA. Distances or clusters of samples are shown.','FontSize',9);

print('-dpdf',out_file);
close();

%%
% Labels for percentage/ratio colorbars used in heatmaps
ratio_labels = {'0.1%','0.2%','0.4%','0.8%','1.6%','3.125%','6.25%','12.5%','25%','50%','100%'};

%% 11---------------------
% HeatMap (delegated to helper function) — uses "histone_ratios_single_PTMs.mat" if present
heatmap_histone(fig_path,info_ratio,raw_names,peptides,c_npeps,ratio_labels);

%% 12---------------------
% Unsupervised clustering (disabled by default; left as in the original)
% clustering_histone(fig_path,info_ratio,raw_names,peptides,ratio_labels);

% Finally, write a text file with figure legends for Figures 01–03 and 10
write_txt(fig_path);

function heatmap_histone(fig_path,info_ratio,raw_names,peptides,c_npeps,ratio_labels)%#ok
%%
% Builds heatmaps if the auxiliary file "<layout>/histone_ratios_single_PTMs.mat"
% exists. That file must define:
%   targets: row labels (e.g., specific PTMs)
%   sratios: matrix (PTMs x samples) of ratios in linear scale
% The function log2-transforms sratios (scaled) and draws:
%   - a "ratio" heatmap with absolute levels
%   - a "z-score" heatmap (row-wise) for pattern emphasis

% 10---------------------
% HeatMap of single PTMs
mat_file = fullfile(fileparts(fig_path),'histone_ratios_single_PTMs.mat');
if 0~=exist(mat_file,'file')
    load(mat_file);% targets, sratios
    sratios = log2((sratios+1e-10)*2^5);%#ok  % rescale and log2 like info_ratio
    
    % heatmap — absolute (log2-scaled) ratios
    cur_fig = '10_heatmap_ratio_single_PTMs';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    set(gcf,'visible','off');
    
    cur_ratio = sratios;
    % HeatMap/clustergram are Bioinformatics Toolbox components
    hmo = HeatMap(cur_ratio,'RowLabels',targets,'ColumnLabels',raw_names,'DisplayRange',5);
    addTitle(hmo,'HeatMap of ratios for single PTMs','FontSize',10);
    plot(hmo);
    colorbar('Position',[0.91 0.226 0.01 0.698],'YTick',1:11,'YTickLabel',ratio_labels);
    
    print('-dpdf','-r300',out_file);
    close();
    delete(hmo);
    
    % clustergram — z-score per PTM for pattern comparison
    cur_fig = '10_heatmap_zscore_single_PTMs';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    set(gcf,'visible','off');
    
    cur_ratio = zscore(sratios,[],2); % z-score per PTM (row)
    cgo = clustergram(cur_ratio,'RowLabels',targets,'ColumnLabels',raw_names,'Cluster',1,'DisplayRatio',[1e-6 0.2]);
    addTitle(cgo,'HeatMap of zscores for single PTMs','FontSize',10);
    plot(cgo);
    
    print('-dpdf','-r300',out_file);
    close();
    delete(cgo);
end;

% 11---------------------
%{
% (Original block left commented as-is; would draw per-subclass heatmaps)
all_figs = {'11_heatmap_ratio_H3_1','11_heatmap_ratio_H3_2','11_heatmap_ratio_H3_3','11_heatmap_ratio_H4','11_heatmap_ratio_HH1','11_heatmap_ratio_HH2A_1','11_heatmap_ratio_HH2A_2','11_heatmap_ratio_HH2A_3','11_heatmap_ratio_HH2B'};
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
    
    if 4==ino && 0==sum(c_npeps(ino+1:end))
        return;
    end;
end;
%}

function clustering_histone(fig_path,info_ratio,raw_names,peptides,ratio_labels)%#ok
%%
% K-means clustering of peptides by their z-scored ratio profiles across samples.
% Produces:
%   (12) a multi-panel figure with raw members (top) and cluster centroids (bottom)
%   (13) one heatmap per cluster with row labels = peptide IDs

zscore_ratio = zscore(info_ratio,[],2);

% 12---------------------
cur_fig = '12_clustering_peptides';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
set(gcf,'visible','off');

nK = 6; % number of clusters
[cidx,ctrs] = kmeans(zscore_ratio,nK,'rep',5); % 5 replicates for more stable solution
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
% Writes a human-readable legend file describing Figures 01–03 and 10.
% This is useful for reports and for non-experts to interpret the outputs.

txt_file = fullfile(fig_path,'Figure Legends.txt');

fp = fopen(txt_file,'w');
if -1==fp
    fprintf(1,'can not open: %s\r\n',txt_file);
    return;
end;

function epiplants_superplots_qc(fig_path, info_auc, info_ratio, zscore_ratio, raw_names, peptides)
% Superplots QC — crea PDFs adicionales, no sustituye a tus 01–03:
%  S01_sorted_barh_identifications.pdf  — barra horizontal ordenada (#IDs + %missing)
%  S02_boxswarm_intensity.pdf           — boxchart + swarm submuestreado (log10 Area)
%  S03_pca_enhanced.pdf                 — PCA con % var explicada y etiquetas
%  S04_corr_heatmap_spearman.pdf        — heatmap Spearman reordenado por clustering
%  S05_missingness_bar.pdf              — % missing por muestra
%  S06_pca_scree.pdf                    — scree plot de PCA

    if 0==exist(fig_path,'dir'), mkdir(fig_path); end

    % -------- helpers locales --------
    function save_pdf(fig, name)
        out = fullfile(fig_path, name);
        try exportgraphics(fig,out,'ContentType','vector','BackgroundColor','white'); catch, print(fig,'-dpdf',out); end
        close(fig);
    end
    function cols = nice_colors(n)
        base = [33 158 188; 251 133 0; 39 125 161; 233 196 106; 244 162 97; ...
                231 111 81; 87 117 144; 47 62 70; 69 123 157; 29 53 87]/255;
        if n<=size(base,1), cols=base(1:n,:); else, cols=repmat(base,ceil(n/size(base,1)),1); cols=cols(1:n,:); end
    end

    % Preparación
    [npeps, nsamp] = size(info_auc);
    logI = log10(info_auc+1);
    ids  = sum(info_auc>0,1)';         % # identificados por muestra
    miss = mean(info_auc<=0,1)'*100;   % % missing

    % ---- S01: barra horizontal ordenada ----
    [ids_sorted, ord] = sort(ids,'descend');
    names_sorted = raw_names(ord); miss_sorted = miss(ord);
    f1 = figure('Visible','off','Color','w','Position',[100 100 900 700]);
    ax = axes(f1); barh(ax, ids_sorted,'FaceColor',[0.13 0.55 0.80]); grid(ax,'on'); box(ax,'off');
    set(ax,'YTick',1:nsamp,'YTickLabel',names_sorted,'FontName','Arial');
    xlabel(ax,'# péptidos identificados (Area>0)'); ylabel(ax,'muestras (ordenadas)');
    xlim(ax,[0 max(ids_sorted)*1.08]);
    for k=1:nsamp
        text(ax, ids_sorted(k)+max(ids_sorted)*0.01, k, sprintf('%.1f%% missing',miss_sorted(k)), ...
            'HorizontalAlignment','left','FontSize',9);
    end
    title(ax,sprintf('Identificaciones por muestra (N teórico = %d)',npeps));
    save_pdf(f1,'S01_sorted_barh_identifications.pdf');

    % ---- S02: box + swarm submuestreado ----
    max_pts = 2000;
    f2 = figure('Visible','off','Color','w','Position',[100 100 1100 700]);
    ax = axes(f2); hold(ax,'on');
    X = repelem(1:nsamp, npeps)'; Y = logI(:);
    try
        boxchart(ax, X, Y, 'BoxWidth',0.5,'MarkerStyle','none');
    catch
        boxplot(ax, Y, X, 'Symbol','.');
    end
    rng(1);
    for s=1:nsamp
        yi = logI(:,s); xi = s*ones(size(yi));
        if numel(yi)>max_pts, idx = randperm(numel(yi),max_pts); yi=yi(idx); xi=xi(idx); end
        try swarmchart(ax, xi, yi, 6, 'filled','MarkerFaceAlpha',0.25);
        catch scatter(ax, xi, yi, 6, 'filled','MarkerFaceAlpha',0.25); end
    end
    set(ax,'XLim',[0.3 nsamp+0.7],'XTick',1:nsamp,'XTickLabel',raw_names,'XTickLabelRotation',90);
    ylabel(ax,'log_{10}(Area+1)'); title(ax,'Intensidad por muestra (box + swarm)');
    grid(ax,'on'); box(ax,'off');
    save_pdf(f2,'S02_boxswarm_intensity.pdf');

    % ---- S03: PCA mejorada ----
    [coeff,score,latent,~,explained] = pca(zscore_ratio'); %#ok<ASGLU>
    f3 = figure('Visible','off','Color','w','Position',[100 100 950 750]);
    ax = axes(f3); hold(ax,'on');
    scatter(ax, score(:,1), score(:,2), 40, nice_colors(nsamp), 'filled');
    for j=1:nsamp, text(ax, score(j,1), score(j,2), sprintf(' %d', j),'Color',[0.2 0.2 0.2],'FontSize',9); end
    xlabel(ax, sprintf('PC1 (%.1f%%)',explained(1))); ylabel(ax, sprintf('PC2 (%.1f%%)',explained(2)));
    title(ax,'PCA de muestras (ratios z-score por péptido)'); grid(ax,'on'); box(ax,'off'); axis(ax,'equal'); axis(ax,'tight');
    save_pdf(f3,'S03_pca_enhanced.pdf');

    % ---- S04: heatmap de correlación (Spearman) con reordenación ----
    try
        R = corr(logI,'Type','Spearman','Rows','pairwise');
        try
            D = 1 - R; D(1:nsamp+1:end) = 0;
            Z = linkage(squareform((D+D')/2),'average'); ord2 = optimalleaforder(Z, squareform((D+D')/2));
        catch
            ord2 = 1:nsamp;
        end
        f4 = figure('Visible','off','Color','w','Position',[100 100 950 850]);
        h = heatmap(f4, raw_names(ord2), raw_names(ord2), R(ord2,ord2), 'CellLabelColor','none');
        try h.Colormap = turbo; catch, h.Colormap = parula; end
        h.ColorLimits = [0 1]; h.XLabel = 'muestras'; h.YLabel = 'muestras'; h.Title = 'Spearman (log_{10}(Area+1))';
        exportgraphics(f4, fullfile(fig_path,'S04_corr_heatmap_spearman.pdf'), 'ContentType','vector','BackgroundColor','white');
        close(f4);
    catch ME
        fprintf(1,'[superplots] corr/heatmap: %s\n', ME.message);
    end

    % ---- S05: %missing por muestra ----
    f5 = figure('Visible','off','Color','w','Position',[100 100 900 600]);
    ax = axes(f5); bar(ax, miss, 'FaceColor',[0.85 0.10 0.12]); grid(ax,'on'); box(ax,'off');
    set(ax,'XTick',1:nsamp,'XTickLabel',raw_names,'XTickLabelRotation',90);
    ylabel(ax,'% missing (Area==0)'); title(ax,'Porcentaje de valores ausentes por muestra');
    ylim(ax,[0 max(5, ceil(max(miss)/5)*5)]);
    save_pdf(f5,'S05_missingness_bar.pdf');

    % ---- S06: scree plot PCA ----
    f6 = figure('Visible','off','Color','w','Position',[100 100 900 600]);
    ax = axes(f6); plot(ax, cumsum(explained), '-o','LineWidth',1.5,'MarkerSize',5); grid(ax,'on'); box(ax,'off');
    xlabel(ax,'# componentes'); ylabel(ax,'Varianza explicada acumulada (%)');
    title(ax,'Scree plot (PCA sobre ratios z-score)'); xlim(ax,[1 max(3, min(nsamp,10))]);
    save_pdf(f6,'S06_pca_scree.pdf');
end


function epiplants_superplots_embeddings(fig_path, zscore_ratio, raw_names)
% Superplots embeddings y selección de k (añadidos):
%  S07_k_silhouette_curve.pdf   — curva silhouette media para elegir k (péptidos)
%  S07b_silhouette_bestk.pdf    — diagrama silhouette en k*
%  S08_pca_conf_ellipses.pdf    — PCA con elipses 95% por grupo inferido
%  S09_tsne_enhanced.pdf        — t-SNE de muestras (si existe tsne)
%  S10_umap_enhanced.pdf        — UMAP de muestras (si existe run_umap/umap)

    if 0==exist(fig_path,'dir'), mkdir(fig_path); end

    % -------- helpers --------
    function save_pdf(fig, name)
        out = fullfile(fig_path, name);
        try exportgraphics(fig,out,'ContentType','vector','BackgroundColor','white'); catch, print(fig,'-dpdf',out); end
        close(fig);
    end
    function cols = nice_colors(n)
        base = [33 158 188; 251 133 0; 39 125 161; 233 196 106; 244 162 97; ...
                231 111 81; 87 117 144; 47 62 70; 69 123 157; 29 53 87]/255;
        if n<=size(base,1), cols=base(1:n,:); else, cols=repmat(base,ceil(n/size(base,1)),1); cols=cols(1:n,:); end
    end
    function [grp, grp_names] = infer_groups_from_names(names)
        n = numel(names); grp = strings(n,1);
        for i=1:n
            nm = string(names{i});
            p = strfind(nm, ','); if ~isempty(p), nm = extractAfter(nm, p(1)); end
            parts = split(nm,'-'); grp(i) = strtrim(parts(1));
            if strlength(grp(i))==0, grp(i) = "grp1"; end
        end
        grp_names = unique(grp,'stable');
    end
    function draw_cov_ellipse(ax, mu, Sigma, chi2val, edgeColor)
        if any(~isfinite(Sigma), 'all') || rank(Sigma)<2, Sigma = Sigma + 1e-6*eye(2); end
        [V,L] = eig(Sigma,'vector'); U = V * diag(sqrt(max(L,0)));
        t = linspace(0,2*pi,200); E = sqrt(chi2val) * U * [cos(t); sin(t)];
        plot(ax, mu(1)+E(1,:), mu(2)+E(2,:), 'Color',edgeColor, 'LineWidth',1.5);
    end

    % ===== 07) elección de k con silhouette (clustering de péptidos) =====
    X = zscore_ratio;   % (n_pep x n_samp)
    Ks = 2:min(10, max(2, size(X,1)-1));
    mean_s = nan(size(Ks)); bestk=NaN; bestS=-Inf; bestC=[]; bestDist='correlation';
    for ii=1:numel(Ks)
        k = Ks(ii);
        try
            cidx = kmeans(X, k, 'Replicates',5,'MaxIter',1000,'OnlinePhase','off','Distance','correlation');
            s = silhouette(X, cidx, 'correlation');
        catch
            cidx = kmeans(X, k, 'Replicates',5,'MaxIter',1000,'OnlinePhase','off','Distance','sqeuclidean');
            s = silhouette(X, cidx, 'sqeuclidean'); bestDist='sqeuclidean';
        end
        mean_s(ii) = mean(s, 'omitnan');
        if mean_s(ii)>bestS, bestS=mean_s(ii); bestk=k; bestC=cidx; end
    end
    f7 = figure('Visible','off','Color','w','Position',[80 80 900 600]);
    plot(Ks, mean_s,'-o','LineWidth',1.6); grid on; box off;
    xlabel('k (clusters de péptidos)'); ylabel('Silhouette media');
    title(sprintf('Selección de k — %s (k*=%d, S=%.3f)', bestDist, bestk, bestS));
    save_pdf(f7,'S07_k_silhouette_curve.pdf');

    try
        f7b = figure('Visible','off','Color','w','Position',[80 80 900 700]);
        silhouette(X, bestC, bestDist);
        title(sprintf('Silhouette en k*=%d (%s)', bestk, bestDist));
        save_pdf(f7b,'S07b_silhouette_bestk.pdf');
    catch
        % si silhouette no está disponible, se omite
    end

    % ===== 08) PCA con elipses 95% por grupo inferido de nombres =====
    [grp, grp_names] = infer_groups_from_names(raw_names);
    [~,score,~,~,expl] = pca(X');   % muestras como observaciones
    f8 = figure('Visible','off','Color','w','Position',[80 80 950 750]);
    ax = axes(f8); hold(ax,'on'); C = nice_colors(numel(grp_names));
    for g=1:numel(grp_names)
        idx = grp==grp_names(g);
        scatter(ax, score(idx,1), score(idx,2), 50, C(g,:), 'filled', 'MarkerFaceAlpha',0.85);
        if nnz(idx)>=2
            mu = mean(score(idx,1:2),1); Sigma = cov(score(idx,1:2));
            draw_cov_ellipse(ax, mu(:), Sigma, 5.991, C(g,:)); % chi2 df=2 (95%)
        end
    end
    for j=1:size(score,1), text(ax, score(j,1), score(j,2), sprintf(' %d', j),'Color',[0.1 0.1 0.1],'FontSize',8); end
    xlabel(ax, sprintf('PC1 (%.1f%%)', expl(1))); ylabel(ax, sprintf('PC2 (%.1f%%)', expl(2)));
    legend(ax, cellstr(grp_names), 'Location','bestoutside'); grid(ax,'on'); box(ax,'off'); axis(ax,'equal'); axis(ax,'tight');
    title(ax,'PCA con elipses de confianza 95% (grupos inferidos)');
    save_pdf(f8,'S08_pca_conf_ellipses.pdf');

    % ===== 09) t-SNE (si existe) =====
    if exist('tsne','file')==2
        try
            rng(1);
            perplex = max(5, min(30, floor((size(X,2)-1)/3)));
            Y = tsne(X','Distance','correlation','Standardize',true,'Perplexity',perplex,'NumDimensions',2);
            f9 = figure('Visible','off','Color','w','Position',[80 80 900 700]);
            ax = axes(f9); hold(ax,'on'); C = nice_colors(numel(grp_names));
            for g=1:numel(grp_names)
                idx = grp==grp_names(g);
                scatter(ax, Y(idx,1), Y(idx,2), 50, C(g,:), 'filled', 'MarkerFaceAlpha',0.85);
            end
            for j=1:size(Y,1), text(ax, Y(j,1), Y(j,2), sprintf(' %d', j),'Color',[0.15 0.15 0.15],'FontSize',8); end
            xlabel(ax,'t-SNE 1'); ylabel(ax,'t-SNE 2'); grid(ax,'on'); box(ax,'off');
            legend(ax, cellstr(grp_names), 'Location','bestoutside');
            title(ax, sprintf('t-SNE (perplexity=%d, distancia=correlation)', perplex));
            save_pdf(f9,'S09_tsne_enhanced.pdf');
        catch ME
            fprintf(1,'[superplots] t-SNE: %s\n', ME.message);
        end
    end

    % ===== 10) UMAP (si existe run_umap/umap) =====
    if exist('run_umap','file')==2 || exist('umap','file')==2
        try
            rng(1); nn = max(10, min(15, size(X,2)-1));
            if exist('run_umap','file')==2
                out = run_umap(X','metric','correlation','n_neighbors',nn,'min_dist',0.3,'verbose','none');
                if isstruct(out) && isfield(out,'embedding'), U = out.embedding; else, U = []; end
            else
                U = umap(X','n_neighbors',nn); % alternativa
            end
            if ~isempty(U)
                f10 = figure('Visible','off','Color','w','Position',[80 80 900 700]);
                ax = axes(f10); hold(ax,'on'); C = nice_colors(numel(grp_names));
                for g=1:numel(grp_names)
                    idx = grp==grp_names(g);
                    scatter(ax, U(idx,1), U(idx,2), 50, C(g,:), 'filled', 'MarkerFaceAlpha',0.85);
                end
                for j=1:size(U,1), text(ax, U(j,1), U(j,2), sprintf(' %d', j),'Color',[0.15 0.15 0.15],'FontSize',8); end
                xlabel(ax,'UMAP 1'); ylabel(ax,'UMAP 2'); grid(ax,'on'); box(ax,'off');
                legend(ax, cellstr(grp_names), 'Location','bestoutside');
                title(ax, sprintf('UMAP (n\\_neighbors=%d, metric=correlation)', nn));
                save_pdf(f10,'S10_umap_enhanced.pdf');
            end
        catch ME
            fprintf(1,'[superplots] UMAP: %s\n', ME.message);
        end
    end
end

% === SUPERPLOTS (añadidos, no sustituyen nada) ============================
try
    epiplants_superplots_qc(fig_path, info_auc, info_ratio, zscore_ratio, raw_names, peptides);
    epiplants_superplots_embeddings(fig_path, zscore_ratio, raw_names);
catch ME
    fprintf(1,'[superplots] %s\n', ME.message);
end

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
end