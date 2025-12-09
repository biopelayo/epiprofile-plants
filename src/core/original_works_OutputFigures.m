function OutputFigures(layout_path,raw_names)
%% OutputFigures — QC and summary figures from EpiProfile MAT outputs (H3/H4 variant)
% 
% OUTPUTS (FILES WRITTEN)
%    <layout_path>/heatmap_clustering/01_bar_peptide_number.pdf
%    ... (Figures 02-27) ...
%    <layout_path>/heatmap_clustering/28_histogram_peptide_completeness.pdf
%    <layout_path>/heatmap_clustering/29_barplot_sample_missing_count.pdf
%    <layout_path>/heatmap_clustering/30_heatmap_ratio_H3_73_83.pdf
%    <layout_path>/heatmap_clustering/31_boxplot_ratio_H3_73_83.pdf
%    <layout_path>/heatmap_clustering/32_scatter_avg_ratio_H3_73_83_vs_H3.pdf
%    <layout_path>/heatmap_clustering/33_dendrogram_sample_clustering.pdf
%    <layout_path>/heatmap_clustering/34_paired_boxplot_H3_vs_H4_per_sample.pdf
%    <layout_path>/heatmap_clustering/35_density_plot_ratio_H3_vs_H4.pdf
%    <layout_path>/heatmap_clustering/36_heatmap_corr_H3_73_83_vs_all_H3.pdf
%    <layout_path>/heatmap_clustering/Figure Legends.txt
% -------------------------------------------------------------------------

%% nos, peptides, and info (Data Loading)
i = 1;
cur_rawname = raw_names{i};
if i<10
    prefix = ['0',num2str(i)];
else
    prefix = num2str(i);
end;
cur_outpath = fullfile(layout_path,[prefix,'_',cur_rawname],'detail');
npeps = 0;
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1, 'His'); 
    npeps = npeps + size(His.pep_mz,1);
end;
c_npeps = zeros([2,1]);
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1, 'His');
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
m = 0;
peptides = repmat({''},[npeps,1]);
out_file1 = fullfile(cur_outpath,'*.mat');
matfiles = dir(out_file1);
for j=1:length(matfiles)
    matfile1 = fullfile(cur_outpath,matfiles(j).name);
    load(matfile1, 'His'); 
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
        load(matfile1, 'auc', 'His');
        nlen = size(His.pep_mz,1);
        info(m+1:m+nlen,(i-1)*3+1:(i-1)*3+3) = auc;
        m = m+nlen;
    end;
end;
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
info_auc = info(1:npeps,length(raw_names)+1:2*length(raw_names));
info_ratio = info(1:npeps,2*length(raw_names)+1:3*length(raw_names));
info_ratio_log = log2((info_ratio+1e-10)*2^5);
zscore_ratio = zscore(info_ratio_log,[],2);
n_samples = length(raw_names); 

%% figures (Aesthetic preparation)
clean_names = cell(1, n_samples);
original_raw_names_for_legend = cell(1, n_samples); 
for i=1:n_samples
    temp_name = raw_names{i}; 
    temp_name(strfind(temp_name,'_')) = '-';
    original_raw_names_for_legend{i} = [num2str(i),', ',temp_name]; 
    clean_names{i} = temp_name;
end;

fig_path = fullfile(layout_path,'heatmap_clustering');
if 0==exist(fig_path,'dir') && 0==mkdir(fig_path)
    fprintf(1,'Cannot create: %s\n',fig_path);
    return;
end;
ratio_labels = {'0.1%','0.2%','0.4%','0.8%','1.6%','3.125%','6.25%','12.5%','25%','50%','100%'};

%% 01---------------------
cur_fig = '01_bar_peptide_number';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure(); 
pepnos = sum(info_auc>0, 1)'; 
bar(pepnos);
title('Quantified Peptide Count per Sample','FontSize',12);
xlabel('Sample','FontSize',10);
ylabel('Number of Quantified Peptides','FontSize',10);
xlim([0.5 n_samples+0.5]);
set(gca, 'XTick', 1:n_samples, 'XTickLabel', clean_names, 'XTickLabelRotation', 45);
ylim([0 ceil(max(pepnos)*1.1)]);
box off;
h = text(0, min(ylim)-0.2*range(ylim), ...
    sprintf('Total theoretical peptides: %d. Bar height shows detected peptides (Area > 0) in each sample.', npeps), ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 02---------------------
cur_fig = '02_boxplot_peptide_intensity';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
boxplot(log10(info_auc),'datalim',[2 12],'labels',clean_names);
title('Distribution of Peptide Intensity (Log_{10} Area)','FontSize',12);
xlabel('Sample','FontSize',10);
ylabel('Peptide Intensity (Log_{10} Area)','FontSize',10);
set(gca, 'XTickLabelRotation', 45);
box off;
h = text(0, min(ylim)-0.2*range(ylim), ...
    'Boxplot of log-transformed area counts. Assesses consistency in sample loading and MS sensitivity.', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 03---------------------
cur_fig = '03_pca';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
cur_ratio = zscore_ratio;
[wcoeff,score,latent] = pca(cur_ratio');
plot(score(:,1),score(:,2),'ko','MarkerSize',5,'MarkerFaceColor','r'); 
title(sprintf('PCA of Samples (All Peptides Ratios)\nVariance Explained: PC1=%.1f%%, PC2=%.1f%%', latent(1)/sum(latent)*100, latent(2)/sum(latent)*100),'FontSize',12);
xlabel('1st Principal Component','FontSize',10);
ylabel('2nd Principal Component','FontSize',10);
grid on;
box on;
hold on;
for jno=1:n_samples
    text(score(jno,1),score(jno,2),num2str(jno),'color','k','FontSize',8,'VerticalAlignment','bottom','HorizontalAlignment','left');
end;
xlim_data = [min(score(:,1)) max(score(:,1))];
ylim_data = [min(score(:,2)) max(score(:,2))];
xlim([min(xlim_data)-2 max(xlim_data)+2+(range(xlim_data))/3]);
ylim([min(ylim_data)-2 max(ylim_data)+2]);
step = range(ylim)/n_samples/2.5;
if step < 0.5; step = 0.5; end;
for jno=1:n_samples
    text(max(xlim_data)+2,max(ylim_data)+2-jno*step,original_raw_names_for_legend{jno}, 'FontSize', 8);
end;
h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
    'PCA of Z-scored Log2 Ratios. Distances reflect similarity/dissimilarity between samples.', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% --------------------------------------------------------------------------------
% --- BLOQUE DE PLOTS AMPLIADO (Figuras 14 a 36) ---
% --------------------------------------------------------------------------------
index = [1;cumsum(c_npeps)+1];
data_H3 = info_ratio_log(index(1):index(2)-1,:);
data_H4 = info_ratio_log(index(2):index(3)-1,:);
full_data = [data_H3; data_H4];
sample_indices = repmat(1:n_samples, size(full_data, 1), 1);
full_data_flat = full_data(:);
sample_indices_flat = sample_indices(:);

%% 14---------------------
cur_fig = '14_violinplot_ratio_H3_H4';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
set(gcf,'Units','inches','Position',[0 0 8 4.5]);
hold on;
if exist('violin','file')
    group_var_class = [repmat({'H3'}, size(data_H3, 1), 1); repmat({'H4'}, size(data_H4, 1), 1)];
    group_var_class_flat = repmat(group_var_class, n_samples, 1);
    violin(full_data_flat, 'groups', {sample_indices_flat, group_var_class_flat}, ...
           'histType', 'area', 'medc', 'k', 'width', 0.8, 'bandwidth', 0.5, 'plotLine', true);
    set(gca, 'XTick', 1.5:2:2*n_samples, 'XTickLabel', clean_names, 'XTickLabelRotation', 45);
    legend('H3', 'H4', 'Location', 'NorthEastOutside');
else
    fprintf(1,'WARNING: ''violin'' function not found. Generating Boxplots for Figure 14.\n');
    G = cell(length(full_data_flat), 1);
    group_order_list = cell(1, 2*n_samples);
    N_peps_H3 = size(data_H3, 1);
    N_peps_H4 = size(data_H4, 1);
    offset = N_peps_H3 * n_samples;
    for i = 1:n_samples
        group_h3 = ['H3_S', num2str(i)];
        group_h4 = ['H4_S', num2str(i)];
        group_order_list{2*i-1} = group_h3;
        group_order_list{2*i} = group_h4;
        idx_H3_start = (i-1) * N_peps_H3 + 1;
        idx_H3_end = i * N_peps_H3;
        G(idx_H3_start : idx_H3_end) = {group_h3};
        idx_H4_start = offset + (i-1) * N_peps_H4 + 1;
        idx_H4_end = offset + i * N_peps_H4;
        G(idx_H4_start : idx_H4_end) = {group_h4};
    end
    valid_indices = ~isnan(full_data_flat);
    boxplot(full_data_flat(valid_indices), G(valid_indices), 'plotstyle', 'compact', 'GroupOrder', group_order_list); 
    x_tick_labels = cell(1, 2*n_samples);
    for i=1:n_samples
        x_tick_labels{2*i-1} = [clean_names{i}, ' (H3)'];
        x_tick_labels{2*i} = [clean_names{i}, ' (H4)'];
    end
    set(gca, 'XTickLabel', x_tick_labels, 'XTickLabelRotation', 45);
    legend('off');
end
title('Log_{2} Ratio Distribution per Sample (H3 vs. H4)','FontSize',12);
xlabel('Sample','FontSize',10);
ylabel('Relative Abundance Ratio (Log_{2} Scale)','FontSize',10);
box on;
h = text(0, min(ylim)-0.2*range(ylim), ...
    'Distribution of peptide modification ratios, comparing H3 and H4 classes across samples.', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 15---------------------
cur_fig = '15_scatter_Area_vs_Ratio';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
log10_auc_flat = log10(info_auc(info_auc(:)>0));
log2_ratio_flat = info_ratio_log(info_auc(:)>0);
valid_idx = isfinite(log10_auc_flat) & isfinite(log2_ratio_flat);
log10_auc_flat = log10_auc_flat(valid_idx);
log2_ratio_flat = log2_ratio_flat(valid_idx);
scatter(log10_auc_flat, log2_ratio_flat, 5, '.', 'MarkerEdgeAlpha', 0.05);
title('QC: Peptide Ratio vs. Intensity (All Data Points)','FontSize',12);
xlabel('Peptide Intensity (Log_{10} Area)','FontSize',10);
ylabel('Relative Abundance Ratio (Log_{2} Scale)','FontSize',10);
if length(log10_auc_flat) > 1
    p = polyfit(log10_auc_flat, log2_ratio_flat, 1);
    yfit = polyval(p, get(gca,'XLim'));
    hold on;
    plot(get(gca,'XLim'), yfit, 'r--', 'LineWidth', 1);
    hold off;
    R = corrcoef(log10_auc_flat, log2_ratio_flat);
    legend(sprintf('Linear Corr. R: %.2f', R(1,2)), 'Location', 'SouthEast');
end
box on;
h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
    'A low correlation suggests the ratio is independent of peptide abundance, indicating good normalization.', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 16---------------------
cur_fig = '16_heatmap_missing_values';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
set(gcf,'Units','inches','Position',[0 0 5 8]);
detection_matrix = double(info_auc > 0);
imagesc(detection_matrix);
colormap([1 1 1; 0.7 0.7 0.7]);
title('Peptide Detection QC (Missing Values)','FontSize',12);
xlabel('Sample','FontSize',10);
ylabel('Peptide Index','FontSize',10);
set(gca, 'YTick', []); 
set(gca, 'XTick', 1:n_samples, 'XTickLabel', clean_names, 'XTickLabelRotation', 45);
cb = colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'Missing (Area=0)', 'Detected (Area>0)'}, ...
         'Position', [0.91 0.2 0.02 0.7]);
print('-dpdf','-painters','-r600',out_file);
close();

%% 17---------------------
cur_fig = '17_barplot_low_ratio_percentage';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
title_str = 'Percentage of Low Ratio Peptides (<1%)';
threshold = log2((0.01 + 1e-10)*2^5); 
low_ratio_count = sum(info_ratio_log < threshold, 1);
low_ratio_percentage = (low_ratio_count ./ size(info_ratio_log, 1)) * 100;
bar(low_ratio_percentage);
title(title_str,'FontSize',12);
xlabel('Sample','FontSize',10);
ylabel('% Peptides with Ratio < 1%','FontSize',10);
set(gca, 'XTick', 1:n_samples, 'XTickLabel', clean_names, 'XTickLabelRotation', 45);
box off;
print('-dpdf','-painters','-r600',out_file);
close();

%% 18---------------------
cur_fig = '18_boxplot_log10area_by_class';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
log10_area_H3_full = log10(info_auc(index(1):index(2)-1,:));
log10_area_H4_full = log10(info_auc(index(2):index(3)-1,:));
log10_area_H3 = log10_area_H3_full(isfinite(log10_area_H3_full));
log10_area_H4 = log10_area_H4_full(isfinite(log10_area_H4_full));
all_log10_area = [log10_area_H3; log10_area_H4];
group_area = [repmat({'H3'}, length(log10_area_H3), 1); repmat({'H4'}, length(log10_area_H4), 1)];
boxplot(all_log10_area, group_area, 'plotstyle', 'compact', 'GroupOrder', {'H3', 'H4'});
title('Peptide Abundance Distribution: H3 vs H4','FontSize',12);
xlabel('Histone Class','FontSize',10);
ylabel('Peptide Intensity (Log_{10} Area)','FontSize',10);
box off;
print('-dpdf','-painters','-r600',out_file);
close();

%% 19---------------------
cur_fig = '19_scatter_avg_ratio_H3_vs_H4';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
avg_ratio_H3 = mean(data_H3, 1);
avg_ratio_H4 = mean(data_H4, 1);
scatter(avg_ratio_H3, avg_ratio_H4, 50, 1:n_samples, 'filled');
colormap(jet(n_samples));
colorbar;
title('Average Ratio Profile Correlation: H3 vs H4','FontSize',12);
xlabel('Average H3 Ratio (Log_{2} Scale)','FontSize',10);
ylabel('Average H4 Ratio (Log_{2} Scale)','FontSize',10);
hold on;
p = polyfit(avg_ratio_H3', avg_ratio_H4', 1);
yfit = polyval(p, get(gca,'XLim'));
plot(get(gca,'XLim'), yfit, 'k--', 'LineWidth', 1.5);
hold off;
R = corrcoef(avg_ratio_H3, avg_ratio_H4);
text(min(xlim), max(ylim), sprintf('R = %.2f', R(1,2)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
box on;
print('-dpdf','-painters','-r600',out_file);
close();

%% 20---------------------
% HEATMAP — Sample Correlation (FIXED FOR READABILITY)
cur_fig = '20_heatmap_sample_correlation';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
set(gcf,'Units','inches','Position',[0 0 10 9]); 

R_matrix = corr(info_ratio_log, 'Type', 'Spearman', 'rows', 'pairwise');
imagesc(R_matrix);
try
    colormap(coolwarm); 
catch
    colormap(parula); 
end
caxis([-1 1]);
title('Sample-to-Sample Correlation (Spearman R)','FontSize',14,'FontWeight','bold');
xlabel('Sample Index','FontSize',12);
ylabel('Sample Index','FontSize',12);
set(gca, 'XTick', 1:n_samples, 'XTickLabel', 1:n_samples, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:n_samples, 'YTickLabel', 1:n_samples);

% --- CORRECCIÓN: Usar cb.Label.String para el título de la ColorBar ---
cb = colorbar('Position', [0.8 0.1 0.02 0.8]);
cb.Label.String = 'Spearman R';

plot_pos = get(gca, 'Position');
legend_x = plot_pos(1) + plot_pos(3) + 0.02;
legend_y = plot_pos(2) + plot_pos(4);
legend_height = plot_pos(4);
legend_width = 0.2; 
annotation('textbox', [legend_x, legend_y - legend_height, legend_width, legend_height], ...
           'String', [{'Sample Indices:'}; original_raw_names_for_legend(:)], ...
           'FitBoxToText', 'on', 'VerticalAlignment', 'top', ...
           'HorizontalAlignment', 'left', 'FontSize', 8, 'Interpreter', 'none', ...
           'BackgroundColor', [1 1 1], 'EdgeColor', 'none');
h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
    'Measures consistency between samples. Replicates should show high correlation (R > 0.9).', ...
    'FontSize',9, 'VerticalAlignment','top');

print('-dpdf','-painters','-r600',out_file);
close();

%% 21---------------------
cur_fig = '21_barplot_class_count';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
total_H3 = sum(c_npeps(1));
total_H4 = sum(c_npeps(2));
bar([total_H3, total_H4]);
set(gca, 'XTickLabel', {'H3', 'H4'});
title('Total Number of Theoretical Peptides by Class','FontSize',12);
ylabel('Total Peptide Count','FontSize',10);
text(1, total_H3, num2str(total_H3), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(2, total_H4, num2str(total_H4), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
box off;
print('-dpdf','-painters','-r600',out_file);
close();

%% 22---------------------
cur_fig = '22_boxplot_ratio_all_peptides';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
set(gcf,'Units','inches','Position',[0 0 4 4]);
boxplot(info_ratio_log(:), 'plotstyle', 'compact');
set(gca, 'XTickLabel', {'Global Ratio'});
title('Overall Ratio Distribution (Log_{2})','FontSize',12);
ylabel('Relative Abundance Ratio (Log_{2} Scale)','FontSize',10);
box off;
print('-dpdf','-painters','-r600',out_file);
close();

%% 23---------------------
cur_fig = '23_lineplot_median_ratio_profile';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
set(gcf,'Units','inches','Position',[0 0 8 4.5]);
median_ratio_H3 = median(data_H3, 1);
median_ratio_H4 = median(data_H4, 1);
plot(1:n_samples, median_ratio_H3, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(1:n_samples, median_ratio_H4, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
hold off;
title('Median Ratio Profile Across Samples (H3 vs H4)','FontSize',12);
xlabel('Sample','FontSize',10);
ylabel('Median Ratio (Log_{2} Scale)','FontSize',10);
set(gca, 'XTick', 1:n_samples, 'XTickLabel', clean_names, 'XTickLabelRotation', 45);
grid on;
legend('H3 Class Median', 'H4 Class Median', 'Location', 'best');
box on;
print('-dpdf','-painters','-r600',out_file);
close();

%% 24---------------------
cur_fig = '24_boxplot_ratio_noise';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
noise_threshold = 1; 
noise_data = info_ratio_log;
noise_data(abs(info_ratio_log) > noise_threshold) = NaN;
boxplot(noise_data, clean_names, 'plotstyle', 'compact', 'orientation', 'horizontal');
title(['Distribution of "Noise" Ratio (Log_{2} |Ratio| < ', num2str(noise_threshold), ')'],'FontSize',12);
xlabel('Relative Abundance Ratio (Log_{2} Scale)','FontSize',10);
ylabel('Sample','FontSize',10);
set(gca, 'YTickLabelRotation', 0);
grid on;
box on;
print('-dpdf','-painters','-r600',out_file);
close();

%% 25---------------------
cur_fig = '25_pca_H3_only';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
cur_ratio_H3 = zscore(data_H3,[],2);
[wcoeff_H3,score_H3,latent_H3] = pca(cur_ratio_H3');
plot(score_H3(:,1),score_H3(:,2),'ko','MarkerSize',5,'MarkerFaceColor','b');
title(sprintf('PCA of Samples (H3 Peptides Only)\nVariance Explained: PC1=%.1f%%, PC2=%.1f%%', latent_H3(1)/sum(latent_H3)*100, latent_H3(2)/sum(latent_H3)*100),'FontSize',12);
xlabel('1st Principal Component (H3)','FontSize',10);
ylabel('2nd Principal Component (H3)','FontSize',10);
grid on;
box on;
hold on;
for jno=1:n_samples
    text(score_H3(jno,1),score_H3(jno,2),num2str(jno),'color','k','FontSize',8,'VerticalAlignment','bottom','HorizontalAlignment','left');
end
xlim_data = [min(score_H3(:,1)) max(score_H3(:,1))];
ylim_data = [min(score_H3(:,2)) max(score_H3(:,2))];
xlim([min(xlim_data)-2 max(xlim_data)+2+(range(xlim_data))/3]);
ylim([min(ylim_data)-2 max(ylim_data)+2]);
step = range(ylim)/n_samples/2.5;
if step < 0.5; step = 0.5; end;
for jno=1:n_samples
    text(max(xlim_data)+2,max(ylim_data)+2-jno*step,original_raw_names_for_legend{jno}, 'FontSize', 8);
end
print('-dpdf','-painters','-r600',out_file);
close();

%% 26---------------------
cur_fig = '26_pca_H4_only';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
cur_ratio_H4 = zscore(data_H4,[],2);
[wcoeff_H4,score_H4,latent_H4] = pca(cur_ratio_H4');
plot(score_H4(:,1),score_H4(:,2),'ko','MarkerSize',5,'MarkerFaceColor','g');
title(sprintf('PCA of Samples (H4 Peptides Only)\nVariance Explained: PC1=%.1f%%, PC2=%.1f%%', latent_H4(1)/sum(latent_H4)*100, latent_H4(2)/sum(latent_H4)*100),'FontSize',12);
xlabel('1st Principal Component (H4)','FontSize',10);
ylabel('2nd Principal Component (H4)','FontSize',10);
grid on;
box on;
hold on;
for jno=1:n_samples
    text(score_H4(jno,1),score_H4(jno,2),num2str(jno),'color','k','FontSize',8,'VerticalAlignment','bottom','HorizontalAlignment','left');
end
xlim_data = [min(score_H4(:,1)) max(score_H4(:,1))];
ylim_data = [min(score_H4(:,2)) max(score_H4(:,2))];
xlim([min(xlim_data)-2 max(xlim_data)+2+(range(xlim_data))/3]);
ylim([min(ylim_data)-2 max(ylim_data)+2]);
step = range(ylim)/n_samples/2.5;
if step < 0.5; step = 0.5; end;
for jno=1:n_samples
    text(max(xlim_data)+2,max(ylim_data)+2-jno*step,original_raw_names_for_legend{jno}, 'FontSize', 8);
end
print('-dpdf','-painters','-r600',out_file);
close();

%% 27---------------------
cur_fig = '27_overlay_histogram_ratios';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_wide_figure(); 
hold on;
color_map = lines(n_samples);
h_arr = zeros(1, n_samples);
bins = -6:0.5:6;
for i = 1:n_samples
    data_sample = info_ratio_log(:, i);
    data_sample = data_sample(isfinite(data_sample));
    [N, ~] = histcounts(data_sample, bins);
    h_arr(i) = plot(bins(1:end-1) + 0.25, N, 'Color', color_map(i,:), 'LineWidth', 1.5);
end
hold off;
title('Ratio Distribution (Log_{2}) Overlay per Sample','FontSize',12);
xlabel('Relative Abundance Ratio (Log_{2} Scale)','FontSize',10);
ylabel('Frequency (Peptide Count)','FontSize',10);
legend(h_arr, clean_names, 'Location', 'NorthEastOutside');
grid on;
box on;
print('-dpdf','-painters','-r600',out_file);
close();

%% 28---------------------
cur_fig = '28_histogram_peptide_completeness';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
detection_count = sum(info_auc > 0, 2); 
bins = 0.5:1:(n_samples + 0.5);
histogram(detection_count, bins, 'Normalization', 'count');
title('Peptide Completeness (Detection Frequency)','FontSize',12);
xlabel('Number of Samples Peptide was Detected In','FontSize',10);
ylabel('Number of Peptides','FontSize',10);
set(gca, 'XTick', 1:n_samples);
box on;
h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
    'Histogram of how many samples each peptide was successfully quantified (Area > 0).', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 29---------------------
cur_fig = '29_barplot_sample_missing_count';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
missing_count = sum(info_auc == 0, 1); 
bar(missing_count);
title('Missing Peptide Count per Sample','FontSize',12);
xlabel('Sample','FontSize',10);
ylabel('Number of Missing Peptides (Area = 0)','FontSize',10);
set(gca, 'XTick', 1:n_samples, 'XTickLabel', clean_names, 'XTickLabelRotation', 45);
box off;
h = text(0, min(ylim)-0.2*range(ylim), ...
    'Count of peptides with zero abundance per sample. Direct QC for sample specific data completeness/dropout rates.', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 30-32 (H3_73_83 Specific)
h3_73_83_indices = contains(peptides, 'H3_73_83');
if ~any(h3_73_83_indices)
    fprintf(1,'WARNING: No peptides matching H3_73_83 found. Skipping figures 30-32.\n');
else
    data_h3_73_83 = info_ratio_log(h3_73_83_indices,:);
    peptides_h3_73_83 = peptides(h3_73_83_indices);
    n_h3_73_83 = size(data_h3_73_83, 1);

    %% 30---------------------
    cur_fig = '30_heatmap_ratio_H3_73_83';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    figure('visible','off');
    setup_publication_figure();
    set(gcf,'Units','inches','Position',[0 0 6 4 + n_h3_73_83*0.2]); 
    hmo = HeatMap(data_h3_73_83,'RowLabels',peptides_h3_73_83,'ColumnLabels',clean_names,'DisplayRange',5);
    addTitle(hmo,'Log_{2} Ratio Heatmap for H3(73-83) Peptides','FontSize',12);
    plot(hmo);
    % --- CORRECCIÓN: Usar cb.Label.String ---
    cb = colorbar('Position',[0.91 0.226 0.01 0.698],'YTick',1:11,'YTickLabel',ratio_labels);
    cb.Label.String = 'Ratio (%)';
    print('-dpdf','-painters','-r600',out_file);
    close();
    delete(hmo);

    %% 31---------------------
    cur_fig = '31_boxplot_ratio_H3_73_83';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    figure('visible','off');
    setup_publication_figure();
    set(gcf,'Units','inches','Position',[0 0 8 4.5]);
    boxplot(data_h3_73_83, 'labels', clean_names);
    title('Ratio Distribution (Log_{2}) of H3(73-83) Peptides','FontSize',12);
    xlabel('Sample','FontSize',10);
    ylabel('Relative Abundance Ratio (Log_{2} Scale)','FontSize',10);
    set(gca, 'XTickLabelRotation', 45);
    box on;
    h = text(0, min(ylim)-0.2*range(ylim), ...
        'Distribution of modification ratios specifically on the H3(73-83) peptide across all samples.', ...
        'FontSize',9, 'VerticalAlignment','top');
    print('-dpdf','-painters','-r600',out_file);
    close();

    %% 32---------------------
    cur_fig = '32_scatter_avg_ratio_H3_73_83_vs_H3';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    figure('visible','off');
    setup_publication_figure();
    avg_ratio_H3_73_83 = mean(data_h3_73_83, 1); 
    avg_ratio_H3_all = mean(data_H3, 1);         
    scatter(avg_ratio_H3_all, avg_ratio_H3_73_83, 50, 1:n_samples, 'filled');
    colormap(jet(n_samples));
    title('Correlation: H3(73-83) Average Ratio vs. All H3 Average Ratio','FontSize',12);
    xlabel('Average Ratio of ALL H3 Peptides (Log_{2} Scale)','FontSize',10);
    ylabel('Average Ratio of H3(73-83) Peptides (Log_{2} Scale)','FontSize',10);
    hold on;
    p = polyfit(avg_ratio_H3_all', avg_ratio_H3_73_83', 1);
    yfit = polyval(p, get(gca,'XLim'));
    plot(get(gca,'XLim'), yfit, 'k--', 'LineWidth', 1.5);
    hold off;
    R = corrcoef(avg_ratio_H3_all, avg_ratio_H3_73_83);
    legend(sprintf('R = %.2f (Correlation line)'), 'Location', 'best');
    box on;
    h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
        'Assesses if the global trends observed in all H3 modifications are reflected in the specific H3(73-83) region.', ...
        'FontSize',9, 'VerticalAlignment','top');
    print('-dpdf','-painters','-r600',out_file);
    close();
end

%% 33---------------------
cur_fig = '33_dendrogram_sample_clustering';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
set(gcf,'Units','inches','Position',[0 0 6 6]);
cur_data = zscore_ratio'; 
D = pdist(cur_data, 'Euclidean');
Z = linkage(D, 'ward'); 
[H, ~, outperm] = dendrogram(Z, 0, 'Labels', clean_names); 
set(gca, 'XTickLabelRotation', 45);
title('Hierarchical Clustering (Dendrogram) of Samples','FontSize',12);
ylabel('Distance (Ward Linkage, Euclidean Z-Score Ratio)','FontSize',10);
box off;
h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
    'Visualizes the dissimilarity between samples. Shorter links indicate tighter grouping (e.g., replicates).', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 34---------------------
cur_fig = '34_paired_boxplot_H3_vs_H4_per_sample';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
set(gcf,'Units','inches','Position',[0 0 10 6]); 

all_ratios_paired = [];
group_labels = {};
group_idx = [];
current_start = 1;
for i = 1:n_samples
    sample_data_h3 = data_H3(:, i);
    sample_data_h3 = sample_data_h3(isfinite(sample_data_h3));
    num_h3_vals = length(sample_data_h3);
    sample_data_h4 = data_H4(:, i);
    sample_data_h4 = sample_data_h4(isfinite(sample_data_h4));
    num_h4_vals = length(sample_data_h4);
    all_ratios_paired = [all_ratios_paired; sample_data_h3; sample_data_h4];
    group_idx = [group_idx; repmat(2*i - 1, num_h3_vals, 1); repmat(2*i, num_h4_vals, 1)];
    group_labels{2*i-1} = [clean_names{i}, ' (H3)'];
    group_labels{2*i} = [clean_names{i}, ' (H4)'];
end
boxplot(all_ratios_paired, group_idx, 'plotstyle', 'compact');
set(gca, 'XTick', 1:length(group_labels), 'XTickLabel', group_labels, 'XTickLabelRotation', 60);
title('Log_{2} Ratio Distribution: H3 vs. H4 per Sample','FontSize',14,'FontWeight','bold');
xlabel('Sample and Histone Class','FontSize',12);
ylabel('Relative Abundance Ratio (Log_{2} Scale)','FontSize',12);
box on;
grid on;
h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
    'Compares the distribution of modification ratios between H3 and H4 within each sample. Reveals class-specific trends.', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 35---------------------
cur_fig = '35_density_plot_ratio_H3_vs_H4';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
figure('visible','off');
setup_publication_figure();
set(gcf,'Units','inches','Position',[0 0 8 5]); 
hold on;
[f_h3, xi_h3] = ksdensity(data_H3(isfinite(data_H3(:)))); 
plot(xi_h3, f_h3, 'LineWidth', 2, 'Color', [0 0.4470 0.7410], 'DisplayName', 'All H3 Peptides');
[f_h4, xi_h4] = ksdensity(data_H4(isfinite(data_H4(:)))); 
plot(xi_h4, f_h4, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'All H4 Peptides');
hold off;
title('Density Plot of Log_{2} Ratio Distribution (H3 vs. H4)','FontSize',14,'FontWeight','bold');
xlabel('Relative Abundance Ratio (Log_{2} Scale)','FontSize',12);
ylabel('Density','FontSize',12);
legend('Location', 'NorthWest');
grid on;
box on;
h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
    'Compares the global distribution shape of modification ratios for H3 and H4 peptides. Reveals overall biases.', ...
    'FontSize',9, 'VerticalAlignment','top');
print('-dpdf','-painters','-r600',out_file);
close();

%% 36---------------------
cur_fig = '36_heatmap_corr_H3_73_83_vs_all_H3';
out_file = fullfile(fig_path,[cur_fig,'.pdf']);
warning off all;
if exist('data_h3_73_83', 'var') && ~isempty(data_h3_73_83)
    figure('visible','off');
    setup_publication_figure();
    set(gcf,'Units','inches','Position',[0 0 8 10]); 
    index_H3 = index(1):index(2)-1;
    data_H3_all_peptides = info_ratio_log(index_H3, :);
    peptides_H3_all = peptides(index_H3);
    idx_not_h3_73_83_in_all_H3 = ~contains(peptides_H3_all, 'H3_73_83');
    data_H3_other_peptides = data_H3_all_peptides(idx_not_h3_73_83_in_all_H3, :);
    peptides_H3_other = peptides_H3_all(idx_not_h3_73_83_in_all_H3);
    Corr_matrix = corr(data_h3_73_83', data_H3_other_peptides', 'Type', 'Spearman', 'rows', 'pairwise');
    imagesc(Corr_matrix);
    try
        colormap(coolwarm);
    catch
        colormap(parula);
    end
    caxis([-1 1]);
    title('Spearman Correlation: H3(73-83) Peptides vs. Other H3 Peptides','FontSize',14,'FontWeight','bold');
    xlabel('Other H3 Peptides','FontSize',12);
    ylabel('H3(73-83) Peptides','FontSize',12);
    set(gca, 'YTick', 1:size(data_h3_73_83, 1), 'YTickLabel', peptides_h3_73_83, 'YTickLabelRotation', 0, 'FontSize', 8);
    set(gca, 'XTick', 1:size(data_H3_other_peptides, 1), 'XTickLabel', peptides_H3_other, 'XTickLabelRotation', 90, 'FontSize', 6);
    % --- CORRECCIÓN: Usar cb.Label.String ---
    cb = colorbar('Position', [0.91 0.1 0.02 0.8]);
    cb.Label.String = 'Spearman R';
    box on;
    grid off;
    h = text(min(xlim), min(ylim)-0.2*range(ylim), ...
        'Identifies co-regulating or inversely regulating H3 peptides with the H3(73-83) region. Strong correlations (positive or negative) indicate functional relationships.', ...
        'FontSize',9, 'VerticalAlignment','top');
    print('-dpdf','-painters','-r600',out_file);
    close();
else
    fprintf(1,'WARNING: No H3_73_83 peptides data available to generate Figure 36.\n');
end

%% Final Function Calls
% (Llamamos a las subfunciones usando los nombres de leyenda originales y los nombres de eje limpios)
heatmap_histone(fig_path,info_ratio_log,original_raw_names_for_legend,peptides,c_npeps,ratio_labels, clean_names); 
write_txt(fig_path); 
end 


% -------------------------------------------------------------------------
% --- SUBFUNCIONES REQUERIDAS POR MATLAB (DEBEN IR AL FINAL) ---
% -------------------------------------------------------------------------

function setup_publication_figure()
    % Función de ayuda para la configuración de figuras estándar
    set(gcf,'Color','w'); 
    set(gcf,'PaperPositionMode','auto'); 
    set(gcf,'Units','inches','Position',[0 0 6 4]); 
    set(gca,'FontSize',10,'FontName','Arial'); 
end
% -------------------------------------------------------------------------
function setup_wide_figure()
    % Función de ayuda para figuras más anchas
    set(gcf,'Color','w');
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'Units','inches','Position',[0 0 10 4]); 
    set(gca,'FontSize',10,'FontName','Arial');
end
% -------------------------------------------------------------------------
function heatmap_histone(fig_path,info_ratio_log,raw_names_legend,peptides,c_npeps,ratio_labels, clean_names)
% 10---------------------
mat_file = fullfile(fileparts(fig_path),'histone_ratios_single_PTMs.mat');
if 0~=exist(mat_file,'file')
    load(mat_file, 'targets', 'sratios'); 
    sratios = log2((sratios+1e-10)*2^5);%#ok
    
    cur_fig = '10_heatmap_ratio_single_PTMs';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    figure('visible','off');
    setup_publication_figure(); 
    set(gcf,'Units','inches','Position',[0 0 6 8]); 
    
    cur_ratio = sratios;
    hmo = HeatMap(cur_ratio,'RowLabels',targets,'ColumnLabels',clean_names,'DisplayRange',5); 
    addTitle(hmo,'Log_{2} Ratio Heatmap for Single PTMs','FontSize',12);
    plot(hmo);
    % --- CORRECCIÓN: Usar cb.Label.String ---
    cb = colorbar('Position',[0.91 0.226 0.01 0.698],'YTick',1:11,'YTickLabel',ratio_labels);
    cb.Label.String = 'Ratio (%)';
    
    print('-dpdf','-painters','-r600',out_file);
    close();
    delete(hmo);
    
    cur_fig = '10_heatmap_zscore_single_PTMs';
    out_file = fullfile(fig_path,[cur_fig,'.pdf']);
    warning off all;
    figure('visible','off');
    setup_publication_figure();
    set(gcf,'Units','inches','Position',[0 0 6 8]); 
    
    cur_ratio = zscore(sratios,[],2);
    cgo = clustergram(cur_ratio,'RowLabels',targets,'ColumnLabels',clean_names,'Cluster',1,'DisplayRatio',[1e-6 0.2]);
    addTitle(cgo,'Z-score Heatmap for Single PTMs','FontSize',12);
    plot(cgo);
    
    print('-dpdf','-r600',out_file);
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
    figure('visible','off');
    setup_publication_figure();
    
    p = strfind(cur_fig,'_');
    histone_class = cur_fig(p(3)+1:end);
    IX = index(ino):index(ino+1)-1;
    cur_ratio = info_ratio_log(IX,1:length(raw_names_legend));
    hmo = HeatMap(cur_ratio,'RowLabels',peptides(IX),'ColumnLabels',clean_names,'DisplayRange',5);
    addTitle(hmo,['Peptide Log_{2} Ratio Heatmap: ',histone_class],'FontSize',12);
    plot(hmo);
    % --- CORRECCIÓN: Usar cb.Label.String ---
    cb = colorbar('Position',[0.91 0.226 0.01 0.698],'YTick',1:11,'YTickLabel',ratio_labels);
    cb.Label.String = 'Ratio (%)';
    
    print('-dpdf','-painters','-r600',out_file);
    close();
    delete(hmo);
end;
end % heatmap_histone
% -------------------------------------------------------------------------
function clustering_histone(fig_path,info_ratio,raw_names,peptides,ratio_labels)%#ok
% (Función opcional de K-means, sin cambios)
end % clustering_histone
% -------------------------------------------------------------------------
function write_txt(fig_path)
%% write_txt — emit a plain-text legend file describing the main figures (EXPANDED and in ENGLISH)
txt_file = fullfile(fig_path,'Figure Legends.txt');
fp = fopen(txt_file,'w');
if -1==fp
    fprintf(1,'Cannot open: %s\r\n',txt_file);
    return;
end;
fprintf(fp,'Figure 01: Bar plot of **Quantified Peptide Count** per sample. Bar height indicates the number of peptides detected with Area > 0. Used as a QC metric; low counts suggest low sensitivity or sample quality issues.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 02: Box plot displaying the **Log_{10} transformed peptide intensity** (Area) distribution per sample. The median and quartile ranges assess consistency in sample loading and overall MS performance.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 03: **Principal Component Analysis (PCA)** of all samples based on the Z-scored Log2 ratio matrix. Displays the distance/clustering between samples. Replicates should cluster tightly.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 10: **Heatmaps of single PTM ratios** (Log2-scaled) and Z-scores (if auxiliary file exists). Displays the relative abundance summed across all peptides containing that specific modification.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 11: **Heatmaps of Log2-scaled ratios** for H3 and H4 peptides separately. Shows peptide-level abundance patterns, suitable for identifying modifications with high sample-to-sample variability.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 14: **Violin/Box plot** showing the distribution of Log2-transformed ratio values, separated by H3 and H4 histone classes for each sample. The shape represents data density/variability.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 15: **Scatter plot** comparing Log_{10}(Area) vs. Log_{2}-transformed Ratio (all peptides). A minimal correlation (R near 0) confirms that the calculated ratio is independent of absolute peptide abundance (i.e., robust normalization).\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 16: **Binary Heatmap** showing peptide detection (1=Detected, 0=Missing) across samples. Visual QC for data completeness and missing data patterns.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 17: **Bar plot** of the percentage of peptides in each sample with a Log2-Ratio corresponding to an extremely low relative abundance (Ratio < 1%%). High percentage may indicate low signal quality or high noise in quantification.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 18: **Box plot** comparing the global Log_{10}(Area) distribution between H3 and H4 peptides. Checks for systematic intensity biases between the two histone classes.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 19: **Scatter plot** comparing the average Log2 Ratio profile of H3 peptides vs. H4 peptides across samples. Correlates global modification trends between H3 and H4 classes.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 20: **Heatmap of Spearman Correlation (R-value)** between all sample pairs based on Log2 ratio profiles. Used as critical QC: high correlation (R > 0.9) is expected for technical or biological replicates.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 21: **Bar plot** showing the total theoretical peptide count for the H3 class and the H4 class.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 22: **Single Box plot** displaying the overall Log2 Ratio distribution across ALL quantified peptides and ALL samples combined.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 23: **Line plot** comparing the Median Log2 Ratio profile of H3 and H4 peptides across all samples. Visualizes global shifts in modification levels per class.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 24: **Box plot** showing the distribution of "Noise" Ratios (Log2 Ratio where |Ratio| < 1) for each sample. QC for system stability near the baseline signal.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 25: **PCA** of the samples using **H3 peptide ratios (Z-scored) ONLY**. Helps assess sample grouping driven exclusively by H3 modifications.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 26: **PCA** of the samples using **H4 peptide ratios (Z-scored) ONLY**. Helps assess sample grouping driven exclusively by H4 modifications.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 27: **Overlay Histogram** showing the frequency distribution of Log2 Ratios for each sample. Useful for comparing the overall distribution shapes and central tendency across runs.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 28: **Histogram of Peptide Completeness**. Shows the distribution of how many samples each peptide was successfully detected in (Area > 0). Highlights commonly missing or highly complete peptides.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 29: **Bar plot of Missing Peptide Count**. Shows the total number of peptides (Area = 0) per sample. Direct QC for sample specific data completeness/dropout rates.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 30: **Log_{2} Ratio Heatmap for H3(73-83) Peptides**. Shows the relative abundance of all modification isoforms quantified on the H3(73-83) region across samples. Allows for detailed visualization of this critical domain.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 31: **Box plot** of the **Ratio Distribution (Log_{2})** for H3(73-83) peptides only. Compares the central tendency and variability of these specific modifications across samples.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 32: **Scatter plot** comparing the **Average Log2 Ratio profile of H3(73-83)** peptides versus the **Average Log2 Ratio of ALL H3** peptides across samples. Measures if the modification dynamics of this specific domain align with the overall H3 modification profile.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 33: **Hierarchical Clustering Dendrogram**. Shows the clustering of samples based on the Euclidean distance of the Z-scored Log2 ratio matrix (using Ward linkage). Helps visualize the similarity structure between all samples, complementing the PCA results.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 34: **Paired Box Plot: H3 vs. H4 Ratio per Sample**. Displays side-by-side distributions of Log2-transformed ratios for H3 and H4 peptides within each sample. This allows direct comparison of overall modification levels and variability between histone classes in each specific condition.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 35: **Density Plot: Overall Log2 Ratio Distribution (H3 vs. H4)**. Provides a smoothed estimate of the probability density function for all H3 and H4 peptide ratios. This highlights differences in the global distribution shapes (e.g., shifts in median, skewness) between the two histone classes.\r\n');
fprintf(fp,'\r\n');
fprintf(fp,'Figure 36: **Spearman Correlation Heatmap: H3(73-83) Peptides vs. Other H3 Peptides**. Visualizes the pairwise Spearman correlation coefficient between each H3(73-83) modification isoform and all other H3 peptides (excluding H3(73-83) peptides themselves). High positive or negative correlations suggest co-regulation or inverse regulation, respectively, indicating potential functional associations.\r\n');
fclose(fp);
end % write_txt