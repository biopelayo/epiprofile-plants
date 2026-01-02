% Suponiendo que YA tienes preparados:
% info_auc, info_ratio_log, zscore_ratio, peptides, clean_names, sample_groups, group_boundaries, c_npeps, ratio_labels

style = epplot_defaults('FontSize',10,'DPI',600);   % estilo único y consistente
outdir = fullfile(layout_path,'heatmap_clustering'); if ~exist(outdir,'dir'), mkdir(outdir); end

% Derivados útiles
[idxH3, idxH4] = split_H3_H4(peptides);
data_H3 = info_ratio_log(idxH3,:);
data_H4 = info_ratio_log(idxH4,:);

% 01–03
fig01_bar_peptide_number(info_auc, clean_names, sample_groups, fullfile(outdir,'01_bar_peptide_number.pdf'), style);
fig02_boxplot_peptide_intensity(info_auc, clean_names, sample_groups, fullfile(outdir,'02_boxplot_peptide_intensity.pdf'), style);
fig03_pca_all(zscore_ratio, sample_groups, fullfile(outdir,'03_pca.pdf'), style);

% 14–20
fig14_violin_ratio_H3_H4(data_H3, data_H4, clean_names, sample_groups, fullfile(outdir,'14_violin_ratio_H3_H4.pdf'), style);
fig15_scatter_area_vs_ratio(info_auc, info_ratio_log, fullfile(outdir,'15_scatter_Area_vs_Ratio.pdf'), style);
fig16_heatmap_missing_values(info_auc, clean_names, sample_groups, fullfile(outdir,'16_heatmap_missing_values.pdf'), style);
fig17_low_ratio_percentage(info_ratio_log, clean_names, sample_groups, fullfile(outdir,'17_low_ratio_percentage.pdf'), style, 1);
fig18_boxplot_log10area_by_class(info_auc, idxH3, idxH4, fullfile(outdir,'18_boxplot_log10area_by_class.pdf'), style);
fig19_scatter_avg_ratio_H3_vs_H4(data_H3, data_H4, sample_groups, fullfile(outdir,'19_scatter_avg_ratio_H3_vs_H4.pdf'), style);
fig20_heatmap_sample_correlation(info_ratio_log, clean_names, sample_groups, fullfile(outdir,'20_heatmap_sample_correlation.pdf'), style);

% 21–26
fig21_bar_class_count(c_npeps, fullfile(outdir,'21_barplot_class_count.pdf'), style);
fig22_boxplot_ratio_all(info_ratio_log, fullfile(outdir,'22_boxplot_ratio_all_peptides.pdf'), style);
fig23_line_median_ratio_profile(data_H3, data_H4, clean_names, sample_groups, fullfile(outdir,'23_lineplot_median_ratio_profile.pdf'), style);
fig24_boxplot_ratio_noise(info_ratio_log, clean_names, sample_groups, fullfile(outdir,'24_boxplot_ratio_noise.pdf'), style, 1);
fig25_pca_H3_only(data_H3, sample_groups, fullfile(outdir,'25_pca_H3_only.pdf'), style);
fig26_pca_H4_only(data_H4, sample_groups, fullfile(outdir,'26_pca_H4_only.pdf'), style);

% 27–29
fig27_overlay_histograms(info_ratio_log, sample_groups, fullfile(outdir,'27_overlay_histogram_ratios.pdf'), style);
fig28_histogram_peptide_completeness(info_auc, fullfile(outdir,'28_histogram_peptide_completeness.pdf'), style);
fig29_bar_missing_count(info_auc, clean_names, sample_groups, fullfile(outdir,'29_barplot_sample_missing_count.pdf'), style);

% 30–32 (ejemplo con región H3_73_83 si la tienes)
ixRegion = contains(peptides,'H3_73_83');
if any(ixRegion)
    data_region = info_ratio_log(ixRegion,:);
    peptides_region = peptides(ixRegion);
    fig30_heatmap_region(data_region, peptides_region, clean_names, sample_groups, fullfile(outdir,'30_heatmap_ratio_H3_73_83.pdf'), style, ratio_labels);
    fig31_boxplot_region(data_region, clean_names, sample_groups, fullfile(outdir,'31_boxplot_ratio_H3_73_83.pdf'), style);
    fig32_scatter_region_vs_H3(data_region, data_H3, sample_groups, fullfile(outdir,'32_scatter_avg_ratio_H3_73_83_vs_H3.pdf'), style);
end

% 33–36
fig33_dendrogram_samples(zscore_ratio, clean_names, sample_groups, fullfile(outdir,'33_dendrogram_sample_clustering.pdf'), style);
% Para fig36 necesitas todos los H3 (excluyendo región si deseas)
fig36_heatmap_corr_region_vs_allH3(data_region, data_H3, peptides_region, peptides(idxH3), fullfile(outdir,'36_heatmap_corr_H3_73_83_vs_all_H3.pdf'), style);

% 37–40
fig37_bar_top_variable_peptides(info_ratio_log, peptides, 25, fullfile(outdir,'37_barplot_top_variable_peptides.pdf'), style);
fig38_heatmap_peptide_correlation(info_ratio_log, fullfile(outdir,'38_heatmap_peptide_correlation.pdf'), style);
fig39_MA_plot_peptides(info_auc, info_ratio_log, fullfile(outdir,'39_MA_plot_peptides.pdf'), style);
% fig40 requiere sample_batches (derívalo de tu phenodata)
% fig40_pca_batch_effect(zscore_ratio, sample_batches, fullfile(outdir,'40_pca_batch_effect.pdf'), style);
