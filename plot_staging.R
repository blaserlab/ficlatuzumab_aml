source("00_bwb20200729.R")
theme_set(theme_cowplot(font_size = 11))

# rename the B/PC partition to PC in all relevant cds's
colData(cds)$aligned_partition <- recode(colData(cds)$aligned_partition, "B/PC" = "PC")
colData(cds_aligned)$partition_assignment <- recode(colData(cds_aligned)$partition_assignment, "B/PC" = "PC")

#qc plots####------------------------------------------------------------------------------------
channels_pct_counts_mito_cds5 <- ggplot(as_tibble(colData(cds5)), aes(x = channel, y = pct_counts_MITO)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 83))

channels_nfeature_rna_cds5 <- ggplot(as_tibble(colData(cds5)), aes(x = channel, y = log10(nFeature_RNA))) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(2.3,3.8))

# make figures of percent mito and n features by best guess assignment
vwsample_pct_counts_mito_cds5 <- ggplot(as_tibble(colData(cds5)), aes(x = VWsample, y = pct_counts_MITO)) +
  geom_violin() +
  coord_cartesian(ylim = c(0, 83))

vwsample_nfeature_rna_cds5 <- ggplot(as_tibble(colData(cds5)), aes(x = VWsample, y = log10(nFeature_RNA))) +
  geom_violin() +
  coord_cartesian(ylim = c(2.3,3.8))

# make figures of percent mito and n features by channel
channels_pct_counts_mito_cds6 <- ggplot(as_tibble(colData(cds6)), aes(x = channel, y = pct_counts_MITO)) +
  geom_violin() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0,83))

channels_nfeature_rna_cds6 <- ggplot(as_tibble(colData(cds6)), aes(x = channel, y = log10(nFeature_RNA))) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(2.3,3.8))

# make figures of percent mito and n features by best guess assignment
vwsample_pct_counts_mito_cds6 <- ggplot(as_tibble(colData(cds6)), aes(x = VWsample, y = pct_counts_MITO)) +
  geom_violin() +
  coord_cartesian(ylim = c(0,83))

vwsample_nfeature_rna_cds6 <- ggplot(as_tibble(colData(cds6)), aes(x = VWsample, y = log10(nFeature_RNA))) +
  geom_violin() +
  coord_cartesian(ylim = c(2.3,3.8))

#survival plots####-------------------------------------------------------------------------------
pfs_plot <- ggsurvplot(
  pfs_fit,
  data = pfs_data,
  risk.table = T,
  pval = F,
  break.time.by = 300,
  risk.table.y.text = F,
  add.all = T,
  ggtheme = theme_cowplot(font_size = 18)
)
pfs_plot$plot <- pfs_plot$plot + labs(x = "Days", title = "Progression Free Survival") + theme(plot.title = element_text(hjust = 0.5))
pfs_plot$table <- pfs_plot$table + labs(x = "Days") + theme(plot.title = element_text(face = "plain")) 
pfs_plot_fixed <- plot_grid(pfs_plot$plot,pfs_plot$table, align = "v", axis = "l", nrow = 2, rel_heights = c(2,1))
pfs_plot_fixed


os_plot <- ggsurvplot(
  os_fit,
  data = os_data,
  risk.table = T,
  pval = T,
  break.time.by = 300,
  risk.table.y.text = F,
  add.all = T
)

os_plot_fixed <- fix_surv_plot(input_plot = os_plot, fontsize = 11, title = NULL, xtitle = "Days")
os_plot_fixed

hsct_os_plot <- ggsurvplot(
  hsct_os_fit,
  data = hsct_os_data,
  risk.table = T,
  pval = T,
  break.time.by = 300,
  risk.table.y.text = F,
  add.all = T
)

hsct_os_plot_fixed <- fix_surv_plot(input_plot = hsct_os_plot, fontsize = 11, title = NULL, xtitle = "Days")
hsct_os_plot_fixed

#categorical plots####----------------------------------------------------------------------------------------
brewer_palette <- brewer.pal(n=8, name = "Set2")

batch_cluster_palette <- c(
  "Early" = brewer_palette[1],
  "T" = brewer_palette[2],
  "Late" =  brewer_palette[3],
  "HLA-DR+" =  brewer_palette[4],
  "B" = brewer_palette[5],
  "Erythrocytic" = brewer_palette[6],
  "B/PC" =  brewer_palette[7],
  "Dividing" = brewer_palette[8],
  "PC" = brewer_palette[7]
)

brewer_palette_patient <- brewer.pal(n = 12, name = "Paired")

patient_palette <- c(
  "E01" = brewer_palette_patient[1],
  "E02" = brewer_palette_patient[2],
  "E03" = brewer_palette_patient[3],
  "E04" = brewer_palette_patient[4],
  "E05" = brewer_palette_patient[5],
  "E06" = brewer_palette_patient[6],
  "E07" = brewer_palette_patient[7],
  "E09" = brewer_palette_patient[8],
  "E10" = brewer_palette_patient[9],
  "E11" = brewer_palette_patient[10],
  "E12" = brewer_palette_patient[11],
  "E13" = brewer_palette_patient[12]
)

# plot aligned data by vwsample
var_plot_patient <- custom_variable_plot(cds = cds_aligned, palette = patient_palette, var = "VWsample", foreground_alpha = 0.4)
var_plot_patient

# plot aligned data by aligned partition
cds_cp_plot_aligned <- custom_cp_plot_min(cds = cds_aligned, 
                                          cp = "partition", 
                                          alpha = 0.2, 
                                          palette = batch_cluster_palette,
                                          overwrite_labels = T)
cds_cp_plot_aligned

# plot unaligned cds by aligned partition
aligned_partition_plot <- custom_cp_plot(
  cds = cds, 
  cp = "aligned_partition",
  var = "aligned_partition", 
  alpha = 0.2,
  overwrite_labels = T,
  palette = batch_cluster_palette) 
aligned_partition_plot

#plot unaligned cds by patient
patient_plot_unaligned <- custom_variable_plot(
  cds = cds,
  var = "VWsample",
  foreground_alpha = 0.2,
  palette = patient_palette
  
)

patient_plot_unaligned

#gene heatmap####--------------------------------------------------------------------------------
heatmap_dotplot_markers <- c("MS4A1","CD38","CD79A","XBP1","CDK6","CD99","CD34","S100A9","LYZ","TUBB","HBA2","HLA-DRB1","CD74","CD3D")

#make the heatmap with gene labels
aligned_partition_heatmap_labeled <- pheatmap(mat = heatmap_mat, scale = "column", show_colnames = T)
# extract the labels
heatmap_markers0 <- aligned_partition_heatmap_labeled[["tree_col"]][["labels"]]
# make a vector of the ones we want
heatmap_markers_only <- heatmap_dotplot_markers
# heatmap_markers_only <- c("MS4A1","CD38","CD79A","XBP1","CDK6","CD99","CD34","S100A9","LYZ","CDKN3","TUBB","HBA2","HLA-DRB1","CD74","CD3D")
# replace the ones we don't want with ""
heatmap_markers_show <- ifelse(heatmap_markers0 %in% heatmap_markers_only, heatmap_markers0, "")
# make a new heatmap with the labels we want
aligned_partition_heatmap_labeled_select_rows <- pheatmap(mat = heatmap_mat, 
                                                          scale = "column", 
                                                          show_colnames = T, 
                                                          labels_col = heatmap_markers_show,
                                                          treeheight_col = 10,
                                                          treeheight_row = 0,
                                                          fontsize = 8,
                                                          angle_col = 90)


heatmap_markers_vicki <- c(heatmap_dotplot_markers,"BANK1","TNFRSF17","AURKB","SOX4","FECH","CD1C","CST3","MNDA","MALAT1","MYBL2","RAB13","EIF3E","CD44","FOS","CD79A","IGJ","IGLL5","CENPU","PTTG1","GMNN","ZWINT","BIRC5","GYPA")

gene_heatmap <- as.ggplot(
  Heatmap(
    scale(heatmap_mat),
    col = rev(brewer.pal(n = 7, name = "RdYlBu")),
    row_names_gp = gpar(fontsize = 8),
    show_column_names = F,
    column_dend_side = "bottom",
    row_dend_width = unit(5, "mm"),
    row_dend_gp = gpar(lwd = 0.5),
    column_dend_height = unit(5, "mm"),
    column_dend_gp = gpar(lwd = 0.5),
    heatmap_legend_param = list(
      title = "Expression",
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    ),
    top_annotation = columnAnnotation(link = anno_mark(
      at = which(colnames(heatmap_mat) %in% heatmap_markers_vicki),
      labels = colnames(heatmap_mat)[colnames(heatmap_mat) %in% heatmap_markers_vicki],
      labels_gp = gpar(fontsize = 8)
    ))
  )
)

#patient heatmap####--------------------------------------------------------------------------------
#make the heatmap with gene labels
aligned_partition_heatmap_labeled_by_pt <- pheatmap(mat = heatmap_mat_by_pt, scale = "column", show_colnames = T)
# extract the labels
heatmap_markers0_by_pt <- aligned_partition_heatmap_labeled_by_pt[["tree_col"]][["labels"]]
# make a vector of the ones we want
heatmap_markers_only_by_pt <- c("HLA-DRB1","ISG15","NPM1","CDK6","JUN","RPS6","CD14","IL1B","IER3","SOX4","CD164")
# replace the ones we don't want with ""
heatmap_markers_show_by_pt <- ifelse(heatmap_markers0_by_pt %in% heatmap_markers_only_by_pt, heatmap_markers0_by_pt, "")
# make a new heatmap with the labels we want
aligned_partition_heatmap_labeled_select_rows_by_pt <- pheatmap(mat = heatmap_mat_by_pt, 
                                                          scale = "column", 
                                                          show_colnames = T, 
                                                          labels_col = heatmap_markers_show_by_pt,
                                                          treeheight_col = 10,
                                                          treeheight_row = 10,
                                                          fontsize = 8,
                                                          angle_col = 90,)


gene_heatmap_by_pt <- ggplotify::as.ggplot(aligned_partition_heatmap_labeled_select_rows_by_pt)+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 50, unit = "pt"))
save_plot(gene_heatmap_by_pt, filename = "figs_out/gene_heatmap_by_pt.pdf", base_height = 5, base_width = 6)

#cluster membership figs####---------------------------------------------------------------------
#plot the percents in geom_tile

# cluster_heatmap <- ggplot(data = cds_aligned_clust_pct, aes(x = patient_day, y = partition_assignment, fill = pct_aligned_cluster_total)) +
#   geom_tile() +
#   scale_fill_gradient(low = "grey90", high = "red", na.value = "transparent") +
#   labs(x = "Patient-Timepoint", y = "Partition",fill = "Percent of\nPartition") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# cluster_heatmap

#plot with pheatmap
mat_colnames <- colnames(cds_aligned_clust_pct_matrix)
col_anno <- data.frame(Patient = str_sub(mat_colnames,1,3),
                       Day = str_replace(mat_colnames,"[:alnum:]*_",""))
rownames(col_anno) <- mat_colnames

col_colors <- list(Patient = patient_palette, Day = annotation_colors2[["Day"]])

cluster_pheatmap <- ggplotify::as.ggplot(
  pheatmap(
    log10(cds_aligned_clust_pct_matrix + 1),
    color = colorRampPalette(c("transparent", "grey80", "red"))(50),
    treeheight_row = 0,
    treeheight_col = 0,
    annotation_col = col_anno,
    annotation_colors = col_colors,
    show_colnames = F,fontsize = 11
  )
)

#plot with all timepoints together
mat_colnames_1 <- colnames(cds_aligned_clust_pct_matrix_1)
col_anno_1 <- data.frame(Patient = mat_colnames_1)
rownames(col_anno_1) <- mat_colnames_1

col_colors_1 <- list(Patient = patient_palette)

cluster_pheatmap_1 <- ggplotify::as.ggplot(
  pheatmap(
    log10(cds_aligned_clust_pct_matrix_1 + 1),
    color = colorRampPalette(c("transparent", "grey80", "red"))(50),
    treeheight_row = 0,
    treeheight_col = 0,
    show_colnames = T,
    fontsize = 11,
    angle_col = 0,cluster_cols = F
  )
)


cluster_anno <-
  HeatmapAnnotation(
    Response = recode(
      colnames(cds_aligned_clust_pct_matrix_1),
      "E01" = "NR",
      "E02" = "NR",
      "E03" = "NR",
      "E04" = "CR",
      "E05" = "CR",
      "E06" = "CR",
      "E07" = "NR",
      "E09" = "CR",
      "E10" = "Normal",
      "E11" = "CR",
      "E12" = "CR",
      "E13" = "CR"
    ),
    col = list(Response = c(
      "Normal" = "white",
      "CR" = "cyan3",
      "NR" = "magenta3"
    )),
    annotation_legend_param = list(title = NULL,ncol = 1,title_gp = gpar(fontsize = 11),border = "grey80"),
    gp = gpar(col = "grey80"),
    show_annotation_name = F
    
  )


cluster_complex_heatmap <- as.ggplot(
  Heatmap(
    matrix = log10(cds_aligned_clust_pct_matrix_1 + 1),
    col = colorRampPalette(c("transparent", "grey80", "red"))(50),
    cluster_rows = F,
    cluster_columns = F,
    column_names_rot = 0,
    column_names_centered = T,
    border = "grey80",
    rect_gp = gpar(col = "grey80"),
    heatmap_legend_param = list(title = NULL, legend_height = unit(2, "cm")),
    bottom_annotation = cluster_anno,
    column_title = "Patient",
    column_title_gp = gpar(fontsize = 11),
    column_title_side = "bottom",
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 9)
  )
)
cluster_complex_heatmap




#plotted with fill color as percent of maximum per partition
cluster_heatmap_percent_max <- ggplot(data = cds_aligned_clust_pct_max, aes(x = patient_day, y = partition_assignment, fill = pct_max)) + 
  geom_tile() + 
  scale_fill_gradient(low = "grey90", high = "red", na.value = "transparent") + 
  labs(x = "Patient-Timepoint", y = "Partition", fill = "Percent of Max") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
cluster_heatmap_percent_max


# make a umap of the aligned cds, faceted in a grid by timepoint and patient, ie sample, colored by partition assignnment from the aligned cds

aligned_umap_sample <- custom_cp_plot(cds_aligned, cp = "partition",overwrite_labels = F, legend_pos = "none") + 
  facet_grid(rows = vars(VWsample), cols = vars(treatment_day1)) + 
  theme(strip.background = element_blank()) 
aligned_umap_sample

#dot plot showing genes expressed in partitions####------------------------------------------------------
#dot_plot_genes <- c("S100A9","CD3D","CD34","HLA-DRA","MS4A1","CDK6") 
#dot_plot_gene_alias <- c("S100A9","CD3D","CD34","HLA-DRA","CD20","CDK6")

dot_plot_genes <- c(heatmap_markers_only,"POU2AF1", "CD4", "CD8A", "GZMB", "CD83", "CD86", "LGALS6", "NCAM1", "KLRD1", "CD40")
dot_plot_gene_alias <- str_replace(dot_plot_genes, "MS4A1", "MS4A1")

names(dot_plot_gene_alias) <- dot_plot_genes
dot_plot_gene_alias
marker_gene_plot <- custom_gene_dotplot(cds = cds,
                    markers = dot_plot_genes,group_cells_by = "aligned_partition",
                    gene_alias_vector = dot_plot_gene_alias,
                    ordering_type = "cluster_row_col") +
  scale_color_gradient(low = "grey90", high =  "red")+
  labs(x = "Cluster", color = "log expression")
marker_gene_plot


plot_genes_by_group(cds = cds,markers = c(dot_plot_genes,"MKI67"), group_cells_by = "aligned_partition")

#AML blast markers####-----------------------------------------------------------------------------------
#represent as dot plots on cds_trimmed
aml_marker_plot <- custom_gene_dotplot(
  cds = cds_trimmed[,colData(cds_trimmed)$aligned_partition %in% c("T","Early","Late","HLA-DR+")],
  markers = c(heatmap_dotplot_markers[heatmap_dotplot_markers %notin% c("HBA2", "XBP1", "CDK6","TUBB","MS4A1","CD79A")],"GATA2","CD33","CD14"),
  # markers = c("CD34", "KIT", "GATA2", "CD14", "CD33", "ITGAM","CD19","MS4A1","CD3D","CD7"),
  # markers = c(dot_plot_genes,"CD34", "KIT", "GATA2", "CD14", "CD33", "ITGAM","CD19","MS4A1","CD3D","CD7"),
  group_cells_by = "aligned_partition",
  gene_alias_vector = dot_plot_gene_alias,
  ordering_type = "cluster_row_col",
  max.size = 4
) +
  scale_color_gradient(low = "grey90", high =  "#DC0000") +
  guides(size = guide_legend(order = 2,title = "Proportion")) +
  labs(color = "Expression",x = NULL,y = NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  # theme(legend.position = "bottom", legend.box = "vertical", legend.justification = "center",legend.key.width = unit(0.4,"in"))
aml_marker_plot

#HGF expression visualized####------------------------------------------------------
#fractionated dotplot
hgf_dot_plot <- 
  ggplot(hgf_expr_data_plot %>% filter(aligned_partition != "B"), 
       aes(x = response, 
           y = aligned_partition, 
           size = prop_exprs, 
           color = mean_exprs)) + 
  geom_point() +
  facet_wrap(facets = vars(treatment_day1),nrow = 1)+
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")+
  # theme(axis.line.y = element_blank())+
  # theme(axis.ticks.y = element_blank())+
  theme(axis.text.y = element_text(angle = 90, size = 8, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(strip.background = element_blank())+
  scale_color_gradient(low = "grey90", high = "#3C5488")+
  scale_size_area(max_size = 4)+
  guides(size = guide_legend(order = 2)) 
hgf_dot_plot

# umap expression for supplement
colData(cds_aligned)$response <- factor(colData(cds_aligned)$response, levels = c("CR","NR","Normal"))
hgf_umap_plot <- custom_variable_plot(cds = cds_aligned,#[,colData(cds_aligned)$response != "Normal"], 
                                      var = "normalized_hgf", 
                                      log10transform = TRUE) +
  facet_grid(cols = vars(treatment_day1), rows = vars(response)) +
  theme(strip.background = element_blank()) +
  labs(color = "Expression", title = NULL)
hgf_umap_plot

# CSF1R faceted by day and response only in blasts####---------------------------------------------------

csf1r_dot_plot <- 
  ggplot(csf1r_expr_data_plot %>% filter(aligned_partition != "B"), 
       aes(x = response, 
           y = aligned_partition, 
           size = prop_exprs, 
           color = mean_exprs)) + 
  geom_point() +
  facet_wrap(facets = vars(treatment_day1), nrow = 1) + 
  labs(x = NULL, y = NULL, size = "Proportion", color = "Expression")+
  # theme(axis.line.y = element_blank())+
  # theme(axis.ticks.y = element_blank())+
  theme(axis.text.y = element_text(angle = 90, size = 8, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 8)) + 
  theme(strip.background = element_blank()) +
  scale_color_gradient(low = "grey90", high = "darkgreen")+
  scale_size_area(max_size = 4)+
  guides(size = guide_legend(order = 2))
csf1r_dot_plot


# umap expression for supplement

csf1r_umap_plot <- custom_variable_plot(cds = cds_aligned,#[,colData(cds_aligned)$response != "Normal"], 
                                        var = "normalized_csf1",
                                        log10transform = TRUE) +
  facet_grid(cols = vars(treatment_day1), rows = vars(response)) +
  theme(strip.background = element_blank()) +
  labs(color = "Expression", title = NULL)
csf1r_umap_plot

#gene module heatmaps####--------------------------------------------------------------------------
#gene modules byt he composite factor response_day
gene_module_rd <-
  as.ggplot(
    pheatmap::pheatmap(
      agg_mat2,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "column",
      clustering_method = "ward.D2",
      fontsize = 8,
      treeheight_row = 10,
      treeheight_col = 10,
      annotation_col = annotation_col2,
      annotation_colors = annotation_colors2,
      show_colnames = F
    )
  )

# now make the gene_module_rd heatmap but only at d0 and d42
gene_module_rd_0_42_44 <- 
  pheatmap::pheatmap(
    agg_mat3,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "column",
    clustering_method = "ward.D2",
    fontsize = 10,
    treeheight_row = 10,
    treeheight_col = 10,
    annotation_col = annotation_col3,
    annotation_colors = annotation_colors3,
    show_colnames = F
  )

#gene modules by partition_assignment_response
gene_module_p_a_r <-
  pheatmap::pheatmap(
    agg_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "column",
    clustering_method = "ward.D2",
    fontsize = 10,#5,
    treeheight_row = 10,
    treeheight_col = 10,
    annotation_col = annotation_col1,
    annotation_colors = annotation_colors1,
    show_colnames = F
  )

#violin genes from the top 10 list from regression analysis (response + patient)-----------------
#d0
d0_genes_to_plot <- response_patient_model_coefs_d0_sig %>% top_n(10,-q_value) %>% pull(gene_short_name)

custom_violin_plot(cds = cds_mod14_d0, variable = "response", genes_to_plot = d0_genes_to_plot)

#violin plots from mod1 aggregate####----------------------------------------------------------
# mod1_agg_data <- plot_genes_violin(cds_subset = new_agg_cds_mod1[rowData(new_agg_cds_mod1)$gene_short_name == "mod1",])[["data"]]

mod1_agg_violins <- ggplot(mod1_agg_data, aes(x = response, y = expression, fill = response)) + 
  geom_jitter(size = 0.025, alpha = 0.4, color = "grey80", fill = "transparent", shape = 21,stroke = 0.2) +
  geom_violin(scale = "area", alpha = 0.4, color = "black") +
  scale_fill_manual(values = c("cyan3","magenta3"),guide = FALSE) +
  facet_wrap(facets = vars(treatment_day1), nrow = 1) +
  theme(strip.background = element_blank()) +
  labs(x = NULL, y = "Expression")
# mod1_agg_violins

#violin plots from mod18 aggregate####----------------------------------------------------------
# mod3_agg_data <- plot_genes_violin(cds_subset = new_agg_cds_mod3[rowData(new_agg_cds_mod3)$gene_short_name == "mod3",])[["data"]]

mod3_agg_violins <- ggplot(mod3_agg_data, aes(x = response, y = expression, fill = response)) + 
  geom_jitter(size = 0.025, alpha = 0.4, color = "grey80", fill = "transparent", shape = 21, stroke = 0.2) +
  geom_violin(scale = "area", alpha = 0.4, color = "black") +
  scale_fill_manual(values = c("cyan3","magenta3"),guide = FALSE) +
  facet_wrap(facets = vars(treatment_day1), nrow = 1) +
  theme(strip.background = element_blank())+
  labs(x = NULL, y = "Expression")
# mod3_agg_violins

#violin plots from mod7 aggregate####----------------------------------------------------------
# mod7_agg_data <- plot_genes_violin(cds_subset = new_agg_cds_mod7[rowData(new_agg_cds_mod7)$gene_short_name == "mod7",])[["data"]]

mod7_agg_violins <- ggplot(mod7_agg_data, aes(x = response, y = expression, fill = response)) + 
  geom_jitter(size = 0.025, alpha = 0.4, color = "grey80", fill = "transparent", shape = 21, stroke = 0.2) +
  geom_violin(scale = "area", alpha = 0.4, color = "black") +
  scale_fill_manual(values = c("cyan3","magenta3"),guide = FALSE) +
  facet_wrap(facets = vars(treatment_day1), nrow = 1) +
  theme(strip.background = element_blank())+
  labs(x = NULL, y = "Expression")
# mod7_agg_violins

#violin plots from mod18 aggregate####----------------------------------------------------------
# mod18_agg_data <- plot_genes_violin(cds_subset = new_agg_cds_mod18[rowData(new_agg_cds_mod18)$gene_short_name == "mod18",])[["data"]]

mod18_agg_violins <- ggplot(mod18_agg_data, aes(x = response, y = expression, fill = response)) + 
  geom_jitter(size = 0.025, alpha = 0.4, color = "grey80", fill = "transparent", shape = 21, stroke = 0.2) +
  geom_violin(scale = "area", alpha = 0.4, color = "black") +
  scale_fill_manual(values = c("cyan3","magenta3"),guide = FALSE) +
  facet_wrap(facets = vars(treatment_day1), nrow = 1) +
  theme(strip.background = element_blank()) +
  labs(x = NULL, y = "Expression")
# mod18_agg_violins
#volcano plot day 0####--------------------------------------------------------------------------------------------------------
genes_to_highlight_0 <- c("HGF","CSF1R")
volcano_data_0 <- pseudobulk_res_by_day[[1]][[2]] %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_0, gene_short_name, ""))

#volcano_data_0
## Volcano plot
volcano_plot_0 <- ggplot(volcano_data_0,aes(x = log2FoldChange, y = -log10(padj), colour = threshold, fill = threshold, label = text_label)) +
  geom_point(shape = 21) +
  geom_text_repel(color = "black", point.padding = 0)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80","red"))+
  scale_fill_manual(values = c("transparent","red"))
volcano_plot_0


#volcano plot day 42-44####--------------------------------------------------------------------------------------------------------
genes_to_highlight_42_44 <- c("HGF","CSF1R")


volcano_data_42_44 <- pseudobulk_res_by_day[[5]][[2]] %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_42_44, gene_short_name, ""))

## Volcano plot
volcano_plot_42_44 <- ggplot(volcano_data_42_44, aes(x = log2FoldChange, y = -log10(padj), colour = threshold, fill = threshold, label = text_label)) +
  geom_point(shape = 21) +
  geom_text_repel(color = "black", point.padding = 0)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80","red"))+
  scale_fill_manual(values = c("transparent","red"))
volcano_plot_42_44


#revigo bubbles####----------------------------------------------------------------------------------------------------------------
# preview the revigo data and make the plots
# module 1
# read_csv("revigo_data/REVIGO_module1.csv") %>% View()

revigo_mod1 <- revigo_bubbles(revigo.data = read_csv("revigo_data/REVIGO_module1.csv") %>% mutate(description = recode(description, "positive regulation of biological process" = "positive regulation of\nbiological process")), 
               dispensability = 0.2, 
               pval = -10,
               plot_text_size = 2
               ) + 
  theme_cowplot(font_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))
# revigo_mod1


# module 3
# read_csv("revigo_data/REVIGO_module3.csv") %>% View()

revigo_mod3 <- revigo_bubbles(revigo.data = read_csv("revigo_data/REVIGO_module3.csv") %>% mutate(description = recode(description, "cell-cell adhesion via plasma-membrane adhesion molecules" = "cell-cell adhesion via plasma-\nmembrane adhesion molecules")), 
               dispensability = 0, 
               pval = -3,
               plot_text_size = 2
               ) + 
  theme_cowplot(font_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))
# revigo_mod3


# module 7
# read_csv("revigo_data/REVIGO_module7.csv") %>% View()

revigo_mod7 <- revigo_bubbles(revigo.data = read_csv("revigo_data/REVIGO_module7.csv") %>% 
                                mutate(description = recode(description, 
                                                            "regulation of type I interferon production" = "regulation of type I\ninterferon production",
                                                            "negative regulation of viral life cycle" = "negative regulation of\nviral life cycle")), 
               dispensability = 0.15, 
               pval = -5.9,
               plot_text_size = 2
               ) + 
  theme_cowplot(font_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))
# revigo_mod7


# module 18
# read_csv("revigo_data/REVIGO_module18.csv") %>% View()

revigo_mod18 <- revigo_bubbles(revigo.data = read_csv("revigo_data/REVIGO_module18.csv") %>% 
                                mutate(description = recode(description, 
                                                            "nuclear-transcribed mRNA catabolic process, nonsense-mediated decay" = "nonsense-mediated decay",
                                                            "SRP-dependent cotranslational protein targeting to membrane" = "SRP-dependent cotranslational\nprotein targeting to membrane")), 
               dispensability = 0.19, 
               pval = -140,
               plot_text_size = 2
               ) + 
  theme_cowplot(font_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))
# revigo_mod18

#module bar charts####-------------------------------------------------------------------------------------------------------------
gorilla_mod7 <- read_delim("gorilla_data/gorilla_mod7.txt",delim = "\t") %>%
  mutate(neg_log10_q = -1*log10(`FDR q-value`)) %>%
  arrange(desc(neg_log10_q)) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  mutate(Description = factor(Description, levels = rev(levels(Description)))) %>%
  filter(Description %in% c("defense response to virus",
                            "response to other organism",
                            "defense response to other organism",
                            "type I interferon signaling pathway",
                            "cytokine-mediated signaling pathway")) %>%
  mutate(Description_rev = recode(Description, 
                                  "defense response to virus" = "Defense response to virus ",
                                  "response to other organism" = "Response to other organism ",
                                  "defense response to other organism" = "Defense response ",
                                  "type I interferon signaling pathway" = "Type I Interferon signaling ",
                                  "cytokine-mediated signaling pathway" = "Cytokine signaling "))

gorilla_mod7_barplot <- ggplot(gorilla_mod7, aes(x = Description_rev, y = neg_log10_q))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.1), width = 0.8, color = "black", fill = "white")+
  geom_text(aes(label = Description_rev), hjust = "right", size = 3)+
  coord_flip()+
  labs(y = expression(-log[10]*" q value"), x = "GO Term")+
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())
gorilla_mod7_barplot


gorilla_mod18 <- read_delim("gorilla_data/gorilla_mod18.txt",delim = "\t") %>%
  mutate(neg_log10_q = -1*log10(`FDR q-value`)) %>%
  arrange(desc(neg_log10_q)) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  mutate(Description = factor(Description, levels = rev(levels(Description)))) %>%
  filter(Description %in% c("SRP-dependent cotranslational protein targeting to membrane",
                            "nuclear-transcribed mRNA catabolic process, nonsense-mediated decay",
                            "viral transcription",
                            "protein localization to endoplasmic reticulum",
                            "translational initiation")) %>%
  mutate(Description_rev = recode(Description, 
                                  "SRP-dependent cotranslational protein targeting to membrane" = "SRP-dependent targeting to membrane ",
                                  "nuclear-transcribed mRNA catabolic process, nonsense-mediated decay" = "Nonsense-mediated decay ",
                                  "viral transcription" = "Viral transcription ",
                                  "protein localization to endoplasmic reticulum" = "Ptn Localization to ER ",
                                  "translational initiation" = "Translational initiation "
                                  ))
  
gorilla_mod18_barplot <- ggplot(gorilla_mod18, aes(x = Description_rev, y = neg_log10_q))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.1), width = 0.8, color = "black", fill = "white")+
  geom_text(aes(label = Description_rev), hjust = "right", size = 3)+
  coord_flip()+
  labs(y = expression(-log[10]*" q value"), x = "GO Term")+
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())
gorilla_mod18_barplot


# put seurat labels on cells #####
colData(cds_aligned)$seurat_assignment <- recode(colData(cds)$predicted.celltype.l2, "Platelet" = "")

seurat_plot <- custom_cp_plot(cds_aligned, cp = "seurat_assignment",alpha = 0.2)

# show hladr cluster

custom_variable_plot(cds_aligned, var = "partition_assignment", value_to_highlight = "HLA-DR+") + facet_wrap(facets = vars(VWsample)) 

# hgf targets and repressed genes

hgf_induced_dotplot

hgf_repressed_dotplot

# xbp1 irf4 TUBB

custom_gene_dotplot(cds = cds, markers = c("XBP1","IRF4","TUBB"), group_cells_by = "aligned_partition")

# celltype proportions between cytof and scrna

cytof_scrna_cellcount_plot <- ggplot(cytof_scrna_cellcounts, mapping = aes(x = log10(scrna_frac), y = log10(cytof_frac))) +
  geom_smooth(method='lm') +
  geom_jitter(shape = 21, size = 2.5, stroke = 0.5, mapping = aes(color = merged_population, fill = merged_population)) +
  scale_fill_manual(values = alpha(brewer.pal(n = 5, name = "Dark2"), alpha = 0.3)) +
  scale_color_manual(values = brewer.pal(n = 5, name = "Dark2")) +
  labs(x = expression(paste(log[10], "(scRNAseq Fraction)")),
       y = expression(paste(log[10], "(CyTOF Fraction)")),
       fill = "Population", 
       color = "Population")
cytof_scrna_cellcount_plot


#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)
