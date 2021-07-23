source("00_bwb20200729.R")
#start with cds blasts and calculate top 25 markers using builtin monocle functions

marker_test_res_by_pt <- top_markers(cds_blasts, group_cells_by="VWsample", reference_cells=1000, cores=39)

top_25_per_cluster_by_pt <- marker_test_res_by_pt %>% pull(gene_short_name) %>% unique()

cds_heatmap_by_pt <- cds_blasts[rowData(cds_blasts)$gene_short_name %in% top_25_per_cluster_by_pt,]

cds_heatmap_counts_by_pt <- normalized_counts(cds = cds_heatmap_by_pt, norm_method = "log", pseudocount = 1)

cds_heatmap_counts_tbl_by_pt <- broom::tidy(cds_heatmap_counts_by_pt)

cell_assignments_by_pt <- as_tibble(colData(cds_heatmap_by_pt)) %>% select(barcode_sample, VWsample)

cds_heatmap_counts_assigned_by_pt <- left_join(cds_heatmap_counts_tbl_by_pt, cell_assignments_by_pt, by = c("column" = "barcode_sample"))

heatmap_data_long_by_pt  <- cds_heatmap_counts_assigned_by_pt %>% group_by(row, VWsample) %>% summarise(mean_value = mean(value)) %>% ungroup()

heatmap_data_long_by_pt <-
  left_join(
    heatmap_data_long_by_pt,
    as_tibble(rowData(cds_blasts)) %>%
      select(id, gene_short_name),
    by = c("row" = "id")
  ) %>%
  select(gene_short_name, mean_value, VWsample)

heatmap_data_wide_by_pt <- heatmap_data_long_by_pt %>% pivot_wider(names_from = gene_short_name, values_from = mean_value)

heatmap_df_by_pt <- as.data.frame(heatmap_data_wide_by_pt)
heatmap_df_by_pt[is.na(heatmap_df_by_pt)] <- 0
rownames(heatmap_df_by_pt) <- heatmap_df_by_pt[,1]
heatmap_mat_by_pt <- as.matrix(heatmap_df_by_pt[,-1])
sum(is.na(heatmap_mat_by_pt))
#heatmap_mat_scaled <- scale(heatmap_mat, center = T, scale = T)
heatmap_mat_by_pt

# print out the list of specific markers
marker_test_res_by_pt %>% arrange(cell_group) %>% write_csv("data_out/pt_specific_blast_markers.csv")



#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)