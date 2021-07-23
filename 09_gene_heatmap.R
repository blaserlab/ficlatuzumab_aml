source("00_bwb20200729.R")

colData(cds)$barcode_sample <- rownames(colData(cds))

top_25_per_cluster <- marker_test_res_p1_aligned %>% pull(gene_short_name) %>% unique() 

cds_heatmap <- cds[rowData(cds)$gene_short_name %in% top_25_per_cluster,]

cds_heatmap_counts <- normalized_counts(cds = cds_heatmap, norm_method = "log", pseudocount = 1)

cds_heatmap_counts_tbl <- broom::tidy(cds_heatmap_counts)

cell_assignments <- as_tibble(colData(cds_heatmap)) %>% select(barcode_sample,aligned_partition)

cds_heatmap_counts_assigned <- left_join(cds_heatmap_counts_tbl,cell_assignments, by = c("column" = "barcode_sample"))

heatmap_data_long  <- cds_heatmap_counts_assigned %>% group_by(row,aligned_partition) %>% summarise(mean_value = mean(value)) %>% ungroup()

heatmap_data_long <-
  left_join(
    heatmap_data_long,
    as_tibble(rowData(cds)) %>%
      select(id, gene_short_name),
    by = c("row" = "id")
  ) %>%
  select(gene_short_name, mean_value, aligned_partition)

heatmap_data_wide <- heatmap_data_long %>% pivot_wider(names_from = gene_short_name, values_from = mean_value)

heatmap_df <- as.data.frame(heatmap_data_wide)
heatmap_df[is.na(heatmap_df)] <- 0
rownames(heatmap_df) <- heatmap_df[,1]
heatmap_mat <- as.matrix(heatmap_df[,-1])
sum(is.na(heatmap_mat))
#heatmap_mat_scaled <- scale(heatmap_mat, center = T, scale = T)
heatmap_mat



save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)