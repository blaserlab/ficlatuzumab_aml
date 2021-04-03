source("00_bwb20200729.R")

cache_list_mm2 <- c(
  "cds_68k",
  "cds_aggregate_blasts_nonormal",
  "ref_68k_sample",
  "ref_68k",
  "cds_aligned_channel",
  "cds_blasts",
  "cds_blasts_aligned",
  "cds_blasts_aligned_CR",
  "cds_blasts_aligned_NR",
  "cds_blasts_CR",
  "cds_blasts_nonormal",
  "cds_blasts_NR",
  "cds_deseq",
  "cds_heatmap",
  "cds_heatmap_by_pt",
  "cds_hgf_monocle_regression",
  "cds_preprocess",
  "cds_s6_targets",
  "cds_deseq_list",
  "aggr_expr_mat_full",
  "cds_heatmap_counts_assigned_by_pt",
  "cds_heatmap_counts_assigned",
  "cds_heatmap_counts_tbl_by_pt",
  "cds_heatmap_counts_tbl",
  "mod7_deseq_list",
  "vcf_plot",
  "vcf_long",
  "cds_heatmap_counts_by_pt",
  "cds_heatmap_counts",
  "new_agg_cds_mod7",
  "new_agg_cds_mod3",
  "new_agg_cds_mod18",
  "new_agg_cds_mod1",
  "new_agg_cds"
  
  
)

save.pigz(list = cache_list_mm2, file = "ficlatuzumab_cache_list_mm2.RData",n.cores = 8)

rm(list = cache_list_mm2)

current_size_MB <- ls() %>% 
  sapply(. %>% 
           get %>% 
           object.size %>% 
           '/'(10^6)) 

cbind(current_size_MB, "MB") %>% as.data.frame

gc()
#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)