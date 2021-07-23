source("00_bwb20200729.R")

# summarise by the new variable patient_day and by partition assignment.  This is the aligned partition assignment
cds_aligned_tbl_longsummary_1 <- cds_aligned_tbl %>% 
  group_by(VWsample,partition_assignment) %>% 
  summarise(n = n()) 
cds_aligned_tbl_longsummary_1

# join the total number of cells onto the patient_day count data and calculate percent
cds_aligned_clust_pct_1 <- left_join(cds_aligned_tbl_longsummary_1, aligned_cluster_counts) %>% 
  mutate(pct_aligned_cluster_total = n/total*100)
cds_aligned_clust_pct_1

# this is just incase you want to print out the wide form table
cds_aligned_clust_pct_1_wide <- cds_aligned_clust_pct_1 %>%
  select(VWsample,partition_assignment,pct_aligned_cluster_total) %>%
  pivot_wider(names_from = VWsample, values_from = pct_aligned_cluster_total,values_fill = 0)
cds_aligned_clust_pct_1_wide
cds_aligned_clust_pct_matrix_1 <- as.matrix(cds_aligned_clust_pct_1_wide[,-1])
rownames(cds_aligned_clust_pct_matrix_1) <- cds_aligned_clust_pct_1_wide %>% pull(partition_assignment)

#save.image.pigz("ficlatuzumab_memory_managed.RData", n.cores = 37)