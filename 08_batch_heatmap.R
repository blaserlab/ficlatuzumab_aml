source("00_bwb20200729.R")
# goal here is to find out how many cells from each patient and each timepoint are in the aligned clusters
# then make a heatmap to indicate what percent of all cells in each aligned cluster is from each patient-sample
# then use this to justify cleaning up non-representative samples
cols_to_add <- left_join(as_tibble(colData(cds_aligned)),fml_metadata_full, by = c("barcode_channel" = "barcode_channel")) %>% select(barcode_channel, treatment_day1 = treatment_day1.y)
colData(cds_aligned)$sanitycheck  <- cols_to_add$barcode_channel
colData(cds_aligned)$treatment_day1 <- cols_to_add$treatment_day1
#sanity check
sum(colData(cds_aligned)$barcode_channel != colData(cds_aligned)$sanitycheck)#looks good
colData(cds_aligned)$sanitycheck <- NULL

# make a composite variable of patient-treatment day
cds_aligned_tbl <- as_tibble(colData(cds_aligned)) %>%
  mutate(patient_day = paste0(VWsample,"_",treatment_day1))
cds_aligned_tbl
# summarise by the new variable patient_day and by partition assignment.  This is the aligned partition assignment
cds_aligned_tbl_longsummary <- cds_aligned_tbl %>% 
  group_by(patient_day,partition_assignment) %>% 
  summarise(n = n()) 
cds_aligned_tbl_longsummary

# get the total counts per partition assignment in order to calculate the percent of total
aligned_cluster_counts <- cds_aligned_tbl %>%
  group_by(partition_assignment) %>%
  summarise(total = n()) 
aligned_cluster_counts

# get the max counts per partition assignment
max_per_partition <- cds_aligned_tbl_longsummary %>%
  group_by(partition_assignment) %>%
  summarise(max_n = max(n)) 
max_per_partition

# join the total number of cells onto the patient_day count data and calculate percent
cds_aligned_clust_pct <- left_join(cds_aligned_tbl_longsummary, aligned_cluster_counts) %>% 
  mutate(pct_aligned_cluster_total = n/total*100) 
cds_aligned_clust_pct

# now join on the max counts per partition

cds_aligned_clust_pct_max <- left_join(cds_aligned_clust_pct,max_per_partition) %>%
  mutate(pct_max = n/max_n*100)
cds_aligned_clust_pct_max

# this is just incase you want to print out the wide form table
cds_aligned_clust_pct_wide <- cds_aligned_clust_pct %>%
  select(patient_day,partition_assignment,pct_aligned_cluster_total) %>%
  pivot_wider(names_from = patient_day, values_from = pct_aligned_cluster_total,values_fill = 0)
cds_aligned_clust_pct_wide
cds_aligned_clust_pct_matrix <- as.matrix(cds_aligned_clust_pct_wide[,-1])
rownames(cds_aligned_clust_pct_matrix) <- cds_aligned_clust_pct_wide %>% pull(partition_assignment)

col_ord <- hclust(dist(t(cds_aligned_clust_pct_matrix), method = "euclidean"), method = "ward.D" )$order

#make the ordered long dataset
cds_aligned_clust_pct_ord <- cds_aligned_clust_pct
cds_aligned_clust_pct_ord$patient_day <- factor(cds_aligned_clust_pct_ord$patient_day,
                                                levels = rownames(t(cds_aligned_clust_pct_matrix))[col_ord])

# this adds the new composite variable back onto the aligned cds
colData(cds_aligned)$patient_day  <- cds_aligned_tbl %>% pull(patient_day)




#save.image.pigz("ficlatuzumab_memory_managed.RData", n.cores = 37)