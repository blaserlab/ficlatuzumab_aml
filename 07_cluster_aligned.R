source("00_packages_functions.R")

## cluster cells
cds_aligned <- cluster_cells(cds_aligned,verbose = TRUE)
plot_cells(cds_aligned, color_cells_by = "partition",group_cells_by = "partition",group_label_size = 3)

## Find marker genes expressed by each cluster
marker_test_res_p_aligned <- top_markers(cds_aligned, group_cells_by="partition", reference_cells=1000, cores=39)

# annotate clusters
colData(cds_aligned)$partition <- partitions(cds_aligned)
colData(cds_aligned)$partition_assignment = as.character(partitions(cds_aligned))
colData(cds_aligned)$partition_assignment = dplyr::recode(colData(cds_aligned)$partition_assignment,#these all look the same, just in a different order
                                                "1"="Early",#"T",#"CD16+ Monocytic",
                                                "2"="Late",#"Myelomonocytic",#"CD14+ Monocytic",
                                                "3"="T",#"Proliferative",#"Myelomonocytic",
                                                "4"="HLA-DR+",
                                                "5"="B",#"CEBPD+ MCL1+ Monocytic",#"CEBPD+ MCL1+ Monocytic",
                                                "6"="Erythrocytic",#"CD14+ Monocytic",#"Erythrocytic",
                                                "7"="B/PC",#"GATA2+ SOX4+",#"Proliferative",
                                                "8"="Dividing"#"B",#"ENO1+",
                                                )
colData(cds_aligned)$response = colData(cds_aligned)$VWsample
colData(cds_aligned)$response = dplyr::recode(colData(cds_aligned)$response,
                                                  "E01"="NR",
                                                  "E02"="NR",
                                                  "E03"="NR",
                                                  "E04"="CR",
                                                  "E05"="CR",
                                                  "E06"="CR",
                                                  "E07"="NR",
                                                  "E09"="CR",
                                                  "E10"="Normal",
                                                  "E11"="CR",
                                                  "E12"="CR",
                                                  "E13"="CR")

#create new coldata column partition_assigment_response (p_a_r)
colData(cds_aligned)$p_a_r<-paste0(colData(cds_aligned)$partition_assignment,"_",colData(cds_aligned)$response)

partition_lut<-unique(cbind(partitions(cds_aligned),colData(cds_aligned)$partition_assignment))
rownames(partition_lut) <- c()
colnames(partition_lut) <- c("cell_group", "partition_assignment")
partition_lut<-as.data.frame(partition_lut)

if (!dir.exists("data_out")){
  dir.create("data_out")
}

marker_test_res_p1_aligned<-left_join(marker_test_res_p_aligned,partition_lut,"cell_group") %>% 
  write_csv("data_out/partition_markers_aligned.csv")

# add aligned cluster assignments bac onto original cds
cds_aligned_tbl <- as_tibble(colData(cds_aligned)) 

colData(cds)$sanitycheck <- cds_aligned_tbl %>% pull(barcode_channel)
colData(cds)$aligned_partition <- cds_aligned_tbl %>% pull(partition_assignment)

#sanity check
sum(colData(cds)$barcode_channel != colData(cds)$sanitycheck)#looks good
colData(cds)$sanitycheck <- NULL

plot_cells(cds, color_cells_by = "aligned_partition")



save.image.pigz("ficlatuzumab_memory_managed.RData", n.cores = 37)
