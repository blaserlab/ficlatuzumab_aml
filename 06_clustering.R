source("00_bwb20200729.R")

## cluster cells
cds <- cluster_cells(cds,verbose = TRUE)
plot_cells(cds, color_cells_by = "partition",group_cells_by = "partition",group_label_size = 3)

## Find marker genes expressed by each cluster
marker_test_res_p <- top_markers(cds, group_cells_by="partition", reference_cells=1000, cores=39)

# annotate clusters
colData(cds)$partition <- partitions(cds)
colData(cds)$partition_assignment = as.character(partitions(cds))
colData(cds)$partition_assignment = dplyr::recode(colData(cds)$partition_assignment,#these all look the same, just in a different order
                                                "1"="Proliferative",#"T",#"CD16+ Monocytic",
                                                "2"="T",#"Myelomonocytic",#"CD14+ Monocytic",
                                                "3"="Myelomonocytic 1",#"Proliferative",#"Myelomonocytic",
                                                "4"="Myelomonocytic 2",#"CEBPD+ MCL1+ Monocytic",#"CEBPD+ MCL1+ Monocytic",
                                                "5"="CEBPD+ MCL1+ Monocytic",#"CD16+ Monocytic",#"T",
                                                "6"="CD16+ Monocytic",#"CD14+ Monocytic",#"Erythrocytic",
                                                "7"="CD14+ Monocytic",#"GATA2+ SOX4+",#"Proliferative",
                                                "8"="B",#"B",#"ENO1+",
                                                "9"="GATA2+ SOX4+",#"ENO1+",#"B/PC",
                                                "10"="ENO1+",#"Erythrocytic",#"B",
                                                "11"="Erythrocytic",#"B/PC",#"DC",
                                                "12"="B/PC",#"DC")#"GATA2+ SOX4+")
                                                "13"="DC")
colData(cds)$response = colData(cds)$VWsample
colData(cds)$response = dplyr::recode(colData(cds)$response,
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
colData(cds)$p_a_r<-paste0(colData(cds)$partition_assignment,"_",colData(cds)$response)

partition_lut<-unique(cbind(partitions(cds),colData(cds)$partition_assignment))
rownames(partition_lut) <- c()
colnames(partition_lut) <- c("cell_group", "partition_assignment")
partition_lut<-as.data.frame(partition_lut)

if (!dir.exists("data_out")){
  dir.create("data_out")
}

marker_test_res_p1<-left_join(marker_test_res_p,partition_lut,"cell_group") %>% write_csv("data_out/partition_markers.csv")

save.image.pigz("ficlatuzumab_memory_managed.RData", n.cores = 37)
