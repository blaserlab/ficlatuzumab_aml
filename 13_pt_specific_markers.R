source("00_bwb20200729")

# want to find markers that are patient specific within aligned partitions

early_blast_tm_patient <- 
  top_markers(cds_trimmed_nonormal[!str_detect(rowData(cds_trimmed_nonormal)$gene_short_name,"^RP|^MT"),colData(cds_trimmed_nonormal)$aligned_partition == "Early"],
            group_cells_by = "VWsample",
            reference_cells = colData(cds_trimmed_nonormal) %>% as_tibble() %>% filter(aligned_partition == "Early") %>% pull(barcode_sample),
            cores = 39,
            genes_to_test_per_group = 50)
early_tm_patient_sig <- blast1_tm_patient %>% filter(marker_test_q_value<0.05) %>% write_csv("data_out/early_tm_patient_sig.csv")

save.image.pigz("ficlatuzumab_memory_managed.RData")
