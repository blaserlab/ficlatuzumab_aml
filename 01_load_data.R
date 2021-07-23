source("00_bwb20200729.R")

cyan1 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-22-17cyan1",barcode_filtered = TRUE);colData(cyan1)$channel <- "cyan1"
cyan2 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-22-17cyan2",barcode_filtered = TRUE);colData(cyan2)$channel <- "cyan2"

blue1 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-29-17blue1",barcode_filtered = TRUE);colData(blue1)$channel <- "blue1"
blue2 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-29-17blue2",barcode_filtered = TRUE);colData(blue2)$channel <- "blue2"

green1 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-30-17green1",barcode_filtered = TRUE);colData(green1)$channel <- "green1"
green2 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-30-17green2",barcode_filtered = TRUE);colData(green2)$channel <- "green2"

red1 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_9-5-17red1",barcode_filtered = TRUE);colData(red1)$channel <- "red1"
red2 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_9-5-17red2",barcode_filtered = TRUE);colData(red2)$channel <- "red2"

yellow1 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_9-5-17yellow1",barcode_filtered = TRUE);colData(yellow1)$channel <- "yellow1"
yellow2 <- load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/ngs_archive/aml/VW_9-5-17yellow2",barcode_filtered = TRUE);colData(yellow2)$channel <- "yellow2"

cds_list <- list(cyan1,cyan2,blue1,blue2,green1,green2,red1,red2,yellow1,yellow2)

cds0 <- combine_cds(cds_list = cds_list, keep_all_genes = TRUE,cell_names_unique = FALSE)

colData(cds0)$barcode_channel <- paste0(colData(cds0)$barcode,"_",colData(cds0)$channel)
colData(cds0)$barcode_sample <- rownames(colData(cds0))

#save.image.pigz("ficlatuzumab_freemuxlet.RData")