source("00_bwb20200729.R")

#load the data and get it into useable form


cytof_population_means <-
  read_csv("cytof_data/cytof_population_means.csv") %>%
  mutate(population = str_replace(population, " -.*$", "")) %>%
  mutate(sample_id = str_replace(sample_id, "\\.fcs.*$", "")) %>%
  mutate(sample_id = str_sub(sample_id, 18, -1)) %>%
  mutate(treatment = str_extract(sample_id, "_.{2,4}$") %>% str_sub(2, -1)) %>%
  mutate(patient = str_extract(sample_id, "^\\w{3,5}_") %>% str_sub(1, -2)) %>%
  mutate(day = str_extract(sample_id, "_.._") %>% str_sub(2, -2)) %>%
  select(
    c(
      -`175Lu_pS6`,
      -`144Nd_pPLCg2`,
      -`158Gd_pStat3`,
      -`171Yb_pERK`,
      -`156Gd_p-p38`,
      -`165Ho_pCREB`,
      -`146Nd_IgD`,
      -`143Nd_bio_pMET`,
      -`167Er_total_pMET`,
      -`157Gd_CD24`,
      -`169Tm_CD45RA`,
      -`153Eu_STAT1`,
      -`152Sm_Akt`,
      -`164Dy_Ikba`,
      -`150Nd_STAT5`
    )
  ) %>%
  filter(treatment == "UN") %>%
  filter(population != "Singlet") %>%
  filter(sample_id %notin% c("HIMC1_UN", "HIMC2_UN", "HIMC3_UN"))
cytof_population_means

cytof_cell_groups <- cytof_population_means %>% select(population) %>% unique()

cytof_samples <- cytof_population_means %>% pull(sample_id) %>% unique()

#extract the cell ids for cd34+ cd33+ t and b cells from the scrnaseq data
## for t and b used teh aligned clusters

enlist_barcodes <- function(sample_ids, barcode_table,i) {
  sample_ids <- sample_ids[i]
  barcodes <- barcode_table %>%
    filter(sample_id == sample_ids) %>%
    pull(barcode_sample)
  return(barcodes)
}

cd33_cell_barcodes <- 
  plot_cells(cds = cds, genes = "CD33")[["plot_env"]][["data_df"]] %>% 
  as_tibble() %>%
  select(barcode_sample,VWsample,treatment_day1,CD33_value = value) %>%
  mutate(sample_id = paste0(VWsample,"_D",treatment_day1,"_UN")) %>%
  filter(!is.na(CD33_value))

cd33_barcode_list <- lapply(X = seq_along(cytof_samples), FUN = enlist_barcodes, sample_ids = cytof_samples, barcode_table = cd33_cell_barcodes)
names(cd33_barcode_list) <- cytof_samples

cd34_cell_barcodes <- 
  plot_cells(cds = cds, genes = "CD34")[["plot_env"]][["data_df"]] %>% 
  as_tibble() %>%
  select(barcode_sample,VWsample,treatment_day1,CD34_value = value) %>%
  mutate(sample_id = paste0(VWsample,"_D",treatment_day1,"_UN")) %>%
  filter(!is.na(CD34_value))

cd34_barcode_list <- lapply(X = seq_along(cytof_samples), FUN = enlist_barcodes, sample_ids = cytof_samples, barcode_table = cd34_cell_barcodes)
names(cd34_barcode_list) <- cytof_samples

Tcell_barcodes <- 
  as_tibble(colData(cds)) %>%
  filter(aligned_partition == "T") %>%
  select(barcode_sample, VWsample, treatment_day1) %>%
  mutate(sample_id = paste0(VWsample,"_D",treatment_day1,"_UN")) 
Tcell_barcodes  

Tcell_barcode_list <- lapply(X = seq_along(cytof_samples), FUN = enlist_barcodes, sample_ids = cytof_samples, barcode_table = Tcell_barcodes)
names(Tcell_barcode_list) <- cytof_samples

Bcell_barcodes <- 
  as_tibble(colData(cds)) %>%
  filter(aligned_partition == "B") %>%
  select(barcode_sample, VWsample, treatment_day1) %>%
  mutate(sample_id = paste0(VWsample,"_D",treatment_day1,"_UN")) 
Bcell_barcodes  

Bcell_barcode_list <- lapply(X = seq_along(cytof_samples), FUN = enlist_barcodes, sample_ids = cytof_samples, barcode_table = Bcell_barcodes)
names(Bcell_barcode_list) <- cytof_samples

# get size_factor normalized counts for a list of marker genes and average for the cell groups made above.
cytof_genes <- c(
  ENSG00000156738 = "MS4A1",
  ENSG00000185291 = "IL3RA",
  ENSG00000174059 = "CD34",
  ENSG00000173762 = "CD7",
  ENSG00000177455 = "CD19",
  ENSG00000139193 = "CD27",
  ENSG00000140678 = "ITGAX",
  ENSG00000170458 = "CD14",
  ENSG00000203747 = "FCGR3A",
  ENSG00000153563 = "CD8A",
  ENSG00000198851 = "CD3E",
  ENSG00000105383 = "CD33",
  ENSG00000010610 = "CD4",
  ENSG00000081237 = "PTPRC",
  ENSG00000004468 = "CD38",
  ENSG00000168685 = "IL7R",
  ENSG00000134460 = "IL2RA",
  ENSG00000139193 = "CD27",
  ENSG00000196126 = "HLA-DRB1"
)



get_cytof_counts <- function(cds, barcodes, genes, population,i) {
  name <- names(barcodes[i])
  barcodes <- barcodes[[i]]
  cds <- cds[rowData(cds)$gene_short_name %in% genes,colData(cds)$barcode_sample %in% barcodes]
  counts <- as.matrix(normalized_counts(cds = cds))
  row.names(counts) <- recode(row.names(counts),!!!genes)
  mean_counts <- rowMeans(counts,na.rm = T)
  mean_counts <- as_tibble(t(mean_counts))
  sample_name <- tibble(sample_id = name)
  return_data <- bind_cols(sample_name,mean_counts)
  population <- tibble(population = population)
  return_data <- bind_cols(population,return_data) %>% mutate(population_sample = paste0(population,"_",sample_id))
  return(return_data)
}

cd33_sc_data <- bind_rows(
  lapply(
    X = seq_along(cd33_barcode_list),
    FUN = possibly(get_cytof_counts, otherwise = NULL, quiet = F),
    cds = cds,
    barcodes = cd33_barcode_list,
    genes = cytof_genes,
    population = "CD33+"
  )
)


cd34_sc_data <- bind_rows(
  lapply(
    X = seq_along(cd34_barcode_list),
    FUN = possibly(get_cytof_counts, otherwise = NULL, quiet = F),
    cds = cds,
    barcodes = cd34_barcode_list,
    genes = cytof_genes,
    population = "CD34+"
  )
)

Tcells_sc_data <- bind_rows(
  lapply(
    X = seq_along(Tcell_barcode_list),
    FUN = possibly(get_cytof_counts, otherwise = NULL, quiet = F),
    cds = cds,
    barcodes = Tcell_barcode_list,
    genes = cytof_genes,
    population = "Tcells"
  )
)


Bcell_sc_data <- bind_rows(
  lapply(
    X = seq_along(Bcell_barcode_list),
    FUN = possibly(get_cytof_counts, otherwise = NULL, quiet = F),
    cds = cds,
    barcodes = Bcell_barcode_list,
    genes = cytof_genes,
    population = "Bcell"
  )
)

sc_data_full <- bind_rows(
  cd33_sc_data,
  cd34_sc_data,
  Tcells_sc_data,
  Bcell_sc_data
) %>% select(-population,-sample_id)
sc_data_full <- sc_data_full[complete.cases(sc_data_full),]

#add a colum onto the cytof data for matching
cytof_population_means <- cytof_population_means %>% mutate(population_sample = paste0(population,"_",sample_id))

#join the cytof and scrnaseq data together
cytof_scrnaseq_means <- left_join(cytof_population_means,sc_data_full)
cytof_scrnaseq_means <- cytof_scrnaseq_means[complete.cases(cytof_scrnaseq_means),]

make_cor_data <- function(data, cell_populations = NULL, i = NULL) {
  if (!is.null(cell_populations)) {
    cell_populations <- cell_populations[i]
    data <- data %>%
      filter(population == cell_populations) %>%
      select(-population,-sample_id,-treatment,-patient,-day)
  } else {
    data <- data %>%
      select(-population,-sample_id,-treatment,-patient,-day)
  }
  data <- as.data.frame(data)
  row.names(data) <- data$population_sample
  data$population_sample <- NULL
  data <- data[complete.cases(data), ]
  
  #remove columns with no counts
  data <- data[,colSums(data)>0]
  
  #make the correlation matrix
  cor_res <- cor(data)
  
  return(cor_res)
}


cytof_scrnaseq_cor_data <- make_cor_data(data = cytof_scrnaseq_means)

#make the annotation dataframes
# 
# col_anno <- data.frame(
#   row.names = rownames(cytof_scrnaseq_cor_data),
#   "Celltype" = recode(
#     rownames(cytof_scrnaseq_cor_data),
#     "153Eu_STAT1" = "CD33",
#     "147Sm_CD20" = "B",
#     "151Eu_CD123__IL-3R" = "CD34",
#     "154Sm_CD45" = "Pan-Hematopoietic",
#     "148Nd_CD34" = "CD34",
#     "141Pr_CD7" = "T",
#     "142Nd_CD19" = "B",
#     "155Gd_CD27" = "T",
#     "159Tb_CD11c" = "CD33",
#     "160Gd_CD14" = "CD33",
#     "166Er_CD16" = "CD33",
#     "168Er_CD8a" = "T",
#     "170Er_CD3" = "T",
#     "174Yb_HLA-DR" = "B",
#     "163Dy_CD33" = "CD33",
#     "FCGR3A" = "CD33",
#     "PTPRC" = "Pan-Hematopoietic",
#     "CD34" = "CD34",
#     "CD8A" = "T",
#     "STAT1" = "CD33",
#     "CD14" = "CD33",
#     "HLA-DRB1" = "B",
#     "IL3RA" = "CD34",
#     "MS4A1" = "B",
#     "CD3E" = "T",
#     "CD27" = "T",
#     "CD19" = "B",
#     "ITGAX" = "CD33",
#     "CD7" = "T",
#     "CD33" = "CD33"
#   ),
#   "Assay" = recode(
#     rownames(cytof_scrnaseq_cor_data),
#     "153Eu_STAT1" = "CyTOF",
#     "147Sm_CD20" = "CyTOF",
#     "151Eu_CD123__IL-3R" = "CyTOF",
#     "154Sm_CD45" = "CyTOF",
#     "148Nd_CD34" = "CyTOF",
#     "141Pr_CD7" = "CyTOF",
#     "142Nd_CD19" = "CyTOF",
#     "155Gd_CD27" = "CyTOF",
#     "159Tb_CD11c" = "CyTOF",
#     "160Gd_CD14" = "CyTOF",
#     "166Er_CD16" = "CyTOF",
#     "168Er_CD8a" = "CyTOF",
#     "170Er_CD3" = "CyTOF",
#     "174Yb_HLA-DR" = "CyTOF",
#     "163Dy_CD33" = "CyTOF",
#     "FCGR3A" = "scRNA-seq",
#     "PTPRC" = "scRNA-seq",
#     "CD34" = "scRNA-seq",
#     "CD8A" = "scRNA-seq",
#     "STAT1" = "scRNA-seq",
#     "CD14" = "scRNA-seq",
#     "HLA-DRB1" = "scRNA-seq",
#     "IL3RA" = "scRNA-seq",
#     "MS4A1" = "scRNA-seq",
#     "CD3E" = "scRNA-seq",
#     "CD27" = "scRNA-seq",
#     "CD19" = "scRNA-seq",
#     "ITGAX" = "scRNA-seq",
#     "CD7" = "scRNA-seq",
#     "CD33" = "scRNA-seq"
#   )
# )
# 
# 
# col_anno
# 
# # annotation colors
# anno_colors = list(
#   Assay = c(CyTOF = "black", `scRNA-seq` = "white"),
#   Celltype = c(B = "#e12be3", CD33 = "#f8e200", CD34 = "#015179", `Pan-Hematopoietic` = "#760000", `T` = "#00c190")
# )
# 
# 

#plot the diagonally symmetric correlation matrix
cor_heatmap <-
  pheatmap::pheatmap(
    cytof_scrnaseq_cor_data,
    angle_col = 90,
    treeheight_row = 10,
    treeheight_col = 10,
    show_colnames = TRUE#,
    #annotation_col = col_anno,
    #annotation_colors = anno_colors
  )
cor_heatmap <- as.ggplot(cor_heatmap)
cor_heatmap

# replot the cor heatmap with same values but with cytof on one side and scrnaseq on the other
cytof_scrnaseq_cor_data_alt1 <- cytof_scrnaseq_cor_data[rownames(cytof_scrnaseq_cor_data) %>% str_which("^[:digit:]"),
                                                        rownames(cytof_scrnaseq_cor_data) %>% str_which("^[:alpha:]")]
cytof_scrnaseq_cor_data_alt1

#make the annotation dataframes

# col_anno_alt1 <- data.frame(
#   row.names = colnames(cytof_scrnaseq_cor_data_alt1),
#   "Celltype" = recode(
#     colnames(cytof_scrnaseq_cor_data_alt1),
#     "FCGR3A" = "CD33",
#     "CD34" = "CD34",
#     "CD8A" = "T",
#     "STAT1" = "CD33",
#     "CD14" = "CD33",
#     "IL3RA" = "CD34",
#     "MS4A1" = "B",
#     "CD3E" = "T",
#     "CD27" = "T",
#     "CD19" = "B",
#     "ITGAX" = "CD33",
#     "CD7" = "T",
#     "CD33" = "CD33"
#   )
# )
# 
# row_anno_alt1 <- data.frame(
#   row.names = rownames(cytof_scrnaseq_cor_data_alt1),
#   "CyTOF" = recode(
#     rownames(cytof_scrnaseq_cor_data_alt1),
#     "166Er_CD16" = "CD33",
#     "148Nd_CD34" = "CD34",
#     "168Er_CD8a" = "T",
#     "153Eu_STAT1" = "CD33",
#     "160Gd_CD14" = "CD33",
#     "151Eu_CD123__IL-3R" = "CD34",
#     "147Sm_CD20" = "B",
#     "170Er_CD3" = "T",
#     "155Gd_CD27" = "T",
#     "142Nd_CD19" = "B",
#     "159Tb_CD11c" = "CD33",
#     "141Pr_CD7" = "T",
#     "163Dy_CD33" = "CD33"
#   )
# )
# 
# anno_colors_alt1 = list(
#   Celltype = c(B = "#e12be3", CD33 = "#f8e200", CD34 = "#015179", `T` = "#00c190"),
#   CyTOF = c(B = "#e12be3", CD33 = "#f8e200", CD34 = "#015179", `T` = "#00c190")
# )
# 
# 
# 

cytof_scrnaseq_cor_data_alt2 <- cytof_scrnaseq_cor_data_alt1
rownames(cytof_scrnaseq_cor_data_alt2) <- str_replace(rownames(cytof_scrnaseq_cor_data_alt1),"^.*_","")


# column_ha <- HeatmapAnnotation(scrna_markers = anno_block(gp = gpar(fill = c("#fff362", "#ff8fdf", "#61c958",  "#a7ffeb")),
#                                                           labels = c("CD34","CD33","B","T")),
#                                which = "column")
# row_ha <- HeatmapAnnotation(cytof_markers = anno_block(gp = gpar(fill = c("#fff362", "#ff8fdf", "#61c958",  "#a7ffeb")),
#                                                        labels = c("CD34","CD33","B","T"),
#                                                        labels_rot = 0),
#                             which = "row")
set.seed(123)
cytof_scrnaseq <- as.ggplot(Heatmap(cytof_scrnaseq_cor_data_alt2,
        col = colorRampPalette(rev(brewer.pal(
          n = 7, name = "RdYlBu"
        )))(100),
        row_title = "CyTOF Markers",
        row_title_side = "left",
        row_dend_side = "right",
        row_names_side = "left",
        column_title = "scRNA-seq Markers",
        column_title_side = "bottom",column_names_rot = 45,
        # top_annotation = column_ha,
        # column_km = 4,
        # right_annotation = row_ha,
        # row_km = 4,
        # show_parent_dend_line = F,
        heatmap_legend_param = list(title = "Pearson\nCorrelation",
                                    legend_height = unit(4, "cm"),
                                    title_gp = gpar(fontsize = 11)))
)

cytof_scrnaseq

# named and ordered vector of cytof and scrna markers
# cytof = scrna
cytof_scrna_pairs <- c(
  "CD45" = "PTPRC",
  "CD25" = "IL2RA",
  "IL-7R" = "IL7R",
  "CD3" = "CD3E",
  "CD7" = "CD7",
  "CD8a" = "CD8A",
  "CD27" = "CD27",
  "CD4" = "CD4",
  "CD34" = "CD34",
  "IL-3R" = "IL3RA",
  "CD14" = "CD14",
  "CD11c" = "ITGAX",
  "CD33" = "CD33",
  "CD16" = "FCGR3A",
  "HLA-DR" = "HLA-DRB1",
  "CD20" = "MS4A1",
  "CD19" = "CD19",
  "CD38" = "CD38"
)

cytof_scrnaseq_supervised <- as.ggplot(Heatmap(cytof_scrnaseq_cor_data_alt2,
        col = colorRampPalette(rev(brewer.pal(
          n = 7, name = "RdYlBu"
        )))(100),
        row_title = "CyTOF Markers",
        row_title_side = "left",
        row_dend_side = "right",
        row_names_side = "left",
        column_title = "scRNA-seq Markers",
        column_title_side = "bottom",
        column_names_rot = 45, row_order = names(cytof_scrna_pairs), column_order = rev(unname(cytof_scrna_pairs)) ,
        # top_annotation = column_ha,
        # column_km = 4,
        # right_annotation = row_ha,
        # row_km = 4,
        # show_parent_dend_line = F,
        heatmap_legend_param = list(title = "Pearson\nCorrelation",
                                    legend_height = unit(4, "cm"),
                                    title_gp = gpar(fontsize = 11)))
)

cytof_scrnaseq_supervised

dim(cytof_scrnaseq_cor_data_alt2)


# cytof population proportions
# make a table for interleaved bar chart
cytof_scrna_cellcounts <- inner_join(
  read_csv("cytof_data/ravi_analysis/celltype.freq.csv") %>%
    pivot_longer(-X1, names_to = "population", values_to = "cellcount") %>%
    rename(cytof_sample = X1) %>%
    mutate(patient = str_replace(cytof_sample, "_.*", "")) %>%
    mutate(day = str_extract(cytof_sample, "D[:digit:]")) %>%
    filter(patient != "HIMC") %>%
    mutate(
      merged_population = recode(
        population,
        "CD24_myeloid" = "Late",
        "NK" = "T/NK",
        "pcDC" = "HLA-DR+",
        "T" = "T/NK"
      )
    ) %>%
    group_by(patient, day, merged_population) %>%
    summarise(cytof_cellcount = sum(cellcount)) %>%
    mutate(cytof_frac = cytof_cellcount / sum(cytof_cellcount)),
  
  colData(cds) %>%
    as_tibble() %>%
    select(VWsample, treatment_day1, aligned_partition) %>%
    transmute(
      patient = VWsample,
      day = paste0("D", treatment_day1),
      population = aligned_partition
    ) %>%
    filter(day %in% c("D0", "D1", "D2")) %>%
    filter(
      patient %in% c("E01", "E03", "E04", "E05", "E07", "E09", "E12", "E13")
    ) %>%
    filter(population %in% c("Early", "Late", "HLA-DR+", "B", "T")) %>%
    mutate(merged_population = recode(population, "T" = "T/NK")) %>%
    group_by(patient, day, merged_population) %>%
    summarise(scrna_cellcount = n()) %>%
    mutate(scrna_frac = scrna_cellcount / sum(scrna_cellcount))
) %>%
  mutate(pdmp = paste0(patient,"_",day,"_",merged_population)) %>%
  ungroup() %>%
  select(pdmp,cytof_frac,scrna_frac,merged_population, patient, day)



cytof_scrna_cellcounts
cor(tbl_to_matrix(cytof_scrna_cellcounts %>% select(-merged_population, -patient, -day)))

cellcount_fit <- lm(log10(cytof_frac)~log10(scrna_frac), data = cytof_scrna_cellcounts)
summary(cellcount_fit)
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page
plot(cellcount_fit)

#save.image.pigz("ficlatuzumab_memory_managed.RData")
