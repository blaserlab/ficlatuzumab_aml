source("00_bwb20200729.R")

cds_deseq_list <- list(
  cds_blasts_nonormal[,colData(cds_blasts_nonormal)$treatment_day1 == "0"],
  cds_blasts_nonormal[,colData(cds_blasts_nonormal)$treatment_day1 == "1"],
  cds_blasts_nonormal[,colData(cds_blasts_nonormal)$treatment_day1 == "2"],
  cds_blasts_nonormal[,colData(cds_blasts_nonormal)$treatment_day1 == "3"],
  cds_blasts_nonormal[,colData(cds_blasts_nonormal)$treatment_day1 == "42-44"]
)

# make a function to perform deseq on pseudobulk data from input cds
# all comparisons are NR vs CR.  Positive l2fc means up in NR.

nr_cr_pseudobulk <- function(cds_deseq) {
  
  groups <-
    colData(cds_deseq) %>% as_tibble() %>% select(VWsample, response)
  
  # get the aggregate counts
  aggregate_counts <-
    aggregate.Matrix(t(counts(cds_deseq)), groupings = groups, fun = "sum")
  counts_matrix <- as.matrix(t(aggregate_counts))
  
  # make the coldata for deseq
  samples <- colnames(counts_matrix)
  response <- str_replace(colnames(counts_matrix), "[:alnum:]*_", "")
  
  coldata <- data.frame(response, row.names = samples)
  stopifnot(all(rownames(coldata) == colnames(counts_matrix)))
  
  # make the deseq object
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = coldata,
                                design = ~ response)
  
  # do the thing
  dds <- DESeq(dds)
  res <- results(dds)
  result <-
    as.data.frame(res)  %>% rownames_to_column(var = "id") %>% as_tibble() %>%
    left_join(., as_tibble(rowData(cds_deseq)[, c("id", "gene_short_name")]))
  
  #qc
  rld <- rlog(dds, blind = TRUE)
  pca_plot <- DESeq2::plotPCA(rld, intgroup = "response")
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Plot heatmap
  heatmap <- pheatmap(rld_cor, annotation = coldata[, c("response"), drop = F])
  dispersion <- plotDispEsts(dds)
  
  return_list <- list(res@elementMetadata@listData[["description"]],
                      result, 
                      pca_plot, 
                      heatmap, 
                      dispersion
                      )
  return(return_list)
}

pseudobulk_res_by_day <- lapply(
  X = cds_deseq_list, 
  FUN = possibly(nr_cr_pseudobulk, otherwise = NULL)
)

pseudobulk_res_by_day[[1]][[2]] %>% arrange(padj) %>% write_csv("data_out/pseudobulk_day0.csv")
pseudobulk_res_by_day[[5]][[2]] %>% arrange(padj) %>% write_csv("data_out/pseudobulk_day42-44.csv")

names(pseudobulk_res_by_day)<- c("Day0","Day1","Day2","Day3","Day42-44")
