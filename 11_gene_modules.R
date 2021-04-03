source("00_bwb20200729.R")

# # identify gene modules
# pr_graph_test_res <-
#   graph_test(cds,
#              neighbor_graph = "knn",
#              cores = 24,
#              verbose = TRUE)
# pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
# gene_module_df <- find_gene_modules(cds[pr_deg_ids, ], cores = 39, louvain_iter = 10)


gene_module_df_anno <- left_join(gene_module_df,as_tibble(rowData(cds))) %>% write_csv("data_out/all_modules.csv")

# gene modules according to response_day
cell_group_df2 <- tibble::tibble(cell=row.names(colData(cds_blasts)), cell_group=colData(cds_blasts)$response_day)# cds_blasts is comparable to cds with only the early, late and hladr partitions.  It includes the normal control.
agg_mat2 <- aggregate_gene_expression(cds_blasts, gene_module_df, cell_group_df2)
row.names(agg_mat2) = stringr::str_c("Module ", row.names(agg_mat2))
colnames(agg_mat2) = stringr::str_c(colnames(agg_mat2))

#set up the annotation dataframe
annotation_col2 <- data.frame(
  Day = factor(str_replace(colnames(agg_mat2),"^.*_","")),
  Response = factor(str_replace(colnames(agg_mat2),"_.*",""))
)
annotation_col2
rownames(annotation_col2) <-
    paste0(annotation_col2$Response, "_", annotation_col2$Day)
rownames(annotation_col2)[rownames(annotation_col2) == "Normal_day_Control"] <- "Normal"# gotta fix the outlier
annotation_col2

# make the colors easier to tell apart
clustercolors_purples<-brewer.pal(n = 5, name = "Purples")
clustercolors2<-c("darkgreen",clustercolors_purples)
names(clustercolors2)<- sort(factor(unique(annotation_col2$Day), levels = c("Control","0","1","2","3","42-44")))
annotation_colors2 <- list(
  Response = c(Normal ="white",CR = "cyan3",NR = "magenta3"),
  Day = clustercolors2
)
annotation_colors2

# print out the gene moduel components
if (!dir.exists("data_out/gene_modules")){
  dir.create("data_out/gene_modules")
}

modules <- 
  gene_module_df_anno %>% 
  pull(module) %>% 
  unique()

lapply(
  X = modules,
  FUN = function(x) {
    gene_module_df_anno %>%
      filter(module == x) %>%
      write_csv(paste0("data_out/gene_modules/module_", x, ".csv"))
  }
)


#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)
