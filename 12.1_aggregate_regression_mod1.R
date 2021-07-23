source("00_bwb20200729.R")

cds_aggregate_blasts_nonormal<-cds_blasts_nonormal

get_gene_lists <- function(cds_dictionary, gene_names, list_names, i = NULL){
  if (!is.null(i)) {
    gene_names <- gene_names[[i]]
    list_name <- list_names[i]
    
  } else {
    list_name <- list_names
  }
  gene_tibble<- as_tibble(rowData(cds_dictionary)) %>% 
    filter(gene_short_name %in% gene_names) %>% 
    select(id,gene_short_name)
  gene_list<-list(gene_tibble)
  names(gene_list) <- list_name
 return(gene_list) 
}

mod1genes<-gene_module_df_anno %>% filter(module==1) %>% pull(gene_short_name)
mod3genes<-gene_module_df_anno %>% filter(module==3) %>% pull(gene_short_name)
mod7genes<-gene_module_df_anno %>% filter(module==7) %>% pull(gene_short_name)
mod11genes<-gene_module_df_anno %>% filter(module==11) %>% pull(gene_short_name)
mod18genes<-gene_module_df_anno %>% filter(module==18) %>% pull(gene_short_name)


#need to run these as single genes because there is an error when adding more than one composite score onto the new cds
gene_name_list<-list(mod1genes)
list_name_list<-c("mod1")

module_genes <- lapply(X = seq_along(list_name_list),
                       FUN = get_gene_lists,
                       cds_dictionary = cds_aggregate_blasts_nonormal,
                       gene_names = gene_name_list,
                       list_names = list_name_list)

get_agg_counts <- function(cds, gene_list, row_name_list, i = NULL) {
  colData(cds)$barcode_sample <- rownames(colData(cds))
  if (!is.null(i)) {
    gene_metadata <- gene_list[[i]][[1]]
    row_name_list <- row_name_list[i]
  } else {
    gene_metadata <- gene_list[[1]]
  }
  counts_tbl <- broom::tidy(exprs(cds)) %>%
    filter(row %in% gene_metadata$id) %>%
    group_by(column) %>%
    summarise(counts = sum(value)) %>%
    left_join(as_tibble(colData(cds)), ., by = c("barcode_sample" = "column")) %>%
    select(barcode_sample, counts)
  agg_expr_mat <-
    as.matrix(pivot_wider(counts_tbl, names_from = barcode_sample, values_from = counts))
  row.names(agg_expr_mat) <- row_name_list
  return(agg_expr_mat)
}

agg_counts_list <- lapply(X = seq_along(module_genes),
                          FUN = get_agg_counts,
                          cds = cds_aggregate_blasts_nonormal,
                          gene_list = module_genes,
                          row_name_list = list_name_list)

agg_expr_matrix <- do.call(rbind, agg_counts_list)

#sanity check
sum(colnames(agg_expr_matrix) != rownames(colData(cds_aggregate_blasts_nonormal)))
sum(colnames(agg_expr_matrix) != colnames(exprs(cds_aggregate_blasts_nonormal)))
aggr_expr_mat_full<-rbind2(exprs(cds_aggregate_blasts_nonormal),agg_expr_matrix)

agg_gene_metadata<-data.frame("id" = list_name_list, "gene_short_name" = list_name_list)
agg_gene_metadata
row.names(agg_gene_metadata) <- list_name_list
aggr_gene_metadata_full<-rbind2(rowData(cds_aggregate_blasts_nonormal)[,1:2],agg_gene_metadata)
class(aggr_expr_mat_full)

new_agg_cds<-new_cell_data_set(expression_data = aggr_expr_mat_full, 
                               cell_metadata = colData(cds_aggregate_blasts_nonormal), 
                               gene_metadata = aggr_gene_metadata_full)

new_agg_cds_mod1 <- new_agg_cds

mod1_monocle_regression <- lapply(
  X = seq_along(treatment_days),
  FUN = monocle_regression_func,
  cds = new_agg_cds_mod1,
  gene_or_genes = "mod1",
  stratification_variable = "treatment_day1",
  stratification_value = sort(treatment_days),
  form = "~response"
)
mod1_monocle_regression



#save.image.pigz("ficlatuzumab_memory_managed.RData")
