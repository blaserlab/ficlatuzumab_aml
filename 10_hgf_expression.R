source("00_bwb20200729.R")

cols_to_add <- left_join(as_tibble(colData(cds)),fml_metadata_full, by = c("barcode_channel" = "barcode_channel")) %>% select(barcode_channel, treatment_day1 = treatment_day1.y)
colData(cds)$sanitycheck  <- cols_to_add$barcode_channel
colData(cds)$treatment_day1 <- cols_to_add$treatment_day1
#sanity check
sum(colData(cds_aligned)$barcode_channel != colData(cds_aligned)$sanitycheck)#looks good
colData(cds_aligned)$sanitycheck <- NULL


#make new cds just for plotting violins and hgf in dotplots

cds_trimmed <- cds[,colData(cds)$aligned_partition %in% c("Early","Late","HLA-DR+","B","T")]
cds_blasts <- cds[,colData(cds)$aligned_partition %in% c("Early","Late","HLA-DR+")]
cds_blasts_nonormal <- cds_blasts[,colData(cds_blasts)$response != "Normal"]
cds_trimmed_nonormal <- cds_trimmed[,colData(cds_trimmed)$response != "Normal"]

# transform the data to be able to make a fractionated dotplot  
hgf_expr_data <-
  plot_genes_violin(cds_subset = cds_trimmed_nonormal[rowData(cds_trimmed_nonormal)$gene_short_name == "HGF", ], group_cells_by = "response")[["data"]] %>%
  mutate(prtd = paste0(aligned_partition,"_",response,"_",treatment_day1))
hgf_expr_data

binary_hgf_expr <- hgf_expr_data %>% 
  mutate(binary = ifelse(expression == 0, 0, 1)) %>% 
  group_by(aligned_partition,response,treatment_day1) %>% 
  summarise(binary_proportion = sum(binary)/n()) %>%
  mutate(prtd = paste0(aligned_partition,"_",response,"_",treatment_day1)) %>%
  ungroup() %>%
  select(prtd,binary_proportion)

hgf_expr_data_binary <- left_join(hgf_expr_data,binary_hgf_expr)

hgf_expr_data_plot <- hgf_expr_data_binary %>% 
  group_by(prtd) %>% 
  summarise(mean_exprs = mean(expression), 
            prop_exprs = mean(binary_proportion),
            sum_exprs = sum(expression),
            feature_label = unique(feature_label),
            aligned_partition = unique(aligned_partition),
            response = unique(response),
            treatment_day1 = unique(treatment_day1))
hgf_expr_data_plot


# use built in monocle tools to do the regression
colData(cds_trimmed_nonormal)$apd <- paste0(colData(cds_trimmed_nonormal)$aligned_partition,"_day_",colData(cds_trimmed_nonormal)$treatment_day1)
apd_list <- colData(cds_trimmed_nonormal)$apd %>% unique()
treatment_days <- c("0","1","2","3","42-44")
treatment_days <- factor(treatment_days, levels = c("0","1","2","3","42-44"))
aligned_partitions <- unique(colData(cds_trimmed_nonormal)$aligned_partition)

monocle_regression_func <- function(cds, gene_or_genes, stratification_variable, stratification_value, form, i = 1) {
  stratification_value <- stratification_value[i]
  cds <- cds[rowData(cds)$gene_short_name %in% gene_or_genes, colData(cds)[[stratification_variable]] == stratification_value]
  gene_fits <- fit_models(cds = cds, model_formula_str = form, expression_family = "negbinomial")
  fit_coefs <- coefficient_table(gene_fits)
  fit_coefs_return <- fit_coefs %>% filter(term != "(Intercept)") %>% mutate(stratification = stratification_value, formula = form) %>% select(stratification, formula, term, estimate, p_value, q_value)
  return(fit_coefs_return)
}


hgf_monocle_regression2 <- lapply(
  X = seq_along(treatment_days),
  FUN = monocle_regression_func,
  cds = cds_trimmed_nonormal[,colData(cds_trimmed_nonormal)$aligned_partition != "B"],
  gene_or_genes = "HGF",
  stratification_variable = "treatment_day1",
  stratification_value = sort(treatment_days),
  form = "~response+aligned_partition"
)
hgf_monocle_regression2 %>% discard(function(x) nrow(x) == 0) %>% bind_rows() %>% filter(term == "responseNR") %>% write_csv("data_out/hgf_dotplot_qvalues.csv")
hgf_monocle_regression2
# unique_apd <- colData(cds_trimmed_nonormal)$apd %>% unique()
# 
# hgf_monocle_regression4 <- lapply(
#   X = seq_along(unique_apd),
#   FUN = monocle_regression_func,
#   cds = cds_trimmed_nonormal,
#   gene_or_genes = "HGF",
#   stratification_variable = "apd",
#   stratification_value = sort(unique_apd),
#   form = "~response"
# )
# 
# hgf_monocle_regression4 %>% discard(function(x) nrow(x) == 0) %>% bind_rows() %>% View()
# 







# tack on pre-batch corrected expression data onto batch-corrected cds
##extract the counts data
hgf_data <- plot_cells(cds, genes = "HGF")[["plot_env"]][["data_df"]] %>% select(barcode_channel, HGF_raw_counts = value) %>% as_tibble()
#reorder by joining to the cds_aligned coldata 
hgf_data <- left_join(as_tibble(colData(cds_aligned)),hgf_data) %>% select(sanity_check = barcode_channel, HGF_raw_counts)
#tack onto the cds_aligned coldata as a metadatacolumn
colData(cds_aligned)[c("sanity_check","HGF_raw_counts")] <- hgf_data
#sanity check
sum(colData(cds_aligned)$barcode_channel != colData(cds_aligned)$sanity_check)
#remove sanity check
colData(cds_aligned)$sanity_check <- NULL
#normalize by size factor
colData(cds_aligned)$normalized_hgf <- colData(cds_aligned)$HGF_raw_counts/colData(cds_aligned)$Size_Factor


#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)
