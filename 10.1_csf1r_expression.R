source("00_bwb20200729.R")

# transform the data to be able to make a fractionated dotplot  
csf1r_expr_data <-
  plot_genes_violin(cds_subset = cds_trimmed_nonormal[rowData(cds_trimmed_nonormal)$gene_short_name == "CSF1R", ], group_cells_by = "response")[["data"]] %>%
  mutate(prtd = paste0(aligned_partition,"_",response,"_",treatment_day1))
csf1r_expr_data

binary_csf1r_expr <- csf1r_expr_data %>% 
  mutate(binary = ifelse(expression == 0, 0, 1)) %>% 
  group_by(aligned_partition,response,treatment_day1) %>% 
  summarise(binary_proportion = sum(binary)/n()) %>%
  mutate(prtd = paste0(aligned_partition,"_",response,"_",treatment_day1)) %>%
  ungroup() %>%
  select(prtd,binary_proportion)

csf1r_expr_data_binary <- left_join(csf1r_expr_data,binary_csf1r_expr)

csf1r_expr_data_plot <- csf1r_expr_data_binary %>% 
  group_by(prtd) %>% 
  summarise(mean_exprs = mean(expression), 
            prop_exprs = mean(binary_proportion),
            sum_exprs = sum(expression),
            feature_label = unique(feature_label),
            aligned_partition = unique(aligned_partition),
            response = unique(response),
            treatment_day1 = unique(treatment_day1))
csf1r_expr_data_plot


# use built in monocle tools to do the regression
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

# hgf_monocle_regression

csf1r_monocle_regression2 <- lapply(
  X = seq_along(treatment_days),
  FUN = monocle_regression_func,
  cds = cds_trimmed_nonormal[,colData(cds_trimmed_nonormal)$aligned_partition != "B"],
  gene_or_genes = "CSF1R",
  stratification_variable = "treatment_day1",
  stratification_value = sort(treatment_days),
  form = "~response+aligned_partition"
)
csf1r_monocle_regression2 %>% discard(function(x) nrow(x) == 0) %>% bind_rows() %>% filter(term == "responseNR") %>% write_csv("data_out/csf1r_dotplot_qvalues.csv")
csf1r_monocle_regression2

# tack on pre-batch corrected expression data onto batch-corrected cds
##extract the counts data
csf1r_data <- plot_cells(cds, genes = "CSF1R")[["plot_env"]][["data_df"]] %>% select(barcode_channel, CSF1R_raw_counts = value) %>% as_tibble()
#reorder by joining to the cds_aligned coldata 
csf1r_data <- left_join(as_tibble(colData(cds_aligned)),csf1r_data) %>% select(sanity_check = barcode_channel, CSF1R_raw_counts)
#tack onto the cds_aligned coldata as a metadatacolumn
colData(cds_aligned)[c("sanity_check","CSF1R_raw_counts")] <- csf1r_data
#sanity check
sum(colData(cds_aligned)$barcode_channel != colData(cds_aligned)$sanity_check)
#remove sanity check
colData(cds_aligned)$sanity_check <- NULL
#normalize by size factor
colData(cds_aligned)$normalized_csf1r <- colData(cds_aligned)$CSF1R_raw_counts/colData(cds_aligned)$Size_Factor

save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)
