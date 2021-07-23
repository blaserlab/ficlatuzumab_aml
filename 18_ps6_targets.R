source("00_bwb20200729.R")

ps6_genes <-
  c(
    "NOP56",
    "NOP14",
    "GAR1",
    "RRP9",
    "RRP15",
    "RRP12",
    "PWP2",
    "DDX18",
    "EIF4EBP1",
    "EIF4B",
    "EEF2K",
    "POLDIP3"
  )

#start by extracting data for each of these genes from pseudobulk data already calculated
if(!dir.exists("data_out/s6_targets")){
  dir.create("data_out/s6_targets")
}
pseudobulk_res_by_day[[1]][[2]] %>% filter(gene_short_name %in% ps6_genes) %>% arrange(padj)# %>% write_csv("data_out/s6_targets/s6_pseudobulk_day0.csv")#positive l2fc indicates up in NR
pseudobulk_res_by_day[[2]][[2]] %>% filter(gene_short_name %in% ps6_genes) %>% arrange(padj)# %>% write_csv("data_out/s6_targets/s6_pseudobulk_day1.csv")#positive l2fc indicates up in NR
pseudobulk_res_by_day[[3]][[2]] %>% filter(gene_short_name %in% ps6_genes) %>% arrange(padj)# %>% write_csv("data_out/s6_targets/s6_pseudobulk_day2.csv")#positive l2fc indicates up in NR
pseudobulk_res_by_day[[4]][[2]] %>% filter(gene_short_name %in% ps6_genes) %>% arrange(padj)# %>% write_csv("data_out/s6_targets/s6_pseudobulk_day3.csv")#positive l2fc indicates up in NR
pseudobulk_res_by_day[[5]][[2]] %>% filter(gene_short_name %in% ps6_genes) %>% arrange(padj)# %>% write_csv("data_out/s6_targets/s6_pseudobulk_day42_44.csv")#positive l2fc indicates up in NR

# now run the builtin monocle functions on cds_blasts_nonormal.  This is unaligned data
cds_s6_targets <- cds_blasts_nonormal[rowData(cds_blasts_nonormal)$gene_short_name %in% ps6_genes,]
s6_gene_fits <- fit_models(cds_s6_targets, model_formula_str = "~response")
s6_fit_coefs <- coefficient_table(s6_gene_fits)
s6_fit_coefs %>% filter(term != "(Intercept)") %>% select(gene_short_name, term, q_value, estimate) %>% arrange((q_value)) %>% write_csv("data_out/s6_targets/s6_monocle_regression_all.csv")

# subfractionate by treatment day.  Positive indicates up in NR compared to CR
s6_gene_fits_d0 <- fit_models(cds_s6_targets[,colData(cds_s6_targets)$treatment_day1 == "0"], model_formula_str = "~response")
s6_fit_coefs_d0 <- coefficient_table(s6_gene_fits_d0)
s6_fit_coefs_d0 %>% filter(term != "(Intercept)") %>% select(gene_short_name, term, q_value, estimate) %>% arrange((q_value))# %>% write_csv("data_out/s6_targets/s6_monocle_regression_d0.csv")

s6_gene_fits_d1 <- fit_models(cds_s6_targets[,colData(cds_s6_targets)$treatment_day1 == "1"], model_formula_str = "~response")
s6_fit_coefs_d1 <- coefficient_table(s6_gene_fits_d1)
s6_fit_coefs_d1 %>% filter(term != "(Intercept)") %>% select(gene_short_name, term, q_value, estimate) %>% arrange((q_value))# %>% write_csv("data_out/s6_targets/s6_monocle_regression_d1.csv")

s6_gene_fits_d2 <- fit_models(cds_s6_targets[,colData(cds_s6_targets)$treatment_day1 == "2"], model_formula_str = "~response")
s6_fit_coefs_d2 <- coefficient_table(s6_gene_fits_d2)
s6_fit_coefs_d2 %>% filter(term != "(Intercept)") %>% select(gene_short_name, term, q_value, estimate) %>% arrange((q_value))# %>% write_csv("data_out/s6_targets/s6_monocle_regression_d2.csv")

s6_gene_fits_d3 <- fit_models(cds_s6_targets[,colData(cds_s6_targets)$treatment_day1 == "3"], model_formula_str = "~response")
s6_fit_coefs_d3 <- coefficient_table(s6_gene_fits_d3)
s6_fit_coefs_d3 %>% filter(term != "(Intercept)") %>% select(gene_short_name, term, q_value, estimate) %>% arrange((q_value))# %>% write_csv("data_out/s6_targets/s6_monocle_regression_d3.csv")

s6_gene_fits_d42_44 <- fit_models(cds_s6_targets[,colData(cds_s6_targets)$treatment_day1 == "42-44"], model_formula_str = "~response")
s6_fit_coefs_d42_44 <- coefficient_table(s6_gene_fits_d42_44)
s6_fit_coefs_d42_44 %>% filter(term != "(Intercept)") %>% select(gene_short_name, term, q_value, estimate) %>% arrange((q_value))# %>% write_csv("data_out/s6_targets/s6_monocle_regression_d42_44.csv")
