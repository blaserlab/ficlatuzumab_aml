source("00_bwb20200729.R")

dir.create("combined_goterm_output")

harvest_goterm_data <- function(revigo, gorilla, outfile) {
  revigo_tbl <-
    read_csv(paste0("revigo_data/", revigo)) %>% 
    mutate(module = paste0("module_", str_extract(revigo, "[:digit:]+")))
  gorilla_tbl <-
    read_delim(paste0("gorilla_data/", gorilla),delim = "\t") %>% 
    mutate(module = paste0("module_", str_extract(gorilla, "[:digit:]+")))
  joined <- full_join(gorilla_tbl, revigo_tbl, by = c("GO Term" = "term_ID")) %>% 
    select(Module = module.x, 
           GO_term = `GO Term`, 
           Description, 
           gorilla_pval = `P-value`, 
           gorilla_qval = `FDR q-value`,
           revigo_log10p = `log10 p-value`,
           uniqueness,
           dispensability,
           representative,
           eliminated,
           plot_X,
           plot_Y,
           Genes) %>% 
    write_csv(outfile)
  
}

pmap(.l = list(
  revigo = list.files("revigo_data"),
  gorilla = list.files("gorilla_data"),
  outfile = paste0("combined_goterm_output/",str_replace(list.files("revigo_data"),"^[:alpha:]*_",""))
),
.f = harvest_goterm_data)


#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)