source("00_packages_functions.R")

#make a new column with just the treatment day
colData(cds)$treatment_day <-
  colData(cds)$response_day_babis %>% 
  str_replace("[:alpha:]{1,}_", "")

#create the outupot directory
if (!dir.exists("plots_out/hgf_induced_targets")){
  dir.create(("plots_out/hgf_induced_targets"))
}

mclapply(
  X = seq_along(hgf_induced),
  FUN = plot_cells_alt_2,
  cds = cds[, colData(cds)$treatment_day != "day_NA"],
  gene_or_genes = hgf_induced,
  col_var = "treatment_day",
  row_var = "response_babis",
  title = hgf_induced,
  h = 3.75,
  w = 7,
  cell_size = 0.25,
  alpha = 0.5,
  outfile = paste0("plots_out/hgf_induced_targets/",hgf_induced,".pdf"),
  mc.preschedule = T,
  mc.cores = 39
)


chunk_func <- function(thing,chunksize){
  numchunks <- ceiling(length(thing)/chunksize)
  chunklist <- vector(mode = "list", length = numchunks)
  for (i in 1:numchunks) {
    chunklist[[i]] <- thing[1:chunksize]
    thing <- thing[chunksize+1:length(thing)]
  }
  return(chunklist)
}
hgf_repressed_list <- chunk_func(thing = hgf_repressed, chunksize = 10)

if (!dir.exists("plots_out/hgf_repressed_targets")){
  dir.create("plots_out/hgf_repressed_targets")
}

for (j in 1:length(hgf_repressed_list)) {
  hgf_repressed_chunked <- hgf_repressed_list[[j]]
  mclapply(
    X = seq_along(hgf_repressed_chunked),
    FUN = plot_cells_alt_2,
    cds = cds[, colData(cds)$treatment_day != "day_NA"],
    gene_or_genes = hgf_repressed_chunked,
    col_var = "treatment_day",
    row_var = "response_babis",
    title = hgf_repressed_chunked,
    h = 3.75,
    w = 7,
    cell_size = 0.25,
    alpha = 0.5,
    outfile = paste0(
      "plots_out/hgf_repressed_targets/",
      hgf_repressed_chunked,
      ".pdf"
    ),
    mc.preschedule = T,
    mc.cores = 10
  )
  
}

vicki_colon_genes <- c("ABCC4","API5","CAP1","CAPRIN1","CCNH","CCT8","CCZ1","CLCN3","CNOT8","CSE1L","HAT1","IFITM1","ITGAV","KIF2A","PRPS2","RAB11A","RNASE4","SEPTIN2","SFPQ","TM9SF3")

vicki_colon_list <- chunk_func(thing = vicki_colon_genes, chunksize = 1)

if (!dir.exists("plots_out/vicki_targets")){
  dir.create("plots_out/vicki_targets")
}

lapply(
  X = seq_along(vicki_colon_genes),
  FUN = possibly(plot_cells_alt_2, otherwise = NULL, quiet = FALSE),
  cds = cds[, colData(cds)$treatment_day != "day_NA"],
  gene_or_genes = vicki_colon_genes,
  col_var = "treatment_day",
  row_var = "response_babis",
  title = vicki_colon_genes,
  h = 3.75,
  w = 7,
  cell_size = 0.25,
  alpha = 0.5,
  outfile = paste0("plots_out/vicki_targets/",
                   vicki_colon_genes,
                   ".pdf")
)
  



save.image.pigz("ficlatuzumab_aml.RData",n.cores = 30)
