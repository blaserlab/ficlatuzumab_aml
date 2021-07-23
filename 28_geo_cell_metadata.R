source("00_bwb20200729.R")

left_join(fml_metadata_full, colData(cds) %>% as_tibble()) %>%
  select(channel, barcode, patient = VWsample, treatment_day, response) %>%
  write_delim("local_data/ficlatuzumab_aml_geo/cell_metadata.tsv",delim = "\t")
