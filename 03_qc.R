source("00_20200729.R")

#get mito genes for qc
is.mito <- as_tibble(rowData(cds2.2)) %>% filter(str_detect(gene_short_name,"^MT-")==TRUE) %>% pull(id) %>% as.character()

# run basic qc stats on the whole dataset
cds3 <- scater::calculateQCMetrics(object = cds2.2, feature_controls = list(MITO=is.mito))

# get a vector of batches (channels)
channels <- as_tibble(colData(cds3)) %>% pull(channel) %>% unique()

# extract the cds as a tibble to make it easier to work with
cds_qc_tbl <- as_tibble(colData(cds3))

# make a qc function
qc_func <- function(cds_tbl, batch, i) {
  return(cds_tbl %>% 
           filter(channel == batch[i]) %>% 
           mutate(qc.detected = isOutlier(total_features_by_counts,log = TRUE,type = "lower"), 
                  qc.mito = isOutlier(pct_counts_MITO, type = "higher")))
  
}

qc_attr_feature_func <- function(cds_tbl, batch, i) {
  cds_tbl_filtered <- cds_tbl %>% 
    filter(channel == batch[i])
  qc_attr_data <- isOutlier(cds_tbl_filtered$total_features_by_counts, log = TRUE,type = "lower")
return(attr(qc_attr_data, "thresholds"))
  
}

qc_attr_mito_func <- function(cds_tbl, batch, i) {
  cds_tbl_filtered <- cds_tbl %>% 
    filter(channel == batch[i])
  qc_attr_data <- isOutlier(cds_tbl_filtered$pct_counts_MITO, type = "higher")
return(attr(qc_attr_data, "thresholds"))
  
}

# run the function to calculate identify outliers by channel
if (!dir.exists("qc_data_out")){
  dir.create("qc_data_out")
}

qc_list <- lapply(X = seq_along(channels), FUN = qc_func, cds_tbl = cds_qc_tbl, batch = channels)

qc_attr_mito <- lapply(X = seq_along(channels), FUN = qc_attr_mito_func, cds_tbl = cds_qc_tbl, batch = channels)
names(qc_attr_mito) <- channels
qc_thresholds_mito <- bind_rows(qc_attr_mito,.id = "channel") %>% write_csv("qc_data_out/qc_thresholds_mito.csv")

qc_attr_features <- lapply(X = seq_along(channels), FUN = qc_attr_feature_func, cds_tbl = cds_qc_tbl, batch = channels)
names(qc_attr_features) <- channels
qc_thresholds_features <- bind_rows(qc_attr_features,.id = "channel") %>% write_csv("qc_data_out/qc_thresholds_features.csv")

# make a table of just the barcode_channel and qc calls
cds_qc_tbl_thin <- bind_rows(qc_list) %>% 
  select(channel, barcode_channel,qc.detected,qc.mito) %>% 
  mutate(qc.any = qc.detected | qc.mito)

#summarize and print out qc stats
qc_stats <- cds_qc_tbl_thin %>% 
  group_by(channel) %>%
  summarise(qc.detected.sum = sum(qc.detected), qc.mito.sum = sum(qc.mito), qc.any.sum = sum(qc.any)) %>% 
  write_csv("qc_data_out/qc_stats.csv")

# bind the qc calls back onto the main cds
cds4 <- cds3
colData(cds4)[53:56] <- cds_qc_tbl_thin[2:5]
#sanity check
sum(colData(cds4)$barcode_channel != colData(cds4)$barcode_channel.1) #looks good.  I guess all the cells stayed in order
#sanity check sanity check
all.equal(colData(cds4)$barcode_channel,cds_qc_tbl_thin$barcode_channel) # still looks good so I guess this didn't need to be a join

#keep only the columns we really want at this point
cds5 <- cds4
colData(cds5) <- colData(cds5)[,c(1:3,5:12,14:16,48,54:56)]

# filter out cells that were flagged with qc.any = TRUE
cds6 <- cds5[,colData(cds5)$qc.any == FALSE]

#make a summary table of all cell counts to make sure we aren't losing any
cds_qc_list <- list(cds0,cds1,cds2,cds2.1,cds2.2,cds3,cds4,cds5,cds6)
cds_vec <- c("cds0","cds1","cds2","cds2.1","cds2.2","cds3","cds4","cds5","cds6")
description_vec <- c("made from combine_cds",
                    "added metadata columns",
                    "filtered out cells without freemuxlet data",
                    "filtered for singlets in doubletCalls",
                    "filtered out ambiguous cells from cyan 1 and cyan 2",
                    "calculated QC metrics",
                    "added metadata columns",
                    "removed metadata columns",
                    "removed low quality cells") 
nrow_vec <- sapply(X = cds_qc_list, FUN = function(x){nrow(colData(x))})
cds_cell_numbers <- tibble(cds_vec,description_vec,nrow_vec) %>% 
  select(cds = cds_vec, cell_number = nrow_vec) %>%
  mutate(user_friendly = recode(cds,
                                "cds0" = "starting number of cells",
                                "cds1" = "starting number of cells",
                                "cds2" = "after removing cells with no freemuxlet data",
                                "cds2.1" = "after removing inter- and intra-patient doublets",
                                "cds2.2" = "after removing freemuxlet ambiguous cells",
                                "cds3" = "after removing freemuxlet ambiguous cells",
                                "cds4" = "after removing freemuxlet ambiguous cells",
                                "cds5" = "after removing freemuxlet ambiguous cells",
                                "cds6" = "final dataset after removing low qualtiy cells")) %>%
  group_by(cell_number, user_friendly) %>%
  summarise(cds = dplyr::first(cds)) %>%
  arrange(desc(cell_number)) %>%
  ungroup() %>%
  mutate(difference = lag(cell_number,n = 1, default = cell_number[1])-cell_number) %>%
  mutate(percent_removed = difference/cell_number[1]*100)
  write_csv("qc_data_out/cds_cell_numbers.csv")

