source("00_bwb20200729.R")

# load the freemuxlet cell metadata
fml_data <- list.files("~/network/X/Labs/Blaser/Brad/collaborations/vicki/freemuxlet_data", full.name = TRUE)
fml_data_list <- lapply(X = fml_data, FUN=read_delim, delim = "\t",skip = 1,col_names = c("barcode","channel","nCount_RNA","nFeature_RNA","FML.DROPLET.TYPE","FNL.BEST.GUESS","doubletCalls","FML.BEST.GUESS.FORMATTED"))
fml_data_tbl <- bind_rows(fml_data_list) %>% mutate(barcode_channel = paste0(barcode,"_",channel))

# load up the sample/color information from the latin square
latin_square <- read_csv("~/network/X/Labs/Blaser/Brad/collaborations/vicki/ucsf_data_20200723/shared/shared/bwb_latin_square.csv") %>%
  mutate(sample_pool = paste0(ori.sample.name,"_",pool))

# load up ravi's cluster assignments
rp_clusters <- read_delim("~/network/X/Labs/Blaser/Brad/collaborations/vicki/ucsf_data_20200723/shared/shared/fmlclust_to_sample.txt", delim = "\t") %>%
  mutate(sample_pool = paste0(ori.sample.name,"_",str_replace(fmlclust,":CLUST[:digit:]*","")))

sample_ids <- left_join(rp_clusters,latin_square %>% select(-ori.sample.name), by = "sample_pool")  %>% 
  mutate(treatment_day1 = recode(treatment_day, "42" = "42-44","43" = "42-44","44" = "42-44"))  %>%
  na.omit()

fml_metadata_meh <- left_join(fml_data_tbl,sample_ids, by = c("FML.BEST.GUESS.FORMATTED" = "fmlclust"))

fml_metadata_full <- left_join(as_tibble(colData(cds0)), fml_metadata_meh %>% select(-barcode,-channel), by = "barcode_channel") %>% 
  select(-sample) #need to do this because a small number of droplets didn't get a freemuxlet score and the data frame needs to be the right dimensions before adding back into the cds

# print out a summary descriptive table of all cells after freemuxlet only
if (!dir.exists("freemuxlet_data_out")) {
  dir.create("freemuxlet_data_out")
}

fml_metadata_summary1 <- fml_metadata_full %>% 
  group_by(doubletCalls,FNL.BEST.GUESS) %>% 
  summarise(n = n()) %>%
  mutate(doubletCalls = replace_na(doubletCalls,"no freemuxlet data or ambiguous")) %>%
  mutate(FNL.BEST.GUESS = replace_na(FNL.BEST.GUESS, "no freemuxlet data or ambiguous")) %>%
  ungroup() %>% 
  write_csv("freemuxlet_data_out/fml_metadata_summary1.csv")

fml_metadata_summary2 <- fml_metadata_summary1 %>%
  group_by(doubletCalls) %>%
  summarise(total = sum(n)) %>%
  mutate(user_friendly = recode(doubletCalls, 
                                "DoubletFinder.DBL" = "Intra-patient doublets",
                                "Freemuxlet.DBL" = "Inter-patient doublets", 
                                "Singlet" = "Singlets")) %>%
  write_csv("freemuxlet_data_out/fml_metadata_summary2.csv")


# add fml data into cds object
cds1 <- cds0

colData(cds1)[,6:18] <- fml_metadata_full[,4:16]
# sanity check
sum(colData(cds1)$barcode_channel != colData(cds1)$barcode_channel.1) #looks good

# remove extra columns not needed at this point
colData(cds1)$barcode_channel.1 <- NULL
colData(cds1)$barcode.1 <- NULL
colData(cds1)$channel.1 <- NULL
colData(cds1)$Size_Factor.1 <- NULL

# filter the main cds object to include only singlets
cds2 <- cds1[,!is.na(colData(cds1)$doubletCalls)]#filters out cells that didn't get freemuxlet data calls
cds2.1  <- cds2[,colData(cds2)$doubletCalls == "Singlet"]#selects only freemuxlet singlets
cds2.2 <- cds2.1[,!is.na(colData(cds2.1)$pool)]#filters out ambiguous clusters from cyan 1 and cyan 2


#sanity check
colData(cds2.2)$doubletCalls %>% unique()# looks good
colData(cds2.2)$channel %>% unique()# looks good
colData(cds2.2)$pool %>% unique()# looks good
sum(colData(cds2.2)$channel != colData(cds2.2)$pool)# looks good

#remove the extra column
colData(cds2.2)$pool <- NULL

# make a descriptive table of the cells remaining after freemuxlet and doublet finder
singlet_only_summary <- as_tibble(colData(cds2.2)) %>% 
  group_by(VWsample) %>% 
  summarise(n = n()) %>% 
  mutate(running_total = cumsum(n)) %>% 
  write_csv("freemuxlet_data_out/singlet_only_summary.csv")

#save.image.pigz("ficlatuzumab_freemuxlet.RData")