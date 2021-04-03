source("00_bwb20200729.R")

freemuxlet_matrix_data <-
  read_csv("freemuxlet_data_out/fml_metadata_summary1.csv") %>%
  filter(doubletCalls != "DoubletFinder.DBL") %>%
  mutate(best_guess1 = fct_relevel(str_extract(FNL.BEST.GUESS, "^[:digit:]+"), c("0","1","2","3","4","5","6","7","8","9","10","11"))) %>%
  mutate(best_guess2 = fct_relevel(str_extract(FNL.BEST.GUESS, "[:digit:]+$"), c("0","1","2","3","4","5","6","7","8","9","10","11")))



fml_metadata_summary1 %>%
  group_by(doubletCalls) %>%
  summarise(n = sum(n))
fml_metadata_full %>%
  group_by(doubletCalls) %>%
  summarise(n = n())

colData(cds_aligned)
