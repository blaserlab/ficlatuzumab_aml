read_csv("~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-22-17cyan1/outs/metrics_summary.csv")

pipestance_paths <- c(
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-22-17cyan1",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-22-17cyan2",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-29-17blue1",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-29-17blue2",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-30-17green1",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_8-30-17green2",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_9-5-17red1",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_9-5-17red2",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_9-5-17yellow1",
  "~/network/X/Labs/Blaser/ngs_archive/aml/VW_9-5-17yellow2"
)

harvest_metrics <- function(pipestance_path) {
  metrics <- read_csv(paste0(pipestance_path, "/outs/metrics_summary.csv")) %>%
    mutate(channel = str_extract(pipestance_path,"cyan1|cyan2|blue1|blue2|green1|green2|red1|red2|yellow1|yellow2"))
}

metrics_10X <- bind_rows(lapply(X = pipestance_paths,
                                FUN = harvest_metrics)) 

