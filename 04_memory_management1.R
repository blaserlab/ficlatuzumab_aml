source("00_bwb20200729.R")

save.image.pigz(file = "ficlatuzumab_freemuxlet_01_to_03.RData",n.cores = 37)

cds_preprocess <- cds6
rm(
  blue1,
  blue2,
  cds_list,
  cds_qc_list,
  cds_qc_tbl,
  cds_qc_tbl_thin,
  cds0,
  cds1,
  cds2,
  cds2.1,
  cds2.2,
  cds3,
  cds4,
  cds5,
  cds6,
  cyan1,
  cyan2,
  green1,
  green2,
  qc_list,
  red1,
  red2,
  yellow1,
  yellow2
)

save.image.pigz(file = "ficlatuzumab_memory_managed.RData", n.cores = 37)
