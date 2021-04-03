source("00_bwb20200729.R")

#system("gunzip -k ~/network/X/Labs/Blaser/ngs_archive/aml/genotypes/imputed.exons.filtered.vcf.gz")

#load the vcf file we want to show

vcf <- read_delim(file = "~/network/X/Labs/Blaser/ngs_archive/aml/genotypes/imputed.exons.filtered.vcf",delim = "\t",skip = 11, col_types = cols(ID = col_character()))

vcf_long <- pivot_longer(data = vcf, cols = matches("._."),names_to = "patient", values_to = "genotype")
vcf_long$genotype <- str_replace(vcf_long$genotype,":.*","")
vcf_long$fill_color <- dplyr::recode(vcf_long$genotype, "0|0" = "transparent","0|1" = "black","1|0" = "black","1|1" = "black")
vcf_long$patient %>% unique()
vcf_long$patient_pretty <- dplyr::recode(vcf_long$patient,
                                         "01_161690" = "E01",
                                         "02_161691" = "E02",
                                         "03_161699" = "E03",
                                         "04_161692" = "E04",
                                         "05_161693" = "E05",
                                         "06_161694" = "E06",
                                         "07_161695" = "E07",
                                         "09_161700" = "E09",
                                         "10_161701" = "E10",
                                         "11_161696" = "E11",
                                         "12_161697" = "E12",
                                         "13_161698" = "E13")

vcf_long$patient_pretty <- factor(vcf_long$patient_pretty, levels = c("E01","E02","E03","E04","E05","E06","E07","E09","E10","E11","E12","E13"))



#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 30)