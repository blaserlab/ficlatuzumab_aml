source("00_bwb20200729.R")

os_data <- read_csv("survival_data/os.csv")
pfs_data <- read_csv("survival_data/pfs.csv")
hsct_os_data <- read_csv("survival_data/hsct_os.csv")

os_fit <- survfit(Surv(Days,Status) ~ Response, data = os_data)
pfs_fit <- survfit(Surv(Days,Status) ~ Response, data = pfs_data)
hsct_os_fit <- survfit(Surv(Days,Status) ~ Response, data = hsct_os_data)

#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 30)
