source("00_bwb20200729.R")

# Normalize and pre-process the data
cds_preprocess <- preprocess_cds(cds_preprocess, num_dim = 100)#generates PC's
plot_pc_variance_explained(cds_preprocess)

#align using batchelor and VWsample variable
cds_aligned <- align_cds(cds = cds_preprocess, alignment_group = "VWsample")

#align using batchelor and channel variable
cds_aligned_channel <- align_cds(cds = cds_preprocess, alignment_group = "channel")

## Reduce dimensionality and visualize cells
cds <- reduce_dimension(cds_preprocess, cores = 39)
cds_aligned <- reduce_dimension(cds_aligned, cores = 39)
cds_aligned_channel <- reduce_dimension(cds_aligned_channel, cores = 39)

# previz cells

previz_noalign <- custom_variable_plot(cds, var = "VWsample",foreground_alpha = 0.2)
previz_align_patient <- custom_variable_plot(cds_aligned, var = "VWsample", foreground_alpha = 0.2)
previz_align_channel <- custom_variable_plot(cds_aligned_channel, var = "VWsample", foreground_alpha = 0.2)

save.image.pigz("ficlatuzumab_memory_managed.RData", n.cores = 37)
