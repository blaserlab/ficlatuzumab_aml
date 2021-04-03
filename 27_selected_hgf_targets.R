source("00_bwb20200729.R")

hgf_induced <- c("HGF","ACOT11","TOMM34","TRMT2A","ZNF136","RAB38","PIBF1","KIAA1524","HEATR5A","DENND1B")

hgf_induced_dotplot <- custom_gene_dotplot(cds = cds[, colData(cds)$aligned_partition %in% c("Early", "Late", "HLA-DR+")],
                    markers = hgf_induced,
                    group_cells_by = "response",max.size = 4) +
  scale_color_gradient(low = "grey90", high =  "#DC0000") +
  guides(size = guide_legend(order = 2,title = "Proportion")) +
  labs(color = "Expression",x = NULL,y = NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
hgf_induced_dotplot


hgf_repressed <- c("HGF","ALDH1A1","APOL3","ARHGEF3","CASP1","CD3G","CD40","CORO7","DHRS3","ENGASE","GBP2","HLA-F","HMGB1","METTL7A","OS9","PDK3","PDK4","PDLIM5","PLAC8","RNF207","TMTC1","TNFRSF1b","TXNIP","WSB1","ZCCHC24","ZNF467")

hgf_repressed_dotplot <- custom_gene_dotplot(cds = cds[, colData(cds)$aligned_partition %in% c("Early", "Late", "HLA-DR+")],
                    markers = hgf_repressed,
                    group_cells_by = "response") +
  scale_color_gradient(low = "grey90", high =  "#3C5488") +
  guides(size = guide_legend(order = 2,title = "Prop.")) +
  labs(color = "Expr.",x = NULL,y = NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
hgf_repressed_dotplot


#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)