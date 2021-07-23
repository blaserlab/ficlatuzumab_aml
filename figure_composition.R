source("00_bwb20200729.R")

#fig3####-------------------------------------------------------------------
figure_3 <- plot_grid(
  plot_grid(
    var_plot_patient, cds_cp_plot_aligned, blank_plot, align = "h", axis = "b",nrow = 1,labels = c("A","B",""), rel_widths = c(10,8.5,2)
  ),
  plot_grid(gene_heatmap, labels = "C"),
  plot_grid(aml_marker_plot, hgf_dot_plot,align = "h", axis = "b", ncol = 2, rel_widths = c(1.1,1),labels = c("D","E")),
  nrow = 3,
  align = "v", 
  axis = "l"
)

save_plot(
  figure_3,
  filename = "figs_out/figure_3.pdf",
  base_width = 7.5,
  base_height = 8
)

save_plot(
  figure_3,
  filename = "figs_out/figure_3.png",
  base_width = 7.5,
  base_height = 8,
  dpi = 300
)
# figure 4####---------------------------------------------------------------------------

fig4bottom <- plot_grid(
  mod3_agg_violins, revigo_mod3,
  mod18_agg_violins, revigo_mod18, 
  mod1_agg_violins, revigo_mod1, 
  mod7_agg_violins, revigo_mod7, 
  nrow = 4,
  align = "hv", axis = "lb", labels = c("B","C","D","E","F","G","H","I"), hjust = 0
)

fig4top <- plot_grid(
  gene_module_rd, blank_plot, rel_widths = c(10,2),nrow = 1
)

figure_4 <- plot_grid(
  fig4top,
  fig4bottom,
  nrow = 2,
  rel_heights = c(0.8,2.4),
  labels = c("A",""), hjust = 0
)

save_plot(
  figure_4,
  filename = "figs_out/figure_4.pdf",
  base_width = 6.5,
  base_height = 10
)

save_plot(
  figure_4,
  filename = "figs_out/figure_4.png",
  base_width = 6.5,
  base_height = 10,
  dpi = 300
)

#supplemental_fig3####---------------------------------------------------------

figure_S3_top <- plot_grid(
  patient_plot_unaligned,
  aligned_partition_plot,
  align = "h",
  axis = "b",
  ncol = 2,
  rel_widths = c(10,8.5),
  labels = c("A","B")
)

figure_S3_middle <- plot_grid(
  seurat_plot,
  blank_plot,
  align = "h",
  axis = "b",
  ncol= 2,
  rel_widths = c(10,6),
  labels = c("C","")
)

figure_S3 <- plot_grid(
  figure_S3_top,
  figure_S3_middle,
  nrow = 2,
  rel_heights = c(10,13),
  align = "v",
  axis = "l"

)

save_plot(
  figure_S3,
  filename = "figs_out/figure_S3.pdf",
  base_width = 7.5,
  base_height = 7
)

save_plot(
  figure_S3,
  filename = "figs_out/figure_S3.png",
  base_width = 7.5,
  base_height = 7,
  dpi = 300
)

# supplemental figure S4####--------------------------------------------------
save_plot(
  cytof_scrnaseq_supervised,
  filename = "figs_out/figure_S4_supervised.pdf",
  base_width = 7.5,
  base_height = 6.5
)

save_plot(
  cytof_scrnaseq_supervised,
  filename = "figs_out/figure_S4_supervised.png",
  base_width = 7.5,
  base_height = 6.5,
  dpi = 300
)

# supplemental figure 5####---------------------------------------------------
S5_plots <- cowplot::align_plots(
  hgf_umap_plot,
  csf1r_umap_plot,
  csf1r_dot_plot,
  align = "v",
  axis = "l"
)

S5_bottom <- plot_grid(
  S5_plots[[3]],
  hgf_induced_dotplot,
  align = "h",
  axis = "b",
  labels = c("","D"),
  label_x = c(0,-0.05)
)


figure_S5 <- plot_grid(
  # hgf_umap_plot,
  S5_plots[[1]],
  S5_plots[[2]],
  S5_bottom,
  nrow = 3,
  labels = "AUTO"
)




save_plot(
  figure_S5,
  filename = "figs_out/figure_S5.pdf",
  base_width = 7.5,
  base_height = 6.5
)

save_plot(
  figure_S5,
  filename = "figs_out/figure_S5.png",
  base_width = 7.5,
  base_height = 8.5,
  dpi = 300
)


# cytof scrna celltype proportions #####--------------------------------------------------------------
save_plot(cytof_scrna_cellcount_plot, filename = "figs_out/cytof_scrna_cellcounts.pdf", base_height = 3.0, base_width = 4.125)


#save.image.pigz("ficlatuzumab_memory_managed.RData")
