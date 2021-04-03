source("00_bwb20200729.R")
library(Seurat)
library(SeuratDisk)
library(patchwork)

# load the reference dataset

seurat_reference <- LoadH5Seurat("~/workspace_pipelines/OSUBlaserlab_analyses/brad/scrnaseq_analysis/preprocessing/references/pbmc_multimodal.h5seurat")

DimPlot(object = seurat_reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# extract the count data from the main cds (unaligned) and rename with gene names
# the following block may be incorrect and might not run.  The correct script got lost when R crashed, but we have the data objects stored separately below.
# basically the idea here is to extract the counts table, replace ensdarg numbers with gene names and then make a seurat object out of this

exprs_renamed_long <- broom::tidy(exprs(cds)) %>%
  as_tibble() %>%
  left_join(.,rowData(cds) %>% as_tibble(), by = c("row" = "id")) %>%
  select(gene_short_name, column, value) %>%
  select(row = gene_short_name, column, value)
exprs_renamed_df <- as.data.frame(exprs_renamed_long)
exprs_renamed_rows <- exprs_renamed_df %>% pull(row) %>% unique() %>% data.frame()
exprs_renamed_cols <- exprs_renamed_df %>% pull(column) %>% unique() %>% data.frame()
exprs_renamed_df$rowIdx <- match(exprs_renamed_df$row, exprs_renamed_rows$row)
exprs_renamed_df$colIdx <- match(exprs_renamed_df$column, exprs_renamed_cols$column)



exprs_spars <- sparseMatrix(
  i =exprs_renamed_df$rowIdx,
  j = exprs_renamed_df$coldx,
  x = exprs_renamed_df$value,
  dimnames = list(exprs_renamed_df$row,exprs_renamed_df$column)
) 

#seurat_cds <- CreateSeuratObject(monocle3::exprs(cds),)
#seurat_cds <- SCTransform(seurat_cds, verbose = T)

anchors <- FindTransferAnchors(
  reference = seurat_reference,
  query = seurat_cds,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50, 
  features = intersect(rownames(x = seurat_reference), VariableFeatures(object = seurat_cds)),
  reference.assay = "SCT",
  query.assay = "SCT",
  verbose = T
)


seurat_cds <- TransferData(
  anchorset = anchors,
  reference = seurat_reference,
  query = seurat_cds,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  )
)

seurat_cds <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = seurat_reference,
  query = seurat_cds,
  new.reduction.name = "ref.spca"
)

seurat_cds <- ProjectUMAP(
  query = seurat_cds,
  query.reduction = "ref.spca",
  reference = seurat_reference,
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

seurat1 <- DimPlot(seurat_cds, 
             reduction = "ref.umap", 
             group.by = "predicted.celltype.l1", 
             label = TRUE, 
             label.size = 3, 
             repel = TRUE) + 
  NoLegend()

seurat1

seurat2 <- 
  DimPlot(seurat_cds, 
          reduction = "ref.umap", 
          group.by = "predicted.celltype.l2", 
          label = TRUE, 
          label.size = 3 ,
          repel = TRUE) + 
  NoLegend()

seurat_metadata <- seurat_cds@meta.data

save.image.pigz(file = "ficlatuzumab_data_with_seurat.RData", n.cores = 39)
save.pigz(
    seurat_reference,
    seurat_cds,
    seurat1,
    seurat2,
    anchors,
    exprs_renamed_cols,
    exprs_renamed_df,
    exprs_renamed_long,
    exprs_renamed_rows,
    exprs_spars,
  file = "seurat_objects_only.RData",
  n.cores = 39
)

rm(
  seurat_reference,
  seurat_cds,
  anchors,
  exprs_renamed_cols,
  exprs_renamed_df,
  exprs_renamed_long,
  exprs_renamed_rows,
  exprs_spars
)

seurat_metadata <- as_tibble(seurat_metadata,rownames = "colData_rownames")

seurat_colData <-
  left_join(colData(cds) %>% as_tibble(rownames = "colData_rownames"),
            seurat_metadata,
            by = "colData_rownames")

seurat_colData

colData(cds)$predicted.celltype.l2.score <- seurat_colData$predicted.celltype.l2.score
colData(cds)$predicted.celltype.l2 <- seurat_colData$predicted.celltype.l2




plot_cells(cds, color_cells_by = "predicted.celltype.l2")

save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)
