library(tidyverse)
library(cowplot)
library(monocle3)
library(Seurat)
library(fastSave)
library(RColorBrewer)
library(Matrix)
library(igraph)
library(synchrony)
library(ggplotify)
library(foreach)
library(survminer)
library(survival)
library(broom)
library(scater)
library(DescTools)
library(Matrix.utils)
library(magrittr)
library(textshape)
#library(DESeq2)
library(pheatmap)
library(MASS)
library(ggrepel)
library(grid)
library(ComplexHeatmap)
library(patchwork)
library("Hmisc")
theme_set(theme_cowplot(font_size = 11))
mutate <- dplyr::mutate
filter <- dplyr::filter
arrange <- dplyr::arrange
pull <- dplyr::pull
select <- dplyr::select

#custom operators
`%notin%` <- Negate(`%in%`)
empty_data <- data.frame(
  horiz = c(1,2),
  vert = c(3,4)
)
# blank plot to fill in plot_grid
blank_plot<-ggplot(data = empty_data, aes(x = horiz, y = vert))+
  geom_point(color = "transparent")+
  theme(axis.ticks = element_line(color = "transparent"),
        axis.title = element_text(color = "transparent"),
        axis.line = element_line(color = "transparent"),
        axis.text = element_text(color = "transparent"))




#functions
tbl_to_matrix <- function(data) {
  data <- data %>%
    as.data.frame()
  rownames(data) <- data[,1]
  return(as.matrix(data[,-1]))
}

my.aggregate.Matrix = function (x, groupings = NULL, form = NULL, fun = "sum", ...)
{
  if (!methods::is(x, "Matrix"))
    x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  if (fun == "count")
    x <- x != 0
  groupings2 <- data.frame(A=as.factor(groupings))
  if (is.null(form))
    form <- stats::as.formula("~0+.")
  form <- stats::as.formula(form)
  mapping <- Matrix.utils::dMcast(groupings2, form)
  colnames(mapping) <- substring(colnames(mapping), 2)
  result <- Matrix::t(mapping) %*% x
  if (fun == "mean")
    result <- result/as.numeric(table(groupings)[rownames(result)])
  attr(result, "crosswalk") <- grr::extract(groupings, match(rownames(result),
                                                             groupings2$A))
  return(result)
}



remake_gxp_vicki<-function(cds,gene, title,alpha,cellsize,h,w,legend_position,outdir){
  p<-plot_cells_alt(cds = cds, 
                 genes = gene, 
                 reduction_method = "UMAP",
                 group_cells_by = "partition",
                 show_trajectory_graph = FALSE,
                 label_cell_groups = FALSE,
                 label_branch_points = FALSE,
                 label_roots = FALSE,
                 label_leaves = FALSE,
                 cell_size = cellsize,alpha = alpha)+
    theme_cowplot(font_size = 8)+
    theme(legend.position = legend_position)+
    theme(strip.background = element_rect(fill = "transparent"))+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  save_plot(plot = p, filename = paste0(outdir,"/",gene,".pdf"),base_height = h, base_width = w)
}

# custom scrnaseq functions

plot_cells_alt <-
  function (cds,
            gene_or_genes,
            cell_size = 1,
            alpha = 1 ,
            ncol = NULL,
            plot_title = NULL,
            legend_position = "right") {
    data <- plot_cells(cds = cds, genes = gene_or_genes)[["plot_env"]][["data_df"]]
    data$gene_short_name <-
      factor(data$gene_short_name, levels = gene_or_genes)
    background_data <- data %>% filter(is.na(value))
    foreground_data <- data %>% filter(!is.na(value))
    p <- ggplot() +
      geom_point(
        data = background_data,
        aes(x = data_dim_1, y = data_dim_2),
        color = "grey80",
        shape = 1,
        size = cell_size,
        stroke = 0.25
      ) +
      geom_point(
        data = foreground_data,
        aes(
          x = data_dim_1,
          y = data_dim_2,
          color = value/Size_Factor
        ),
        shape = 16,
        size = cell_size,
        alpha = alpha
      ) +
      scale_color_viridis_c(end = 0.8) +
      labs(
        x = "UMAP 1",
        y = "UMAP 2",
        color = "Expr.",
        title = plot_title
      ) +
      theme(legend.position = legend_position) +
      theme(plot.title = element_text(hjust = 0.5))
    if (length(gene_or_genes)>1){
      p <- p +
        facet_wrap(facets = vars(feature_label)) +
        theme(strip.background = element_blank())
    }
    return(p)
  }

custom_variable_plot <- function(cds,
                                 i = NULL,
                                 var,
                                 value_to_highlight = NULL,
                                 foreground_alpha = 1,
                                 legend_pos = "right",
                                 cell_size = 0.5,
                                 legend_title = NULL,
                                 plot_title = NULL,
                                 palette = NULL,
                                 legend_ncol = 1,
                                 log10transform = FALSE) {
  if (!is.null(i)) {
    value_to_highlight <- value_to_highlight[[i]]
    plot_title <- plot_title[[i]]
    outfile <- paste0("plots_out/", outfile[[i]], ".pdf")
  }
  
  data <- plot_cells(cds)[["data"]]
  data_long <-
    data %>% pivot_longer(cols = matches(var), names_to = "var")
  plot <- ggplot()
  if (!is.null(value_to_highlight)) {
    data_background <-
      data_long %>% filter(value %notin% value_to_highlight)
    data_long <- data_long %>% filter(value %in% value_to_highlight)
    plot <- plot +
      geom_point(
        data = data_background,
        aes(x = data_dim_1, y = data_dim_2),
        stroke = 0.25,
        shape = 1,
        size = cell_size,
        color = "grey80"
      )
  }
  if (class(data_long$value)=="numeric"){
    data_background <- data_long %>% filter(is.na(value))
    data_long <- data_long %>% filter(!is.na(value))
    plot <- plot + 
      geom_point(
        data = data_background,
        aes(x = data_dim_1, y = data_dim_2),
        stroke = 0.25,
        shape = 1, 
        size = cell_size,
        color = "grey80"
      )
  }
  if (log10transform == TRUE) {
    data_long <- data_long %>% mutate(value = log10(value))
  }
  plot <- plot +
    geom_point(
      data = data_long,
      aes(
        x = data_dim_1,
        y = data_dim_2,
        fill = value,
        color = value
      ),
      stroke = 0.25,
      shape = 21,
      alpha = foreground_alpha,
      size = cell_size
    )
  if (class(data_long$value) == "numeric") {
    plot <- plot + 
      scale_fill_viridis_c(guide = F) + 
      scale_color_viridis_c()
  } else if (!is.null(palette)) {
    if (palette[1] == "viridis") {
      plot <- plot +
        scale_fill_viridis_d(begin = 0.1, end = 0.9) +
        scale_color_viridis_d(begin = 0.1,
                              end = 0.9,
                              guide = F)
    } else if (palette[1] == "rcolorbrewer") {
      colourCount = length(unique(data_long$value))
      getPalette = colorRampPalette(brewer.pal(12, "Paired"))
      plot <- plot +
        scale_color_manual(values = getPalette(colourCount), guide = F) +
        scale_fill_manual(values = getPalette(colourCount))
    }  else if (length(palette)>1) {
      plot <- plot +
        scale_color_manual(values = palette, guide = F) +
        scale_fill_manual(values = palette)
    }
    
  }
  else {
    plot <- plot +
      scale_color_discrete(guide = F) +
      scale_fill_discrete()
  }
  
  if (class(data_long$value) != "numeric") {
    plot <-
      plot + guides(fill = guide_legend(
        ncol = legend_ncol,
        override.aes = list(
          size = 2,
          alpha = 1,
          color = "transparent"
        )
      ))
  }
  plot <-
    plot + labs(
      x = "UMAP 1",
      y = "UMAP 2",
      title = plot_title,
      fill = legend_title
    ) + theme(plot.title = element_text(hjust = 0.5))
  plot <- plot + theme(legend.position = legend_pos)#+coord_fixed()
  return(plot)
  
}


custom_cp_plot <- function(cds,
                           var = NULL,
                           alpha = 1,
                           cp ,
                           overwrite_labels = T,
                           legend_pos = "none",
                           cell_size = 0.5,
                           legend_title = NULL,
                           plot_title = NULL,
                           group_label_size = 3,
                           value_to_highlight = NULL,
                           alt_color_var = NULL,
                           facet_by = NULL,
                           palette = NULL) {
    #extract the data from the input cds
  data <- plot_cells(cds)[["data"]]
  #convert to long format
  data_long <- pivot_longer(data = data,
                            cols = cp,
                            names_to = "cp")
  # generate text data frame for cluster/partition labels
  text_df <- data_long %>% group_by(value)
  median_coord_df <-
    data_long %>% group_by(value) %>% summarise(
      fraction_of_group = n(),
      text_x = median(data_dim_1),
      text_y = median(data_dim_2)
    )
  text_df <- left_join(text_df, median_coord_df) %>%
    mutate(label = value)
  text_df <-
    text_df %>% group_by(label) %>% summarise(text_x = dplyr::first(text_x),
                                              text_y = dplyr::first(text_y))
  # if highlighting a categorical variable, generate background data and keep data_long as foreground
  if (!is.null(value_to_highlight)) {
    background_data_long <- data_long %>% filter((!!sym(var)) != value_to_highlight)
    data_long <- data_long %>% filter((!!sym(var)) %in% value_to_highlight)
  }
  #initialize the plot
  plot <- ggplot()
  # lay down the background plot if using
  if(!is.null(value_to_highlight)){
    plot<-plot+
      geom_point(data = background_data_long,
                 aes(x = data_dim_1, y = data_dim_2),
                 stroke = 0.25,
                 shape = 1,
                 size = cell_size,
                 color = "grey80")
  }

  # make the main colored plot from data_long
  if (!is.null(alt_color_var)){
    plot <- plot +
      geom_point(
        data = data_long,
        aes(
          x = data_dim_1,
          y = data_dim_2,
          fill = (!!sym(alt_color_var)),
          color = (!!sym(alt_color_var))),
        stroke = 0.25,
        shape = 21,
        alpha = alpha,
        size = cell_size
        )
  } else {
    plot <- plot +
      geom_point(
        data = data_long,
        aes(
          x = data_dim_1,
          y = data_dim_2,
          fill = value,
          color = value),
        stroke = 0.25,
        shape = 21,
        alpha = alpha,
        size = cell_size)
  }

  plot<-plot+
    scale_color_discrete(guide = F) +
    scale_fill_discrete()+
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      title = plot_title,
      fill = legend_title) +
    theme(plot.title = element_text(hjust = 0.5))
  # overwrite labels if you want to
  if (overwrite_labels == T) {
    plot <- plot +
      theme(legend.position = legend_pos) +
      ggrepel::geom_text_repel(
        data = text_df,
        mapping = aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size,min.segment.length = 1,
        box.padding = 0.33,
        ) +
      guides(fill = guide_legend(override.aes = list(
        size = 2,
        alpha = 1,
        color = "transparent"
      )))
  } else {
    plot <- plot +
      theme(legend.position = legend_pos) +
      guides(fill = guide_legend(override.aes = list(
        size = 2,
        alpha = 1,
        color = "transparent"
      )))}
  #option to facet
  if (!is.null(facet_by)) {
    plot <- plot +
      facet_wrap(facets = vars(!!sym(facet_by)),) +
      theme(strip.background = element_blank())
  }
  if (!is.null(palette) && palette[1] =="rcolorbrewer") {
    colourCount = length(unique(data_long$value))
    getPalette = colorRampPalette(brewer.pal(12, "Paired"))
    plot<-plot+
      scale_color_manual(values = getPalette(colourCount), guide = F)+
      scale_fill_manual(values = getPalette(colourCount))
  }

  if (!is.null(palette) && length(palette)>1) {
    plot <- plot +
      scale_color_manual(values = palette, guide = F) +
      scale_fill_manual(values = palette)
  }

    return(plot)
}



custom_violin_plot <-
  function(cds,
           variable,
           genes_to_plot,
           outfile = NULL,
           pseudocount = 0,
           include_jitter = FALSE,
           plot_title = NULL,
           w,
           h,
           rows = 1,
           show_x_label = TRUE,
           legend_pos = "none",
           jitter_alpha = 0.2,
           facet_var = NULL,
           violin_scale = c("count", "area", "width"),
           log_values = c(TRUE, FALSE))
  {
    data_to_plot <-
      plot_genes_violin(cds_subset = cds[rowData(cds)$gene_short_name %in% genes_to_plot,], group_cells_by = variable)[["data"]]
    if (log_values == FALSE) {
      p1 <-
        ggplot(data = data_to_plot, aes(x = !!as.name(variable),
                                        y = expression)) #expression already normalized when data extracted by violin plot function
    } else {
      p1 <-
        ggplot(
          data = data_to_plot,
          aes(
            x = !!as.name(variable),
            y = log10(expression + pseudocount))) #expression already normalized when data extracted by violin plot function
    } 
    p1 <- p1 +
      geom_violin(scale = violin_scale)
    if (include_jitter == T) {
      p1 <-
        p1 + geom_jitter(
          shape = 16,
          size = 1,
          color = "black",
          alpha = jitter_alpha,
          width = 0.2
        )
    }
    if (!is.null(facet_var)) {
      p1 <- p1 +
        facet_wrap(facets = vars(!!sym(facet_var)), nrow = rows) +
        theme(strip.background = element_rect(fill = "transparent"))
    }
    if (show_x_label == F) {
      p1 <- p1 + theme(axis.text.x = element_blank())
    }
    if (!is.null(outfile)) {
      save_plot(
        p1,
        filename = outfile,
        base_width = w,
        base_height = h
      )
    }
    return(p1)
  }
           
plot_genes_in_pseudotime_alt<-function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL, 
                                        ncol = 1, panel_order = NULL, color_cells_by = "pseudotime", 
                                        trend_formula = "~ splines::ns(pseudotime, df=3)", label_by_short_name = TRUE, 
                                        vertical_jitter = NULL, horizontal_jitter = NULL) 
{
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  tryCatch({
    pseudotime(cds_subset)
  }, error = function(x) {
    stop(paste("No pseudotime calculated. Must call order_cells first."))
  })
  colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
  if (!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  if (!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)), 
                            msg = paste("When label_by_short_name = TRUE,", "rowData must have a column of gene", 
                                        "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                "partition") | color_cells_by %in% names(colData(cds_subset)), 
                          msg = paste("color_cells_by must be a column in the", 
                                      "colData table."))
  if (!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in% rowData(cds_subset)$gene_short_name))
    }
    else {
      assertthat::assert_that(all(panel_order %in% row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100, 
                          msg = paste("cds_subset has more than 100 genes -", "pass only the subset of the CDS to be", 
                                      "plotted."))
  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))
  assertthat::assert_that("pseudotime" %in% names(colData(cds_subset)), 
                          msg = paste("pseudotime must be a column in", "colData. Please run order_cells", 
                                      "before running", "plot_genes_in_pseudotime."))
  if (!is.null(min_expr)) {
    assertthat::assert_that(assertthat::is.number(min_expr))
  }
  assertthat::assert_that(assertthat::is.number(cell_size))
  assertthat::assert_that(!is.null(size_factors(cds_subset)))
  if (!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)), 
                            msg = paste("When label_by_short_name = TRUE,", "rowData must have a column of gene", 
                                        "names called gene_short_name."))
  }
  assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                "partition") | color_cells_by %in% names(colData(cds_subset)), 
                          msg = paste("color_cells_by must be a column in the", 
                                      "colData table."))
  if (!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in% rowData(cds_subset)$gene_short_name))
    }
    else {
      assertthat::assert_that(all(panel_order %in% row.names(rowData(cds_subset))))
    }
  }
  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100, 
                          msg = paste("cds_subset has more than 100 genes -", "pass only the subset of the CDS to be", 
                                      "plotted."))
  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[, is.finite(colData(cds_subset)$pseudotime)]
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  if (is.null(min_expr)) {
    min_expr <- 0
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", 
                     by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_colData, by.x = "Cell", 
                     by.y = "row.names")
  cds_exprs$adjusted_expression <- cds_exprs$expression
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  new_data <- data.frame(pseudotime = colData(cds_subset)$pseudotime)
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)
  model_expectation <- model_predictions(model_tbl, new_data = colData(cds_subset))
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), 
                             function(x) {
                               data.frame(expectation = model_expectation[x$f_id, 
                                                                          x$Cell])
                             })
  cds_exprs <- merge(cds_exprs, expectation)
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (!is.null(panel_order)) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
                                      levels = panel_order)
  }
  q <- ggplot(aes(pseudotime, expression), data = cds_exprs)
  if (!is.null(color_cells_by)) {
    q <- q + geom_point(aes_string(color = color_cells_by), 
                        size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                        vertical_jitter))
    if (class(colData(cds_subset)[, color_cells_by]) == "numeric") {
      q <- q + viridis::scale_color_viridis(option = "C")
    }
  }
  else {
    q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter, 
                                                                        vertical_jitter))
  }
  q <- q + geom_line(aes(x = pseudotime, y = expectation), 
                     data = cds_exprs)
  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
                                        ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }
  q <- q + ylab("Expression")
  q <- q + xlab("pseudotime")
  q <- q + monocle_theme_opts()
  q
}


custom_gene_dotplot <- function (cds, markers, group_cells_by = "cluster", reduction_method = "UMAP", 
                                 norm_method = c("log", "size_only"), lower_threshold = 0, 
                                 max.size = 10, ordering_type = c("cluster_row_col", "maximal_on_diag", 
                                                                  "none"), axis_order = c("group_marker", "marker_group"), 
                                 flip_percentage_mean = FALSE, pseudocount = 1, scale_max = 3, 
                                 scale_min = -3,gene_alias_vector = NULL) 
{
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  if (!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% c("cluster", 
                                                  "partition") | group_cells_by %in% names(colData(cds)), 
                            msg = paste("group_cells_by must be a column in", 
                                        "the colData table."))
  }
  norm_method = match.arg(norm_method)
  gene_ids = as.data.frame(fData(cds)) %>% tibble::rownames_to_column() %>% 
    dplyr::filter(rowname %in% markers | gene_short_name %in% 
                    markers) %>% dplyr::pull(rowname)
  if (length(gene_ids) < 1) 
    stop(paste("Please make sure markers are included in the gene_short_name\",\n               \"column of the fData!"))
  if (flip_percentage_mean == FALSE) {
    major_axis <- 1
    minor_axis <- 2
  }
  else if (flip_percentage_mean == TRUE) {
    major_axis <- 2
    minor_axis <- 1
  }
  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids, ]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c("Cell", "Gene", "Expression")
  exprs_mat$Gene <- as.character(exprs_mat$Gene)
  if (group_cells_by == "cluster") {
    cell_group <- tryCatch({
      clusters(cds, reduction_method = reduction_method)
    }, error = function(e) {
      NULL
    })
  }
  else if (group_cells_by == "partition") {
    cell_group <- tryCatch({
      partitions(cds, reduction_method = reduction_method)
    }, error = function(e) {
      NULL
    })
  }
  else {
    cell_group <- colData(cds)[, group_cells_by]
  }
  if (length(unique(cell_group)) < 2) {
    stop(paste("Only one type in group_cells_by. To use plot_genes_by_group,", 
               "please specify a group with more than one type. "))
  }
  names(cell_group) = colnames(cds)
  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>% 
    dplyr::summarize(mean = mean(log(Expression + pseudocount)), 
                     percentage = sum(Expression > lower_threshold)/length(Expression))
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min, 
                        ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max, 
                        ExpVal$mean)
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, "gene_short_name"]
  res <- reshape2::dcast(ExpVal[, 1:4], Group ~ Gene, value.var = colnames(ExpVal)[2 + 
                                                                                     major_axis])
  group_id <- res[, 1]
  res <- res[, -1]
  row.names(res) <- group_id
  if (ordering_type == "cluster_row_col") {
    row_dist <- stats::as.dist((1 - stats::cor(t(res)))/2)
    row_dist[is.na(row_dist)] <- 1
    col_dist <- stats::as.dist((1 - stats::cor(res))/2)
    col_dist[is.na(col_dist)] <- 1
    ph <- pheatmap::pheatmap(res, useRaster = T, cluster_cols = TRUE, 
                             cluster_rows = TRUE, show_rownames = F, show_colnames = F, 
                             clustering_distance_cols = col_dist, clustering_distance_rows = row_dist, 
                             clustering_method = "ward.D2", silent = TRUE, filename = NA)
    ExpVal$Gene <- factor(ExpVal$Gene, levels = colnames(res)[ph$tree_col$order])
    ExpVal$Group <- factor(ExpVal$Group, levels = row.names(res)[ph$tree_row$order])
  }
  else if (ordering_type == "maximal_on_diag") {
    order_mat <- t(apply(res, major_axis, order))
    max_ind_vec <- c()
    for (i in 1:nrow(order_mat)) {
      tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
      max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
    }
    max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]
    if (major_axis == 1) {
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(markers), 
                                            max_ind_vec))
      ExpVal$Gene <- factor(ExpVal$Gene, levels = dimnames(res)[[2]][max_ind_vec])
    }
    else {
      max_ind_vec <- c(max_ind_vec, setdiff(1:length(unique(exprs_mat$Group)), 
                                            max_ind_vec))
      ExpVal$Group <- factor(ExpVal$Group, levels = dimnames(res)[[1]][max_ind_vec])
    }
  }
  else if (ordering_type == "none") {
    ExpVal$Gene <- factor(ExpVal$Gene, levels = markers)
  }
  if (!is.null(gene_alias_vector)){
    ExpVal$Gene <- recode(ExpVal$Gene, !!!gene_alias_vector)
  }
  if (flip_percentage_mean) {
    g <- ggplot(ExpVal, aes(y = Gene, x = Group)) + geom_point(aes(colour = percentage, 
                                                                   size = mean)) + viridis::scale_color_viridis(name = "proportion") + 
      scale_size(name = "log(mean + 0.1)", range = c(0, 
                                                     max.size))
  }
  else {
    g <- ggplot(ExpVal, aes(y = Gene, x = Group)) + geom_point(aes(colour = mean, 
                                                                   size = percentage)) + viridis::scale_color_viridis(name = "log(mean + 0.1)") + 
      scale_size(name = "proportion", range = c(0, max.size))
  }
  if (group_cells_by == "cluster") {
    g <- g + xlab("Cluster")
  }
  else if (group_cells_by == "partition") {
    g <- g + xlab("Partition")
  }
  else {
    g <- g + xlab(group_cells_by)
  }
  g <- g + ylab("Gene") + theme(axis.text.x = element_text(angle = 30, 
                                                                                  hjust = 1))
  if (axis_order == "marker_group") {
    g <- g + coord_flip()
  }
  g
}

revigo_bubbles <- function(one.data,dispensability, pval){
  p1 <- ggplot( data = one.data );
  p1 <- p1 + geom_point( aes( x = plot_X, y = plot_Y, fill = -log10_p_value, size = plot_size), alpha = I(0.4), shape = 21 )# + scale_size_area(max_size = 1);
  p1 <- p1 + scale_fill_viridis_c();
  p1 <- p1 + scale_size( range=c(3, 15), guide = F);
  ex <- one.data [ one.data$dispensability < dispensability , ];
  ex <- ex [ex$log10_p_value < pval,];
  p1 <- p1 + ggrepel::geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = "black", alpha = 1, size = 5);
  p1 <- p1 + labs (x = "semantic space x", y = "semantic space y", size = "GO term frequency", fill = "-log10 p");
  p1 <- p1 + theme(legend.key = element_blank()) ;
  one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
  one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
  p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/5,max(one.data$plot_X)+one.x_range/5);
  p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/5,max(one.data$plot_Y)+one.y_range/5);
  #p1 <- p1 + theme_cowplot(font_size = 16) 
  return(p1)
}


get_mod_agg_values <- function(cds) {
  mod_values <- broom::tidy(normalized_counts(cds = cds, norm_method = "log", pseudocount = 1)) %>%
    group_by(column) %>% 
    summarise(mod_value = sum(value))
  return(mod_values %>% select(cell_uid = column, mod_value))
}

# replace colData

replace_colData <- function(cds, new_colData_tibble) {
  new_colData_df <- as.data.frame(new_colData_tibble)
  rownames(new_colData_df) <- new_colData_df$colData_rownames
  new_colData_df$colData_rownames <- NULL
  new_cds <- new_cell_data_set(
    expression_data = exprs(cds),
    gene_metadata = rowData(cds),
    cell_metadata = new_colData_df
  )
}
