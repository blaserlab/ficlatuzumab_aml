source("00_bwb20200729.R")


revigo_bubbles <- function(revigo.data,dispensability, pval, plot_text_size){
  revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
  one.data <- data.frame(revigo.data);
  names(one.data) <- revigo.names;
  one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
  one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
  one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
  one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
  one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
  one.data$frequency <- as.numeric( as.character(one.data$frequency) );
  one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
  one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
  
  p1 <- ggplot( data = one.data );
  p1 <- p1 + geom_point( aes( x = plot_X, y = plot_Y, fill = -log10_p_value, size = plot_size), alpha = I(0.4), shape = 21 )# + scale_size_area(max_size = 1);
  p1 <- p1 + scale_fill_viridis_c();
  p1 <- p1 + scale_size( range=c(1, 3), guide = F);
  ex <- one.data [ one.data$dispensability <= dispensability , ];
  ex <- ex [ex$log10_p_value < pval,];
  p1 <- p1 + ggrepel::geom_text_repel(data = ex,
                                      aes(plot_X, plot_Y, label = description),
                                      colour = "black",
                                      size = plot_text_size,
                                      min.segment.length = 0,
                                      segment.size = 0.25);
  p1 <- p1 + labs (x = "semantic space x", y = "semantic space y", size = "GO term frequency", fill = "-log10 p");
  p1 <- p1 + theme(legend.key = element_blank()) ;
  one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
  one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
  p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/5,max(one.data$plot_X)+one.x_range/5);
  p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/5,max(one.data$plot_Y)+one.y_range/5);
  #p1 <- p1 + theme_cowplot(font_size = 16) 
  return(p1)
}


#save.image.pigz("ficlatuzumab_memory_managed.RData",n.cores = 39)

