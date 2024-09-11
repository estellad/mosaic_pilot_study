plotSpots_deconmax <- function(sce, 
                               annotate, 
                               pt.size,
                               palette){
  plt_df <- cbind.data.frame(colData(sce), spatialCoords(sce))
  p <- ggplot(plt_df, aes_string(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres", color = annotate)) + 
    geom_point(size = pt.size) + 
    coord_fixed() + 
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right") +
    scale_color_manual(values = palette) 
  
  p
}