# disease = "breast"
# disease = "lung"
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/00_GeoMx_Paths.R")
save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_test_norm/")
fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/GeoMx_Normalization_Benchmark"

# Plotting ----------------------------------------------------------------
norm_titles <- c("Raw", 
                 expression(paste("logCPM -  log(",y,"/",10^{6}," + 1)")), 
                 "log1p -  log(y + 1)", 
                 "Scran size factor -  log(y/s + 1)", 
                 "Mean size factor -  log(y/s + 1)", 
                 "TMM", "Q3", "Quantile", "scTransform", "GeoDiff")

## In total 8 slides across three indications 
library(RColorBrewer)
n <- 8
color_slides <- brewer.pal(n, "Set2")
pie(rep(1,n), col=color_slides)

###############################################################################
#                                 With GeoDiff                                #
###############################################################################
spe <- readRDS(file.path(save_path, paste0(disease, "_geodiff_included.rds")))

palette <- list(
  breast = color_slides[1:2],
  lung = color_slides[3:5],
  dlbcl = color_slides[6:8]
)

## Color by slide 
for(i in 1:length(names(assays(spe)))){
  p <- plotRLExpr(spe, assay = names(assays(spe))[i], color = slide_name) + 
    ggtitle(norm_titles[i]) + xlab("Sample") + 
    scale_color_manual(values = palette[[disease]])
  if(!(i %in% c(6, 7, 8, 9, 10))){p <- p + xlab("")}
  if(!(i %in% c(1,6))){p <- p + ylab("")}
  if(i != 5){p <- p + theme(legend.position = "none")}
  assign(paste0("p", i), p)
}

p_norm_bench <- (p1 | p3 |p2 | p4 | p5) /
  (p6 | p7 | p8 |p9 |p10) + 
  plot_annotation(title = paste0("GeoMx normalization method benchmarking (18676 genes, ", disease, ", ", ncol(spe), " samples)"),
                  theme = theme(plot.title = element_text(hjust=0.5, size = 16)))

pdf(file.path(fig_path, paste0(disease, "_geodiff_included_benchmark_colbyslide.pdf")),
    height = 8, width = 20)
print(p_norm_bench)
dev.off()
