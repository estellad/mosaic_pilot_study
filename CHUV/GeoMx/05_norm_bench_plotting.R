###############################################################################
#                                 With GeoDiff                                #
###############################################################################
save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_test_norm/")
spe <- readRDS(file.path(save_path, paste0(disease, "_geodiff_included.rds")))

# Plotting ----------------------------------------------------------------
norm_titles <- c("Raw", 
                 "log(y/10^6 + 1) - logCPM", 
                 "log(y + 1) - log1p", 
                 "log(y/s + 1) - Scran size factor", 
                 "log(y/s + 1) - Mean size factor", 
                 "TMM", "Q3", "Quantile", "scTransform", "GeoDiff")

## Color by slide 
for(i in 1:length(names(assays(spe)))){
  p <- plotRLExpr(spe, assay = names(assays(spe))[i], color = slide_name) + ggtitle(norm_titles[i]) 
  if(i != 5){p <- p + theme(legend.position = "none")}
  assign(paste0("p", i), p)
}

p_norm_bench <- (p1 | p3 |p2 | p4 | p5) /
  (p6 | p7 | p8 |p9 |p10) + 
  plot_annotation(title = paste0("GeoMx normalization method benchmarking (18676 genes, ", ncol(spe), " samples)"),
                  theme = theme(plot.title = element_text(hjust=0.5)))

# pdf(file.path(save_path, paste0(disease, "_geodiff_included_benchmark_colbyslide.pdf")),
#     height = 11, width = 6)
# print(p_norm_bench)
# dev.off()


## Color by patient 
for(i in 1:length(names(assays(spe)))){
  p <- plotRLExpr(spe, assay = names(assays(spe))[i], color = patient_id) + ggtitle(norm_titles[i]) # Breast lung
  if(i != 5){p <- p + theme(legend.position = "none")}
  assign(paste0("p", i), p)
}

p_norm_bench <- (p1 | p3 |p2 | p4 | p5) /
  (p6 | p7 | p8 |p9 |p10) + 
  plot_annotation(title = paste0("GeoMx normalization method benchmarking (18676 genes, ", ncol(spe), " samples)"),
                  theme = theme(plot.title = element_text(hjust=0.5)))

# pdf(file.path(save_path, paste0(disease, "_geodiff_included_benchmark_colbypt.pdf")),
#     height = 11, width = 6)
# print(p_norm_bench)
# dev.off()



###############################################################################
#                              Without GeoDiff                                #
###############################################################################
spe <- readRDS(file.path(save_path, paste0(disease, ".rds")))

## Color by slide
for(i in 1:length(names(assays(spe)))){
  p <- plotRLExpr(spe, assay = names(assays(spe))[i], color = slide_name) + ggtitle(norm_titles[i]) 
  if(i != 5){p <- p + theme(legend.position = "none")}
  assign(paste0("p", i), p)
}

p_norm_bench <-patchwork::wrap_plots(p1, p3, p2, p4, p5, 
                                    p6, p7, p8, p9, 
                                    ncol = 5) + 
  plot_annotation(title = paste0("GeoMx normalization method benchmarking (18676 genes, ", ncol(spe), " samples)"),
                  theme = theme(plot.title = element_text(hjust=0.5)))

## Color by patient
for(i in 1:length(names(assays(spe)))){
  p <- plotRLExpr(spe, assay = names(assays(spe))[i], color = patient_id) + ggtitle(norm_titles[i]) 
  if(i != 5){p <- p + theme(legend.position = "none")}
  assign(paste0("p", i), p)
}

p_norm_bench <-patchwork::wrap_plots(p1, p3, p2, p4, p5, 
                                     p6, p7, p8, p9, 
                                     ncol = 5) + 
  plot_annotation(title = paste0("GeoMx normalization method benchmarking (18676 genes, ", ncol(spe), " samples)"),
                  theme = theme(plot.title = element_text(hjust=0.5)))
