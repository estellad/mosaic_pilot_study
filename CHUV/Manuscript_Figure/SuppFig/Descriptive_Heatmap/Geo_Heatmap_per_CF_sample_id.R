library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")

# disease = "breast"
# disease = "lung"
# disease = "dlbcl"

## Geo
datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"

disease_list <- c("breast", "lung", "dlbcl")
for(disease in disease_list){
  spe_ruv <- readRDS(file.path(datapath, paste0(disease, "_spe_ruv.rds")))
  spe_ruv <- spe_ruv[, spe_ruv$cell_fraction != "PanCK-"]
  spe_ruv$cell_fraction <- ifelse(spe_ruv$cell_fraction == "Macro", "Macrophage", spe_ruv$cell_fraction)
  
  assign(paste0("spe_ruv_", disease), spe_ruv)
}


plotHeatmap <- function(spe_ruv, label = "cell_fraction", text = "cell-fraction",  
                        disease = "breast", cf_levels, palette = "rocket"){ 
  
  if(disease == "dlbcl"){
    df <- as.data.frame(table(spe_ruv$cell_fraction, spe_ruv$patient))
  }else{
    df <- as.data.frame(table(spe_ruv$cell_fraction, spe_ruv$section_id))
  }
  
  df <- df %>%
    dplyr::rename(cf = Var1,
                  sample = Var2,
                  Count = Freq)
  
  # chrom_sample_order <- c("B1", "B2", "B3", "B4", 
  #                         "L1", "L2", "L3", "L4",
  #                         "D1", "D2", "D3", "D4", "D5", "D6")
  geo_sample_order <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1", 
                        "L1_1", "L2_1", "L3_1", "L3_3", "L4_3",
                        "D1", "D2", "D3", "D4", "D5", "D6")
  chrom_sample_size <- data.frame(
    # chrom_sample = chrom_sample_order,
    sample = geo_sample_order,
    SampleSize = c(3937, 3937, 4727, 3477, 991,
          802, 1756, 17804, 17804, 15592,
          8936, 7097, 14375, 5335, 2485, 1485)
  )
  
  # Weight by sample size to compare smaller sample size fairly
  df2_ <- df %>%
    left_join(chrom_sample_size, by = "sample") %>%
    mutate(Weight = 1/SampleSize,
           WeightedCount = Count * Weight) %>%
    group_by(cf) %>%
    mutate(per = round(100 * WeightedCount/sum(WeightedCount), 0))
  df2_$sample <- factor(df2_$sample, levels = rev(levels(factor(df2_$sample))))
  df2_$cf <- factor(df2_$cf, levels = cf_levels)
  
  # Create the heatmap pt per label
  p <- ggplot(df2_, aes(x = cf, y = sample, fill = per)) + 
    geom_tile() + 
    geom_text(aes(x = cf, y = sample, label = paste0(per, "%"))) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left") +
    scale_fill_viridis_c(option = palette) + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5)) + 
    labs(fill = paste0("Per-", text, " Percentage"), title = "", y = "Sample") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
  
  return(p)
}


p_breast_sampling <- plotHeatmap(spe_ruv_breast, label = "cell_fraction", text = "AOI-label", disease = "breast", 
                              cf_levels = c("Malignant", "Other", "T cells", "Macrophage"), palette = "cividis") 
p_lung_sampling <- plotHeatmap(spe_ruv_lung, label = "cell_fraction", text = "AOI-label", disease = "lung", 
                            cf_levels = c("Malignant", "Other", "T cells", "Macrophage"), palette = "cividis")
p_dlbcl_sampling <- plotHeatmap(spe_ruv_dlbcl, label = "cell_fraction", text = "AOI-label", disease = "dlbcl", 
                             cf_levels = c("B cells", "Other", "T cells", "Macrophage"), palette = "cividis") 

p_final <- ((p_breast_sampling + theme(axis.title.x = element_blank(), plot.margin = margin(0, 2, 0, 0, unit = "cm"))) | 
  (p_lung_sampling + theme(axis.title = element_blank(), plot.margin = margin(0, 2, 0, 2, unit = "cm"))) | 
  (p_dlbcl_sampling + theme(axis.title = element_blank(), plot.margin = margin(0, 0, 0, 2, unit = "cm")))) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 8),        
        legend.title = element_text(size = 10))

p_final <- p_final +
  plot_annotation(title = "GeoMx AOI Sampling") & 
  theme(plot.title = element_text(hjust = 0.53))

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Descriptive_Heatmap_Geo_Chrom"
pdf(file.path(figpath, "Geo_Pt_AOI.pdf"), width = 10, height = 4)
print(p_final)
dev.off()



