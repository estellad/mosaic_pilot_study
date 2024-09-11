library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
library(patchwork)

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Descriptive_Heatmap_Geo_Chrom"

# disease = "breast"; seu <- readRDS(file.path(chrompath, "chrom_breast.rds"))
# disease = "lung"; seu <- readRDS(file.path(chrompath, "chrom_lung.rds"))
# disease = "dlbcl"; seu <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

disease_list <- c("breast", "lung", "dlbcl")

for(disease in disease_list){
  assign(paste0("seu_", disease), readRDS(file.path(chrompath, paste0("chrom_", disease, ".rds"))))
}


# -------------------------------------------------------------------------
plotHeatmap <- function(seu, label = "Harmonised_Level4", text = "cell-type", 
                        disease = "breast", palette = "inferno"){ 
  if(disease == "dlbcl"){
    seu$sample_id <- paste0(substr(seu$sample_id, 1, 1), substr(seu$sample_id, 7, 7))
  }
  
  df <- as.data.frame(table(seu@meta.data[[label]], seu$sample_id))
  df <- df %>%
    dplyr::rename(ct = Var1,
                  sample = Var2,
                  Count = Freq)
  
  # Weight by sample size to compare smaller sample size fairly
  df2_ <- df %>%
    group_by(sample) %>%
    mutate(SampleSize = sum(Count)) %>%
    ungroup() %>%
    mutate(SampleSize = ifelse(sample == "B2", 4727, SampleSize)) %>%
    mutate(Weight = 1/SampleSize,
           WeightedCount = Count * Weight) %>%
    group_by(ct) %>%
    mutate(per = round(100 * WeightedCount/sum(WeightedCount), 0))
  df2_$sample <- factor(df2_$sample, levels = rev(levels(factor(df2_$sample))))
  
  # Create the heatmap pt per label
  p <- ggplot(df2_, aes(x = ct, y = sample, fill = per)) + 
  # df <- as.data.frame(table(seu@meta.data[[label]], seu$sample_id))
  # df2 <- df %>%
  #   group_by(Var1) %>%
  #   mutate(countT= sum(Freq)) %>%
  #   mutate(per=round(100*Freq/countT,0))
  # 
  # data_frame <- data.frame(label = df2$Var1, sample = df2$Var2, per = df2$per)
  # data_frame$label <- factor(data_frame$label, levels = rev(levels(factor(data_frame$label))))
  
  # # Create the heatmap pt per label
  # p <- ggplot(data_frame, aes(x = sample, y = label, fill = per)) + 
    geom_tile() + 
    geom_text(aes(x = ct, y = sample, label = paste0(per, "%")), size = 3) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left") +
    scale_fill_viridis_c(option = palette) + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5)) + 
    labs(fill = paste0("Per-", text, " Percentage"), title = element_blank(), y = "Sample") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(7, "cm"), barheight = unit(0.5, "cm")))
  
}


# -------------------------------------------------------------------------
p_breast_decon <- plotHeatmap(seu_breast, label = "Harmonised_Level4", text = "cell-type", disease = "breast", palette = "inferno")
p_lung_decon <- plotHeatmap(seu_lung, label = "Harmonised_Level4", text = "cell-type", disease = "lung", palette = "inferno")
p_dlbcl_decon <- plotHeatmap(seu_dlbcl, label = "Harmonised_Level4", text = "cell-type", disease = "dlbcl", palette = "inferno")

p_final <- (p_breast_decon + theme(axis.title.x = element_blank())) / 
  (p_lung_decon + theme(axis.title.x = element_blank())) / 
  (p_dlbcl_decon + theme(axis.title.x = element_blank())) + 
  plot_layout(guides = "collect", heights = c(0.5, 0.5, 0.6)) & 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 8),        
        legend.title = element_text(size = 10))

p_final <- p_final +
  plot_annotation(title = "Level 4 Annotation") & 
  theme(plot.title = element_text(hjust = 0.53))

pdf(file.path(figpath, "Chrom_Pt_Level4_Heatmap.pdf"), width = 12, height = 13)
print(p_final)
dev.off()

