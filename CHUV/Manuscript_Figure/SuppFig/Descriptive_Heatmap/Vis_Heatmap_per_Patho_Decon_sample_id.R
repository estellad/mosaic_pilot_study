library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Gallery_Patho_Decon_Pt_Heatmap"

# disease = "breast"
# disease = "lung"
# disease = "dlbcl"

disease_list <- c("breast", "lung", "dlbcl")

for(disease in disease_list){
  # After spotclean -----------------------------------------------------
  datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean/Results")
  save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
  
  assign(paste0("seu_", disease), readRDS(paste0(datapath.spotclean, save_rds_name)))
}


# -------------------------------------------------------------------------
plotHeatmap <- function(seu, label = "Region", text = "region", 
                        disease = "breast", palette = "viridis"){ # label = "Level4_decon_max", text = "cell-type"
  # Heatmap of Region versus sample id ----------------------------------
  if(disease == "dlbcl"){
    seu$sample_id <- paste0(substr(seu$sample_id, 1, 1), substr(seu$sample_id, 7, 7))
  }

  df <- as.data.frame(table(seu@meta.data[[label]], seu@meta.data[["sample_id"]]))
  df2 <- df %>%
    group_by(Var1) %>%
    mutate(countT= sum(Freq)) %>%
    mutate(per=round(100*Freq/countT,0))
  
  data_frame <- data.frame(label = df2$Var1, sample = df2$Var2, per = df2$per)
  data_frame$sample <- factor(data_frame$sample, levels = rev(levels(factor(data_frame$sample))))
  
  # Create the heatmap pt per label
  p <- ggplot(data_frame, aes(x = label, y = sample, fill = per)) + 
    geom_tile() + 
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left") +
    scale_fill_viridis_c(option = palette) + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5)
          ) + 
    labs(fill = paste0("Per-", text, " Percentage"), title = "", y = "Sample") + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                                 barwidth = unit(5, "cm"), barheight = unit(0.3, "cm")))
  
  return(p)
}


# -------------------------------------------------------------------------
p_breast_patho <- plotHeatmap(seu_breast, label = "Region", text = "region", disease = "breast", palette = "viridis")
p_lung_patho <- plotHeatmap(seu_lung, label = "Region", text = "region", disease = "lung", palette = "viridis")
p_dlbcl_patho <- plotHeatmap(seu_dlbcl, label = "Region", text = "region", disease = "dlbcl", palette = "viridis")

p_final <- (p_breast_patho + theme(axis.title.x = element_blank()) | 
              p_lung_patho + theme(axis.title = element_blank()) | 
              p_dlbcl_patho + theme(axis.title = element_blank())) + 
  plot_layout(guides = "collect", widths = c(1.2, 1.1, 0.7)) & 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 8),        
        legend.title = element_text(size = 10))

p_final <- p_final +
  plot_annotation(title = "Pathology Annotation") & 
  theme(plot.title = element_text(hjust = 0.53))

pdf(file.path(figpath, "Vis_Pt_Patho_Heatmap.pdf"), width = 8, height = 5)
print(p_final)
dev.off()


# -------------------------------------------------------------------------
p_breast_decon <- plotHeatmap(seu_breast, label = "Level4_decon_max", text = "cell-type", disease = "breast", palette = "plasma")
p_lung_decon <- plotHeatmap(seu_lung, label = "Level4_decon_max", text = "cell-type", disease = "lung", palette = "plasma")
p_dlbcl_decon <- plotHeatmap(seu_dlbcl, label = "Level4_decon_max", text = "cell-type", disease = "dlbcl", palette = "plasma")

p_final <- (p_breast_decon + theme(axis.title.x = element_blank())) / 
              (p_lung_decon + theme(axis.title.x = element_blank())) / 
              (p_dlbcl_decon + theme(axis.title.x = element_blank())) + 
  plot_layout(guides = "collect", heights = c(0.5, 0.5, 0.6)) & 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 8),        
        legend.title = element_text(size = 10))

p_final <- p_final +
  plot_annotation(title = "Deconvolution majority vote") & 
  theme(plot.title = element_text(hjust = 0.53))

pdf(file.path(figpath, "Vis_Pt_Decon_Heatmap.pdf"), width = 10, height = 14)
print(p_final)
dev.off()


# -------------------------------------------------------------------------
# seu <- seu[, ((colnames(seu) %in% paste0(sce_B1_2$Barcode, "_1")) & (seu$sample_id == "B1_2")) | seu$sample_id %in% c("B1_4", "B2_2", "B3_2", "B4_2")]
# table(seu$Region)
# 
# barcode_order <- data.frame(Barcode = colnames(seu)[seu$Region == "FALSE"])
# 
# CD <- data.frame(Barcode = paste0(sce_B1_2$Barcode, "_1"),
#                  Region = sce_B1_2$Region)
# 
# barcode_order <- barcode_order %>% left_join(CD)
# 
# seu$Region <- ifelse(seu$sample_id == "B1_2", barcode_order$Region, seu$Region)
# saveRDS(seu, paste0(savepath.spotclean, paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean_B1_2_qcd.rds")))
# seu <- readRDS(paste0(savepath.spotclean, paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean_B1_2_qcd.rds")))



