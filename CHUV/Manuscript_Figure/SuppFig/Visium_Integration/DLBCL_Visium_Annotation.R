library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)

disease = "dlbcl"

# After spotclean -----------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
UMAP_name = paste0("/", str_to_title(disease), "-UMAP-")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
savepath = savepath.spotclean

# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(savepath, save_rds_name))

DimPlot(seu, label = TRUE, group.by = "Level4_decon_max") + facet_wrap(~Level4_decon_max) + 
  theme(legend.position = "none")
DimPlot(seu, label = TRUE, group.by = "seurat_clusters") + facet_wrap(~seurat_clusters) + 
  theme(legend.position = "none")

# Plasma cells
FeaturePlot(seu, features = c("IGKC", "IGHG3", "IGHA", "JCHAIN", "IGHG1"))
# Epithelium
FeaturePlot(seu, features = c("MUCL3", "LCN2", "TFF1", "MUC5AC"))
# Stroma
FeaturePlot(seu, features = c("THBS1", "POSTN", "COL3A1", "COL1A2", "COL1A1"), ncol = 3)
# Vessels
FeaturePlot(seu, features = c("ACKR1", "MEOX2", "CCL14", "TSPAN7", "PLA1A", "VWF", "GIMAP4"), ncol = 4)

# Marcophage 
FeaturePlot(seu, features = c("CD68", "CD163", "CD14", "ITGAM", "MRC1", "MSR1", "SIRPA", "ARG1"))
FeaturePlot(seu, features = rownames(seu)[grepl("^MMP", rownames(seu))])
FeaturePlot(seu, features = rownames(seu)[grepl("^TGFB", rownames(seu))])


table# Tu_D1_LMO2
FeaturePlot(seu, features = c("LMO2"))
# Tu_D1_SMIM14
FeaturePlot(seu, features = c("SMIM14"))

# # Tu_D2_mito
# FeaturePlot(seu, features = c("mito?"))

# check marker gene from heatmap for clus 14 and 16
FeaturePlot(seu, features = c("SLC40A1"))

FeaturePlot(seu, features = c("MMP12"))

# Tu_D3_FAM3C
FeaturePlot(seu, features = c("FAM3C"))
# Tu_D4_BCL7A
FeaturePlot(seu, features = c("BCL7A"))
# Tu_D5_CCL22
FeaturePlot(seu, features = c("CCL22"))
# Tu_D6_BCL2
FeaturePlot(seu, features = c("BCL2"))

seu_sub <- seu[, seu$seurat_clusters == "18"]
table(seu_sub$sample_id, seu_sub$Level4_decon_max)


# wide format
test <- table(seu$Level4_decon_max, seu$seurat_clusters)
# long format with counts
test2 <- as.data.frame(test)

# Heatmap (decon x clus) --------------------------------------------------
data_frame <- data.frame(cluster = test2$Var2, decon = test2$Var1, counts = test2$Freq)
data_frame$decon <- factor(data_frame$decon, levels = rev(levels(factor(data_frame$decon))))

# Create the heatmap
ggplot(data_frame, aes(x = cluster, y = decon, fill = counts)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "", x = "Level 4 Deconvolution Majority Vote", y = "Integrated DLBCL Visium Seurat Clusters")


# -------------------------------------------------------------------------
# long format with percentage per-clus
test3 <- test2 %>%
  group_by(Var2) %>%
  mutate(countT= sum(Freq)) %>%
  mutate(per=round(100*Freq/countT,0))

data_frame <- data.frame(cluster = test3$Var2, decon = test3$Var1, per = test3$per)
data_frame$decon <- factor(data_frame$decon, levels = rev(levels(factor(data_frame$decon))))

# Create the heatmap
ggplot(data_frame, aes(x = cluster, y = decon, fill = per)) + 
  geom_tile() + 
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(fill = "Per-cluster\nPercentage", title = "", 
       x = "Integrated DLBCL Visium Seurat Clusters", y = "Level 4 Deconvolution Majority Vote") 


# wide format
library(tidyverse)
test5 <- test3 %>%
  select(Var1, Var2, per) %>%
  spread(Var2, per)



