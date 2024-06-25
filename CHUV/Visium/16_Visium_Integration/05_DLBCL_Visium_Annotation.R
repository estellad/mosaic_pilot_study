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



# Tu_D1_LMO2
FeaturePlot(seu, features = c("LMO2"))
# Tu_D1_SMIM14
FeaturePlot(seu, features = c("SMIM14"))
# Tu_D3_FAM3C
FeaturePlot(seu, features = c("FAM3C"))
# Tu_D5_CCL22
FeaturePlot(seu, features = c("CCL22"))

test <- table(seu$Level4_decon_max, seu$seurat_clusters)
test2 <- as.data.frame(test, ncol = ncol(test)) 
test3 <- test2 %>%
  group_by(Var2) %>%
  mutate(countT= sum(Freq)) %>%
  mutate(per=round(100*Freq/countT,0))

## Does not really work
# test4 <- test3 %>%
#   group_by(Var2) %>%
#   filter(per == max(per))

library(tidyverse)
test5 <- test3 %>%
  select(Var1, Var2, per) %>%
  spread(Var2, per)



