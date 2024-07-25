# -------------------------------------------------------------------------
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5/B3"

# noSpotClean + logNorm + BayesSpace
sample_name = "B3_2"
disease = "breast"
sce <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered.rds") # noSpotClean + logNorm + BayesSpace

source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
set.seed(100)
sce <- runUMAP(sce, assay.type = "logcounts")
set.seed(100)
sce <- runTSNE(sce, assay.type = "logcounts")

sce$spatial.cluster <- as.factor(sce$spatial.cluster)
plotDimRed(sce, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster") |
  plotDimRed(sce, type = "TSNE", annotate = "spatial.cluster", text_by = "spatial.cluster")

sce$pxl_col_in_fullres <- NULL; sce$pxl_row_in_fullres <- NULL

p <- plotSpots(sce, annotate = "spatial.cluster", text_by = "spatial.cluster", y_reverse = TRUE,
               pt.size = 0.5)

plot_title = "Vis_B3_2_Spatial_clus.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 2,
    height = 3.5)
print(p)
dev.off()

table(sce$Region, sce$spatial.cluster)

# Save plots --------------------------------------------------------------
p <- plotDimRed(sce, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster")
plot_title = "Vis_B3_2_UMAP.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 11,
    height = 5)
print(p)
dev.off()


# Merge DE ----------------------------------------------------------------
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/15_Visium_Breast_Volcano/Tumor_Subtype_DE/00_Seurat_DE_heatmap_ClusPlot_by_Cluster_helper.R")
p <- clusterPlot_by_cluster(sce, type = "spot", sample_name = sample_name)

plot_title = "Vis_B3_2_Spatial_clus_by_clus.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 7,
    height = 10)
print(p)
dev.off()



plotDimRed(sce, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster") + facet_wrap(~ spatial.cluster)

sce$spatial.cluster_merge <- ifelse(sce$spatial.cluster %in% c(5, 9), "5_9", sce$spatial.cluster)
sce$spatial.cluster_merge2 <- ifelse(sce$spatial.cluster %in% c(1, 5, 9), "1_5_9", sce$spatial.cluster)

# > table(sce$Tu_Consensus, sce$spatial.cluster159reclus)
# 
#        10  11  12  13  14 159_1 159_2   2   3   4   6   7   8
# FALSE 104 127  59  75  82   258   142  33   1 152 140  31  65
# TRUE    3   0   4   1 143   392   244   3   0  54   0  38   5

sce$spatial.cluster159reclus_tu_consensus <- case_when(sce$spatial.cluster159reclus == "14" & sce$Tu_Consensus ~ "14", 
                                                       sce$spatial.cluster159reclus == "14" & !sce$Tu_Consensus ~ "14_TME", 
                                                       sce$spatial.cluster159reclus == "4" ~ "14_TME",
                                                       sce$spatial.cluster159reclus == "159_1" & sce$Tu_Consensus ~ "159_1", 
                                                       sce$spatial.cluster159reclus == "159_1" & !sce$Tu_Consensus ~ "159_1_TME",
                                                       sce$spatial.cluster159reclus == "159_2" & sce$Tu_Consensus ~ "159_2",
                                                       sce$spatial.cluster159reclus == "159_2" & !sce$Tu_Consensus ~ "159_2_TME",
                                                       .default = sce$spatial.cluster159reclus)

sce$spatial.cluster_tu_consensus <- case_when(sce$spatial.cluster159reclus == "14" & sce$Tu_Consensus ~ "14", 
                                              sce$spatial.cluster159reclus == "14" & !sce$Tu_Consensus ~ "14_TME", 
                                              sce$spatial.cluster159reclus == "4" ~ "14_TME",
                                              sce$spatial.cluster159reclus == "159_1" & sce$Tu_Consensus ~ "159", 
                                              sce$spatial.cluster159reclus == "159_1" & !sce$Tu_Consensus ~ "159_TME",
                                              sce$spatial.cluster159reclus == "159_2" & sce$Tu_Consensus ~ "159",
                                              sce$spatial.cluster159reclus == "159_2" & !sce$Tu_Consensus ~ "159_TME",
                                              .default = sce$spatial.cluster159reclus)

# TME clus 14: Tu consensus FALSE in clus 14 + entire clus 4
# TME clus 159: Tu consensus FALSE in clus 159_1 and clus 159_2

save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B1_2"
DE_result <- seurat_cluster_DE(sce, clusters = c("1_5_9", "14"), cluster_col="spatial.cluster_merge2", n_markers = 40)
contrast <- "1&5&9_14"

# DE_result <- seurat_cluster_DE(sce, cluster_col="spatial.cluster159reclus_tu_consensus", clusters = c("14", "159_1", "159_2"), n_markers = 100)
DE_result <- seurat_cluster_DE(sce, cluster_col="spatial.cluster_tu_consensus", clusters = c("14", "159"), n_markers = 100)
DE_result2 <- seurat_cluster_DE(sce, cluster_col="spatial.cluster_tu_consensus", clusters = c("14_TME", "159_TME"), n_markers = 100)


write.csv(DE_result, file.path(save_path_DE, paste0(contrast, ".csv")))
head(DE_result)
