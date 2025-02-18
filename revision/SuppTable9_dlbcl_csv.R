library(Seurat)
library(dplyr)

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

vispath <- "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/dlbcl/spotclean/Results/"
vis <- readRDS(file.path(vispath, "Dlbcl-merge-SCTpostSpotClean.rds"))

geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/dlbcl_seu_ruv.rds")


# -------------------------------------------------------------------------
## Seurat decontaminated assay SoupX, DE with data slot of RNA
df_save <- data.frame(t(as.matrix(chrom@assays$SoupX@data)), CT = chrom$level1_5_immune_tumor_subtypes)
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppTabS9_chrom.csv")

df_save <- data.frame(t(as.matrix(vis@assays$SCT@data)), CT = vis$annot)
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppTabS9_vis.csv")

df_save <- data.frame(t(as.matrix(geo@assays$originalexp$data)), CT = geo$clusters)
write.csv(df_save, "/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/SourceData/SuppTabS9_geo.csv")
