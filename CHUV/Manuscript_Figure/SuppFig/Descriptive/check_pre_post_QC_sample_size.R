
# GeoMx -------------------------------------------------------------------
geo_breast_raw <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/GeoMx/breast_raw/Data_breast.rds")
geo_lung_raw <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/GeoMx/lung_raw/Data_lung.rds")
geo_dlbcl_raw <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/GeoMx/dlbcl_raw/Data_dlbcl.rds")

geo_breast_ruv <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/breast_spe_ruv.rds")
geo_lung_ruv <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/lung_spe_ruv.rds")
geo_dlbcl_ruv <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/dlbcl_spe_ruv.rds")

# Visium -------------------------------------------------------------------
vis_raw_breast <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/breast_raw/SPE/vis_raw_all_breast.rds")
vis_qcd_breast <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/breast/spotclean/Results/Breast-merge-SCTpostSpotClean.rds")

vis_raw_lung <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/lung_raw/SPE/vis_raw_all_lung.rds")
vis_qcd_lung <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/lung/spotclean/Results/Lung-merge-SCTpostSpotClean.rds")

dim(readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/dlbcl_raw/SPE/DLBCL_1.rds"))
dim(readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/dlbcl_raw/SPE/DLBCL_2.rds"))
dim(readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/dlbcl_raw/SPE/DLBCL_3.rds"))
dim(readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/dlbcl_raw/SPE/DLBCL_4.rds"))
dim(readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/dlbcl_raw/SPE/DLBCL_5.rds"))
dim(readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/dlbcl_raw/SPE/DLBCL_6.rds"))

vis_qcd_dlbcl <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/dlbcl/spotclean/Results/Dlbcl-merge-SCTpostSpotClean.rds")


# Chromium -----------------------------------------------------------------
chrom_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
 
chrom_breast_qcd <- readRDS(file.path(chrom_path, "chrom_breast.rds"))
chrom_lung_qcd <- readRDS(file.path(chrom_path, "chrom_lung.rds"))
chrom_dlbcl_qcd <- readRDS(file.path(chrom_path, "chrom_dlbcl.rds"))

