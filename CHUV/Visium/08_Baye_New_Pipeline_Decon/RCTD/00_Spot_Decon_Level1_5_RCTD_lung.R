args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])
print(i)

library(SpatialExperiment)
library(spacexr)

# disease = "breast"
disease = "lung"
# disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

## Get chromium by indication ------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/" # breast / # lung
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

resultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/RCTD_vis/"

vis_rawmat <- "counts"

# for(i in 1:nsamples){
# Chrom QCed --------------------------------------------------------------
if(disease == "breast"){
  chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
}else if(disease == "lung"){
  chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))
}else{ # DLBCL
  chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
}

cell_types <- as.factor(chrom$level1_5_immune); names(cell_types) <- colnames(chrom)
chrom_mat <- chrom@assays$RNA@counts 

## Get Spot gene expr by sample -----------------------------------------
sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))

# vis: prep coords and count matrix
coords <- data.frame(x = spatialCoords(sce)[, 1], y = spatialCoords(sce)[, 2]); rownames(coords) <- colnames(sce)
counts_vis <- assay(sce, vis_rawmat)

ref <- Reference(chrom_mat, cell_types)
puck <- SpatialRNA(coords, counts_vis)
# RCTD --------------------------------------------------------------------
myRCTD <- create.RCTD(puck, ref, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

norm_weights <- normalize_weights(myRCTD@results$weights)
RCTD_results <- data.frame(as(norm_weights, "matrix"))
write.csv(RCTD_results, paste0(resultpath, save_names[i], "_spot_level1_5_immune_decon_RCTD.csv"))



