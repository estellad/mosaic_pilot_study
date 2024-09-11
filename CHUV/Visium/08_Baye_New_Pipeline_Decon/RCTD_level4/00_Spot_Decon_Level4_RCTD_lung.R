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

vis_rawmat <- "counts"

# i = 4
# for(i in 1:nsamples){
if(disease == "breast"){
  chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
  tumor_classes <- c("Tu_B1_MUCL1", "Tu_B1_MUCL1_necrosis", "Tu_B1_MUCL1_transcription", "Tu_B3_CYP4F8", "Tu_B4_RHOB")
}else if(disease == "lung"){
  chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))
  tumor_classes <- c("Tu_L1_SFTPB",                  "Tu_L2_FXYD2",        "Tu_L3_G0S2",           "Tu_L3_G0S2_immune_signature", 
                     "Tu_L4_KRT17_immune_signature", "Tu_L4_KRT17_mucous", "Tu_L4_KRT17_necrosis", "Tu_L4_KRT17_neutrophil_signature")
}else{ # DLBCL
  chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
  tumor_classes <- c("Tu_D1_LMO2",   "Tu_D1_RGS13",    "Tu_D1_SMIM14", 
                     "Tu_D2_mito",
                     "Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD",
                     "Tu_D4_BCL7A",  "Tu_D4_PNN",      "Tu_D5_CCL22", "Tu_D6_BCL2")
}

print(i)
if(disease %in% c("breast", "lung")){ # (btw B2_2 does not have tumor class)
  chrom_other_pt_tumor <- tumor_classes[!grepl(substr(save_names[i], 1, 2), 
                                               tumor_classes)]
  if(save_names[i] == "B3_2"){
    fibro_remove = "Fibroblast"
  }else{ # all other samples
    fibro_remove = "Fibroblast_B3"
  }
  
  if(save_names[i] == "B2_2"){ # Pull other breast tumor class into tumor
    chrom$Harmonised_Level4 <- ifelse(chrom$Harmonised_Level4 %in% tumor_classes, "Tu_B2", chrom$Harmonised_Level4)
  }
  to_remove = c(chrom_other_pt_tumor, fibro_remove)
}else{ # "dlbcl"
  chrom_other_pt_tumor <- tumor_classes[!grepl(paste0(substr(save_names[i], 1, 1), substr(save_names[i], 7, 7)), 
                                               tumor_classes)]
  to_remove = chrom_other_pt_tumor
}

chrom <- chrom[, !(chrom$Harmonised_Level4 %in% to_remove)]
chrom$Harmonised_Level4 <- droplevels(as.factor(chrom$Harmonised_Level4))
table(chrom$Harmonised_Level4)

cell_types <- as.factor(chrom$Harmonised_Level4); names(cell_types) <- colnames(chrom)
chrom_mat <- chrom@assays$RNA@counts 

## Get Spot gene expr by sample -----------------------------------------
sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))

# vis: prep coords and count matrix
coords <- data.frame(x = spatialCoords(sce)[, 1], y = spatialCoords(sce)[, 2]); rownames(coords) <- colnames(sce)
counts_vis <- assay(sce, vis_rawmat)

ref <- Reference(chrom_mat, cell_types)
puck <- SpatialRNA(coords, counts_vis)
# RCTD --------------------------------------------------------------------
myRCTD <- create.RCTD(puck, ref, max_cores = 2)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

norm_weights <- normalize_weights(myRCTD@results$weights)
RCTD_results <- data.frame(as(norm_weights, "matrix"))

result_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/RCTD"
write.csv(RCTD_results, file.path(result_path, paste0(save_names[i], "_spot_level4_RCTD.csv")))




