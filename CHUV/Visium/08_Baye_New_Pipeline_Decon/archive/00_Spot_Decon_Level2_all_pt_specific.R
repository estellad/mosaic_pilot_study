# disease = "lung"
# disease = "breast"
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
save_bs_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/"

vis_rawmat <- "counts"

i = 1
for(i in 1:nsamples){
  if(disease == "breast"){
    chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
    tumor_classes <- c("Tu_B1", "Tu_B3", "Tu_B4")
  }else if(disease == "lung"){
    chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))
    tumor_classes <- c("Tu_L1", "Tu_L2", "Tu_L3", "Tu_L4")
  }else{ # DLBCL
    chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

    tumor_classes <- c("Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6")
  }
  
  chrom$level1_5 <- chrom$Level2
  chrom$level1_5 <- ifelse(chrom$level1_5 %in% c("Fibro_Muscle", "Vessel"), "Stroma", chrom$level1_5)
  chrom$level1_5 <- as.factor(chrom$level1_5)
  
  chrom$Level2 <- ifelse(chrom$Level2 == "Fibro_Muscle", "Fibro_muscle", chrom$Level2)
  
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
  
  
  ## Get Spot gene expr by sample -----------------------------------------
  sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
  
  save_bs_path_sample <- paste0(save_bs_path, foldername, "/", save_names[i], "/")
  
  ## Decon -------------------------------------------------------------------
  # Ref: prep count matrix and cell anno
  cell_types <- data.frame(
    cellID = colnames(chrom),
    cellType = as.factor(chrom@meta.data[["level1_5"]]), 
    sampleInfo = "sample1")
  rownames(cell_types) <- colnames(chrom)
  # chrom_mat <- chrom@assays$SoupX@counts # decont chromium 
  chrom_mat <- chrom@assays$RNA@counts # raw chromium: use raw chrom for decont
  
  # sce: prep coords and count matrix
  coords <- data.frame(x = spatialCoords(sce)[, 1], y = spatialCoords(sce)[, 2]); rownames(coords) <- colnames(sce)
  counts_sce <- assay(sce, vis_rawmat)
  
  CARD_obj = createCARDObject(
    sc_count = chrom_mat,
    sc_meta = cell_types,
    spatial_count = counts_sce,
    spatial_location = coords,
    ct.varname = "cellType",
    ct.select = unique(cell_types$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5) 
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  
  results = data.frame(CARD_obj@Proportion_CARD)
  write.csv(results, paste0(save_bs_path_sample, save_names[i], "_spot_level1_5_decon.csv"))
}










