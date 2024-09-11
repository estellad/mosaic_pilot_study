disease = "breast"
#disease = "lung"
# disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
if(disease == "breast"){
  chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
}else if(disease == "lung"){
  chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))
}else{ # DLBCL
  chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
}

save_bs_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/"

vis_rawmat <- "counts"

for(i in 1:nsamples){
  print(i)
  
  # chrom_pt_tumor <- chrom_classes[grepl(paste0("Tu_L", i), chrom_classes)]
  # chrom <- chrom[, chrom$Harmonised_Level4 %in% c(chrom_immune_include, chrom_pt_tumor)]
  
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










