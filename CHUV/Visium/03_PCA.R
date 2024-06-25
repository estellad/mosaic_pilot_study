disease_list = c("lung", "breast", "dlbcl")

for(d in 1:3){
  disease = disease_list[d]
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  
  # PCA ---------------------------------------------------------------------
  for(i in 1:nsamples){
    print(save_names[i])
    sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
    
    sce$pxl_col_in_fullres <- spatialCoords(sce)[, 1]
    sce$pxl_row_in_fullres <- spatialCoords(sce)[, 2]
    
    set.seed(42)
    sce <- spatialPreprocess(sce, n.PCs = 50, log.normalize = TRUE)
    
    foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
    save_bs_path <- paste0(baye_savepath, foldername, "/", save_names[i], "/")
    saveRDS(sce, file.path(save_bs_path, paste0(save_names[i], "_PCAed.rds")))
  }
}


