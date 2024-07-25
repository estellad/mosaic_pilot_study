library(SpatialDecon)
library(dplyr)
library(Seurat)
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
sigmatpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5/"
healthy_cells <- c("B", "Epithelia", "Fibro_Muscle", "Granulocyte", "Mast_cell", "Myeloid", "T_NK", "Vessel")

getSigmat <- function(disease){
  if(disease == "breast"){
    chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
    sigmat_name = "sigmat_breast_add_lung_healthy.rds"
  }else if(disease == "lung"){
    chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds")) 
    sigmat_name = "sigmat_lung_add_breast_healthy.rds"
  }else{
    chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds")); chrom$patient <- chrom$sample_id
    sigmat_name = "sigmat_dlbcl.rds"
  }
  
  annot_df <- data.frame(
    CellID = colnames(chrom),
    LabeledCellType = chrom$level1_5
  )
  
  custom_mtx <- create_profile_matrix(mtx = chrom@assays$RNA@counts,            # cell x gene count matrix
                                      cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
                                      cellTypeCol = "LabeledCellType",  # column containing cell type
                                      cellNameCol = "CellID",           # column containing cell ID/name
                                      matrixName = "custom_mini_colon", # name of final profile matrix
                                      outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                      normalize = FALSE,                # Should data be normalized? 
                                      minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                      minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                      scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                      discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.
  
  saveRDS(custom_mtx, file.path(sigmatpath, sigmat_name))
  p <- heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"), labRow = NA, margins = c(10, 5), cexCol = 0.7)
  print(p)
  
  return(custom_mtx)
}


getSigmat("breast")
getSigmat("lung")
getSigmat("dlbcl")

