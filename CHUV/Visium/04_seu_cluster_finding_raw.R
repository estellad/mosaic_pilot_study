# disease = "breast" # raw n_clus: 16 16 16 14 16
# disease = "lung"   # raw n_clus: 11 10 12 15 15
# disease = "dlbcl"  # raw n_clus: 15 16 15 15 16 17 # D1 and D3 has hair, now imputed in tissue with raw and pathology
disease = "dlbcl"    # raw n_clus: 15 16 19 15 16 17
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

n_clus <- c()
for(i in 1:nsamples){  
  save_bs_path <- paste0(baye_savepath, foldername, "/", save_names[i], "/")
  sce <- readRDS(file.path(save_bs_path, paste0(save_names[i], "_PCAed.rds")))
  
  seu <- as.Seurat(sce,
                   counts = "counts",
                   data = NULL,
                   assay = NULL,
                   project = "SingleCellExperiment")

  seu <- FindNeighbors(seu, dims = 1:50, reduction = "PCA")
  seu <- FindClusters(seu, resolution = 2, cluster.name = "seu_spatial",
                      random.seed = 123)
  n_clus <- c(n_clus, nlevels(seu$seu_spatial))
}

print(n_clus)

