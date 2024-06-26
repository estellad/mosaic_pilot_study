read_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_test_norm/")
save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_Final_Qcd/")

# Lung -------------------------------------------------------------------
disease = "lung"
spe <- readRDS(file.path(read_path, paste0(disease, ".rds")))

spe <- spe[, spe$gene_detection_rate >= 0.1] # 122 -> 47
saveRDS(spe, file.path(save_path, paste0(disease, "_qcd.rds")))


# Breast -----------------------------------------------------------------
disease = "breast"
spe <- readRDS(file.path(read_path, paste0(disease, ".rds")))

spe <- spe[, spe$gene_detection_rate >= 0.1] # 117 -> 115
saveRDS(spe, file.path(save_path, paste0(disease, "_qcd.rds")))


# DLBCL ------------------------------------------------------------------
disease = "dlbcl"
spe <- readRDS(file.path(read_path, paste0(disease, ".rds")))

spe <- spe[, spe$gene_detection_rate >= 0.1] # 137 -> 132
saveRDS(spe, file.path(save_path, paste0(disease, "_qcd.rds")))

