# Unify cell fraction names -----------------------------------------------
## GeoMx
read_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_test_norm/"
disease = "lung" # breast
spe <- readRDS(file.path(read_path, paste0(disease, ".rds")))
spe$cell_fraction <- ifelse(spe$cell_fraction == "T_cells", "T cells", spe$cell_fraction)
saveRDS(spe, file.path(read_path, paste0(disease, ".rds")))

disease = "dlbcl"
spe <- readRDS(file.path(read_path, paste0(disease, ".rds")))
spe <- spe[, !is.na(spe$cell_fraction)] # 136 samples 
# spe$cell_fraction_original <- spe$cell_fraction
spe$cell_fraction <- case_when(spe$cell_fraction %in% c("B lymphocytes", "Moslty B lymphocytes", "Mostly B lymphocytes", "Mostly B lymphocytes, with possibly some endothelial cells (vessels)", "Mostly B lymphocytes, with some T lymphocytes", "Mostly B lymphocytes, with endothelial cells (vessels)", "Mostly B lymphocytes, with some endothelial cells (vessels)", "Mostly B lymphocytes, with some endothelial cells (vessels)") ~ "B cells",
                               spe$cell_fraction %in% c("Epithelial cells (glands)", "Endothelial cells (vessel)", "Mostly epithelial cells (glands)") ~ "Other", 
                               spe$cell_fraction %in% c("Macrophages", "Moslty macrophages", "Mostly macrophages") ~ "Macro",
                               spe$cell_fraction %in% c("Mostly T lymphocytes", "T lymphocytes", "Mostly T lymphocytes, with B lymphocytes and endothelial cells (vessels)", 
                                                        "Mostly T lymphocytes, with endothelial cells (vessels)") ~ "T cells")
saveRDS(spe, file.path(read_path, paste0(disease, ".rds")))