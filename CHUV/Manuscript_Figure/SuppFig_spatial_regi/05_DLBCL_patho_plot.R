disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
samples_for_registration <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6") # Visium
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")
patho_color <- dlbcl_patho_color
plot_savepath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/SuppFig_regi_spatial/With_legend/"

for(i in 1:length(samples_for_registration)){
  sce <- readRDS(file.path(save_path_qcd, paste0(samples_for_registration[i], "_qcd.rds")))
  table(sce$Region)
  # sce$Region <- ifelse(sce$Region == "Vessel", "Vessels", sce$Region)
  # sce <- sce[, !is.na(sce$Region)]
  # saveRDS(sce, file.path(save_path_qcd, paste0(samples_for_registration[i], "_qcd.rds")))

  
  # Patho plot --------------------------------------------------------------
  point_size = ifelse(samples_for_registration[i] %in% c("DLBCL_1", "DLBCL_2", "DLBCL_3"), 1.5, 2)
  plot_patho <- plotSpots(sce, in_tissue = NULL, annotate = "Region",
                          x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", pt.size = point_size) +
    scale_color_manual(values = patho_color[names(table(sce$Region))], 
                       na.value = "#d3d3d3") + labs(color = "Pathology")
  
  pdf(file = file.path(plot_savepath, paste0(samples_for_registration[i], "_patho.pdf")),
      width = 8,
      height = 7)
  print(plot_patho)
  dev.off()

}