disease_list = c("breast", "lung")
d = 2

# for(d in 1:length(disease_list)){
disease = disease_list[d]
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")
if(disease == "breast"){
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  samples_for_registration <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2") # Visium
  geo_reg_names <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1")
  patho_color <- breast_patho_color
}else if(disease == "lung"){
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  samples_for_registration <- c("L1_2", "L2_2", "L3_2", "L4_2") # Visium
  geo_reg_names <- c("L1_1", "L2_1", "L3_1", "L4_1")
  patho_color <- lung_patho_color
}

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
regis_savepath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration/"
plot_savepath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/SuppFig_regi_spatial/With_legend/"

for(i in 1:length(samples_for_registration)){
  # sce <- readRDS(file.path(save_path_qcd, paste0(samples_for_registration[i], "_qcd.rds"))) - for patho only
  sce <- readRDS(file.path(paste0(regis_savepath, "/Visium_mapped_spot_obj"), paste0(samples_for_registration[i], "_mapped_spot.rds"))) # spot
  sce_subspot <- readRDS(file.path(paste0(regis_savepath, "/Visium_mapped_subspot_obj"), paste0(samples_for_registration[i], "_mapped_subspot.rds"))) # subspot
  
  if(samples_for_registration[i] != "B1_2"){
    sce$Region_mapped <- ifelse(is.na(sce$Cell_fraction), NA, sce$Region)
    sce_subspot$Region_mapped <- ifelse(is.na(sce_subspot$Cell_fraction), NA, sce_subspot$Region)
  }
  
  y_rev = ifelse(samples_for_registration[i] == "B1_4", FALSE, TRUE)
  # Plotting ----------------------------------------------------------------
  p1 <- plotSpots(sce, in_tissue = NULL, annotate = "Cell_fraction",
                  x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", pt.size = 1, y_reverse = y_rev) + 
    scale_color_manual(values = cell_fraction_color[names(table(sce$Cell_fraction))], 
                       na.value = "#d3d3d3") + labs(color = "Cell fraction")
  
  if(samples_for_registration[i] != "B1_2"){
    p2 <- plotSpots(sce, in_tissue = NULL, annotate = "Region_mapped",
                    x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", pt.size = 1, y_reverse = y_rev) + 
      scale_color_manual(values = patho_color[names(table(sce$Region))], 
                         na.value = "#d3d3d3") + labs(color = "Pathology")
  }
  

  p3 <- plotSpots(sce_subspot, in_tissue = NULL, annotate = "Cell_fraction",
                  x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", y_reverse = y_rev) + 
    scale_color_manual(values = cell_fraction_color[names(table(sce$Cell_fraction))],
                       na.value = "#d3d3d3") + labs(color = "Cell fraction")
  
  if(samples_for_registration[i] != "B1_2"){
    p4 <- plotSpots(sce_subspot, in_tissue = NULL, annotate = "Region_mapped",
                    x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", y_reverse = y_rev) +
      scale_color_manual(values = patho_color[names(table(sce_subspot$Region))], 
                         na.value = "#d3d3d3") + labs(color = "Pathology")
  }
  

  # # Patho plot --------------------------------------------------------------
  # if(samples_for_registration[i] != "B1_2"){
  #   point_size = ifelse(samples_for_registration[i] %in% c("L3_2", "L4_2", "B2_2", "B4_2"), 2, 
  #                       ifelse(samples_for_registration[i] == "B3_2", 2.5, 3))
  #   plot_patho <- plotSpots(sce, in_tissue = NULL, annotate = "Region",
  #                           x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", pt.size = point_size, y_reverse = y_rev) +
  #     scale_color_manual(values = patho_color[names(table(sce$Region))], 
  #                        na.value = "#d3d3d3") + labs(color = "Pathology")
  #   
  #   if(samples_for_registration[i] == "B1_4"){plot_patho <- plot_patho + coord_flip()}
  #   
  #   pdf(file = file.path(paste0(plot_savepath, "Pathology/"),
  #                        paste0(samples_for_registration[i], "_patho.pdf")),
  #       width = 8,
  #       height = 7)
  #   print(plot_patho)
  #   dev.off()
  # }
  

  # Regi plot ---------------------------------------------------------------
  if(samples_for_registration[i] == "B1_2"){
    plot_regi <- p1 / p3  + plot_annotation(title = samples_for_registration[i]) 
    plt_regi_width = 5; plt_regi_height = 6
  }else{
    if(samples_for_registration[i] == "B1_4"){
      plot_regi <- (p1 + coord_flip() | p2 + coord_flip())/ 
        (p3 + coord_flip() | p4 + coord_flip()) + plot_annotation(title = samples_for_registration[i]) 
    }else{
      plot_regi <- (p1 | p2)/ 
        (p3 | p4) + plot_annotation(title = samples_for_registration[i]) 
    }

    plt_regi_width = 10; plt_regi_height = 6
  }
  
  pdf(file = file.path(paste0(plot_savepath, "WithPanCK-"), 
                       paste0(samples_for_registration[i], "_regi.pdf")),
      width = plt_regi_width,
      height = plt_regi_height)
  print(plot_regi)
  dev.off()
  

  # Stats -------------------------------------------------------------------
  print(paste0(samples_for_registration[i], " Percent of registration: ", 
               "spot: ", round(100*(length(which(!is.na(sce$Cell_fraction)))/ncol(sce)), 1),
               "; subspot: ", round(100*(length(which(!is.na(sce_subspot$Cell_fraction)))/ncol(sce_subspot)), 1)
               ))
}


# breast_spot_pct_regi <- c(4,   4.9, 5.9, 4.5, 4.8)
# breast_subspot_pct_regi <- c(4.9, 5.7, 6.4, 6.3, 5.8)
# lung_spot_pct_regi <-   c(3.6, 3.3, 2.3, 0.4)
# lung_subspot_pct_regi <- c(6.2, 5, 3, 0.7)
# 
# sum(c(breast_spot_pct_regi, lung_spot_pct_regi))/9 # spot pct regi
# sum(c(breast_subspot_pct_regi, lung_subspot_pct_regi))/9 # subspot pct regi

# Patho plot for sample L1_4 only
disease = "lung"
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
sce <- readRDS(file.path(save_bs_path, paste0("L1_4", "_qcd.rds")))

plot_patho <- plotSpots(sce, in_tissue = NULL, annotate = "Region",
                        x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", pt.size = 3,
                        y_reverse = FALSE) +
  scale_x_reverse() +
  scale_color_manual(values = patho_color[names(table(sce$Region))], 
                     na.value = "#d3d3d3") + labs(color = "Pathology")

pdf(file = file.path(plot_savepath, paste0("L1_4", "_patho.pdf")),
    width = 8,
    height = 7)
print(plot_patho)
dev.off()


