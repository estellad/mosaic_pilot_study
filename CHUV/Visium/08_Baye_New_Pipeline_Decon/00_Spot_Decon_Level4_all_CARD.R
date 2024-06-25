# disease = "breast"
# disease = "lung"
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

## Get chromium by indication ------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/" # breast / # lung
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

vis_rawmat <- "counts"

i = 1
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
  
  ## Get Spot gene expr by sample -----------------------------------------
  sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
  
  save_bs_path_sample <- paste0(save_bs_path, foldername, "/", save_names[i], "/")
  
  ## Decon -------------------------------------------------------------------
  # Ref: prep count matrix and cell anno
  cell_types <- data.frame(
    cellID = colnames(chrom),
    cellType = as.factor(chrom@meta.data[["Harmonised_Level4"]]), 
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
  results = results[, sort(colnames(results))]
  write.csv(results, paste0(save_bs_path_sample, save_names[i], "_spot_Level4_decon.csv"))
  results <- read.csv(paste0(save_bs_path_sample, save_names[i], "_spot_Level4_decon.csv"), row.names = 1)
  
  # Visualize spatial pie chart -----------------------------------------
  # colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
  #            "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
  #            "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D") # to cell types in alphabetic order
  
  level4_decon_cols = level4_cellcolors[names(level4_cellcolors) %in% sort(colnames(results))]
  
  # colors = level4_cellcolors
  p1 <- CARD.visualize.pie(
    proportion = results,
    spatial_location = CARD_obj@spatial_location, 
    colors = level4_decon_cols, 
    radius = 9) + theme(legend.position = "right") # + # coord_flip() + 
    # scale_y_reverse() +
    # scale_x_reverse()
  print(p1)

  image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
  image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
  
  # # if need to coord_flip an image
  # image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
  # image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
  
  plot_width <- 20
  plot_height <- plot_width * (image_height/image_width) -2

  pdf(file.path(save_bs_path_sample, paste0(save_names[i], "_spot_level4_decon_bycelltype.pdf")),
      height = plot_height, width = plot_width)
  print(p1)
  dev.off()
  
  
  # Decon label ---------------------------------------------------------
  # 1st dominant cell type
  results$max <- colnames(results)[max.col(results, ties.method="first")]
  
  result_max <- NULL
  for(s in 1:nrow(results)){
    results_s <- results[s, ]
    if(results_s[results_s$max] < 0.2 ){
      results_s_max = "Mix"
    }else{
      results_s_max = results_s$max
    }
    result_max <- c(result_max, results_s_max)
  }
  results2 <- results
  results2$max <- result_max
  
  if(nrow(results2) < ncol(sce)){
    sce$Level4_decon_max <- NULL
    sce$barcode <- colnames(sce)
    results2$barcode <- rownames(results2)
    
    CD <- as.data.frame(colData(sce))
    CD <- CD %>% 
      left_join(results2 %>% select(max, barcode)) %>%
      rename(Level4_decon_max = max) %>%
      mutate(Level4_decon_max = factor(Level4_decon_max, levels = c(sort(unique(results2$max))[sort(unique(results2$max)) != "Mix"], "Mix")))
    rownames(CD) <- CD$barcode; CD$barcode <- NULL
    
    colData(sce) <- as(CD, "DFrame")
  }else{
    sce[["Level4_decon_max"]] <- factor(results2$max, levels = c(sort(unique(results2$max))[sort(unique(results2$max)) != "Mix"], "Mix"))
  }
  
  level4_decon_cols_max <- plyr::mapvalues(levels(sce[["Level4_decon_max"]]), from = level4_cellnames, to = level4_cellcolors)
  
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/08_Baye_New_Pipeline_Decon/helpers.R")
  
  sce_plot <- sce[, !is.na(sce$Level4_decon_max)]
  p2 <- plotSpots_deconmax(sce_plot, 
                           annotate = "Level4_decon_max",
                           pt.size = 5,
                           palette = level4_decon_cols_max) # + 
  # coord_flip() + 
  # scale_y_reverse() +
  # scale_x_reverse()
  
  pdf(file.path(save_bs_path_sample, paste0(save_names[i], "_spot_level4_decon_bycelltype_majorityvote.pdf")),
      height = plot_height, width = plot_width)
  print(p2)
  dev.off()
  
  # Decon label (merge top 2 labels if > 20%) ----------------------------------
  # 2nd dominant cell type
  results2nd <- c()
  for(s in 1:nrow(results)){
    results_s <- results[s, !(names(results[s, ]) %in% results[s, ]$max)]
    results_s_2nd <- colnames(results_s)[max.col(results_s[, 1:(ncol(results)-2)], ties.method="first")]
    
    results2nd <- c(results2nd, results_s_2nd)
  }
  results$max2nd <- results2nd
  
  # final plotting label
  results_label <- c()
  for(s in 1:nrow(results)){
    results_s <- results[s, ]
    if(results_s[results_s$max] >= 0.2 & results_s[results_s$max2nd] >= 0.2){ # if 2 or more than 2 >= 0.2
      test <- sort(c(results_s$max, results_s$max2nd))
      results_label <- c(results_label, paste0(test[1], ' - ', test[2]))
    }else if(all(results_s[, -which(names(results_s) %in% c("max", "max2nd"))] < 0.2)){
      results_label <- c(results_label, "Mix")
    }else{ # if only 1 >= 0.2, by max vote
      results_label <- c(results_label, results_s$max)
    }
  }
  results$results_label <- results_label
  
  if(nrow(results) < ncol(sce)){
    sce$Level4_decon_mixtypes <- NULL
    sce$barcode <- colnames(sce)
    results2 <- results
    results2$barcode <- rownames(results2)
    
    CD <- as.data.frame(colData(sce))
    CD <- CD %>% 
      left_join(results2 %>% select(results_label, barcode)) %>%
      rename(Level4_decon_mixtypes = results_label)
    rownames(CD) <- CD$barcode; CD$barcode <- NULL
    
    colData(sce) <- as(CD, "DFrame")
  }else{
    sce[["Level4_decon_mixtypes"]] <- factor(results$results_label) # , levels = c(sort(unique(results_label))[sort(unique(results_label)) != "Mix"], "Mix"))
  }

  # To plot
  
  # Save back the object ------------------------------------------------
  saveRDS(sce, file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
# }


