library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

# -------------------------------------------------------------------------
"Tu_D1" = c("Tu_D1_LMO2", "Tu_D1_RGS13", "Tu_D1_SMIM14");
"Tu_D2" = "Tu_D2_mito";
"Tu_D3" = c("Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD");
"Tu_D4" = c("Tu_D4_BCL7A", "Tu_D4_PNN"); 
"Tu_D5" = "Tu_D5_CCL22";
"Tu_D6" = "Tu_D6_BCL2"

# Decon levels ------------------------------------------------------------
get_result_max <- function(results){
  # Decon label ---------------------------------------------------------
  # 1st dominant cell type
  results$max <- colnames(results)[max.col(results, ties.method="first")]
  
  result_max <- NULL
  for(s in 1:nrow(results)){
    results_s <- results[s, ]
    if(results_s[results_s$max] < 0.5 ){
      results_s_max = "Mix"
    }else{
      results_s_max = results_s$max
    }
    result_max <- c(result_max, results_s_max)
  }
  
  df <- data.frame(Barcode = rownames(results), max = result_max)
  return(df)
}


max_tuunmerged_all <- NULL
max_tumerged_all <- NULL
for(i in 1:6){
  ## Healthy ------------------------------------------------------------
  healthy <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/DLBCL/DLBCL_", i, 
                             "/DLBCL_", i, "_spot_level1_5_decon.csv"))
  
  healthy_sub <- healthy %>%
    mutate(Stroma = Fibro_Muscle + Vessel) %>%
    select(X, B, Epithelia, Myeloid, Stroma, T_NK) %>%
    mutate(X = paste0("DLBCL_", i, "_", X)) 
  
  ## Tumor ------------------------------------------------------------
  tumor <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/DLBCL/DLBCL_", i, 
                           "/DLBCL_", i, "_spot_Level4_decon.csv"))
  # unmerged
  tumor_sub <- tumor %>%
    select(X, all_of(get(paste0("Tu_D", i)))) %>%
    mutate(X = paste0("DLBCL_", i, "_", X)) 
  
  # merged
  tumor_sub_ <- tumor %>%
    select(X, all_of(get(paste0("Tu_D", i)))) %>%
    column_to_rownames("X") %>%
    mutate(Tu_D = rowSums(.)) %>%
    select(Tu_D) %>%
    rownames_to_column("X") %>%
    mutate(X = paste0("DLBCL_", i, "_", X))
  colnames(tumor_sub_)[2] <- paste0(colnames(tumor_sub_)[2], i)
  
  print(all(rownames(healthy_sub) == rownames(tumor_sub)))
  print(all(rownames(healthy_sub) == rownames(tumor_sub_)))
  
  decon_result_tuunmerged <- cbind(healthy_sub %>% column_to_rownames("X"), tumor_sub %>% select(-X)) 
  decon_result_tumerged <- cbind(healthy_sub %>% column_to_rownames("X"), tumor_sub_ %>% select(-X)) 
  
  max_tuunmerged <- get_result_max(decon_result_tuunmerged)
  max_tumerged <- get_result_max(decon_result_tumerged)
  
  max_tuunmerged_all <- rbind(max_tuunmerged_all, max_tuunmerged)
  max_tumerged_all <- rbind(max_tumerged_all, max_tumerged)
}

library(Seurat)
vis <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_Integrated_Owkin/DLBCL_Post_SpotClean/DLBCL-merge-SCTpostSpotClean.rds")
vis$sample_id <- substr(vis$sample_id, 5, nchar(vis$sample_id))

colnames(vis) <- paste0(vis$sample_id, "_", substr(colnames(vis), 1, nchar(colnames(vis)) - 2))


# Unmerged ---------------------------------------------------------------
# Find intersection barcode and filter both
common_barcode <- intersect(colnames(vis), max_tuunmerged_all$Barcode)
vis <- vis[, colnames(vis) %in% common_barcode]
max_tuunmerged_all <- max_tuunmerged_all[max_tuunmerged_all$Barcode %in% common_barcode, ]


max_tuunmerged_all_sorted <- max_tuunmerged_all[order(match(max_tuunmerged_all$Barcode, colnames(vis))), ]

vis$decon_max <- max_tuunmerged_all_sorted$max

# B    Epithelia          Mix      Myeloid       Stroma         T_NK   Tu_D1_LMO2 Tu_D1_SMIM14   Tu_D2_mito 
# 2         1233         7156           66         1023          231         3355           23         3989 
# Tu_D3_FAM3C  Tu_D4_BCL7A  Tu_D5_CCL22   Tu_D6_BCL2 
# 102          394          166          696 

Idents(vis) <- as.factor(vis$decon_max)


# Merged ------------------------------------------------------------------
# Find intersection barcode and filter both
common_barcode <- intersect(colnames(vis), max_tumerged_all$Barcode)
vis <- vis[, colnames(vis) %in% common_barcode]
max_tumerged_all <- max_tumerged_all[max_tumerged_all$Barcode %in% common_barcode, ]


max_tuunmerged_all_sorted <- max_tumerged_all[order(match(max_tumerged_all$Barcode, colnames(vis))), ]

vis$decon_max <- max_tuunmerged_all_sorted$max

# B Epithelia       Mix   Myeloid    Stroma      T_NK     Tu_D1     Tu_D2     Tu_D3     Tu_D4     Tu_D5     Tu_D6 
# 2      1233      4020        45      1006       230      4093      3989      2401       555       166       696 

Idents(vis) <- as.factor(vis$decon_max)


