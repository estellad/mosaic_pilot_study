source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")

datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
sigmatpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level4/"
# deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results_pt_specific"
deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific"

# prev_deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results/"

###########################################################################
#                               Deconvolution                             #
###########################################################################
# disease = "dlbcl"
# disease = "breast"
disease = "lung"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/00_GeoMx_Paths.R")
sample_name <- ifelse(disease == "dlbcl", "patient", "section_id")
spe_ruv <- readRDS(file.path(datapath, paste0(disease, "_spe_ruv.rds")))

for(i in 1:nsamples){
  print(i)
  # Subset spe to each pt, get reference matrix -----------------------------
  section_name = section_id[i] # "D1"
  pt_name = substr(section_name, 1, 2)
  spe_ruv_pt <- spe_ruv[, spe_ruv[[sample_name]] == section_name]
  
  sigmat_pt <- readRDS(file.path(sigmatpath, paste0("sigmat_", pt_name,".rds")))
  
  
  # Prep decon --------------------------------------------------------------
  # NegProbe name
  negProbeName <- rownames(spe_ruv_pt)[which(grepl("Neg", rownames(spe_ruv_pt)))]
  
  ## Normed Batched
  # force the normed, batch corrected count matrix to be not a df, but just matrix
  assay(spe_ruv_pt, "logcounts") <- as.matrix(assay(spe_ruv_pt, "logcounts"))
  assay(spe_ruv_pt, "quantile_batched_linear_scale") <- as.matrix(exp(assay(spe_ruv_pt, "logcounts")) - 1)
  spd_normed_batched_pt <- prepareSpatialDecon(spe_ruv_pt, assay2use = "quantile_batched_linear_scale", negProbeName = negProbeName)
  
  
  # Deconvolution on normalized and batch corrected object ------------------
  library(SpatialDecon)
  res_norm_batch_pt <- spatialdecon(norm = spd_normed_batched_pt$normCount,
                                    bg = spd_normed_batched_pt$backGround,
                                    X = sigmat_pt,
                                    align_genes = TRUE)
  subset_prop_pt <- res_norm_batch_pt$prop_of_all
  
  
  # Gather decon df ---------------------------------------------------------
  df <- data.frame(t(subset_prop_pt))
  df$sample <- rownames(df)
  
  library(tidyr)
  gathered_df <- gather(df, key = "column", value = "value", -sample)
  
  CD <- as.data.frame(colData(spe_ruv_pt)[, c("sample_id2", sample_name, "cell_fraction")])
  CD <- CD %>% dplyr::rename(sample = sample_id2)
  
  library(dplyr)
  gathered_df <- gathered_df %>%
    left_join(CD, by = "sample") %>%
    # filter(cell_fraction != "PanCK-") %>%
    rename(CellType = column,
           Fraction = value) %>%
    mutate(cell_fraction = ifelse(cell_fraction == "Macro", "Macrophage", cell_fraction))
  
  write.csv(gathered_df, row.names = FALSE,
            file.path(paste0(deconresultpath, "/archive"), 
                      paste0(disease, "_batched_decon_long_", section_name, ".csv")))
  
}


# Concatenate all pt specific csv -----------------------------------------
filenames <- list.files(path = paste0(deconresultpath, "/archive"), 
                        pattern = decon_file_pattern, 
                        full.names = TRUE)

# Read and row bind all CSV files
combined_df <- filenames %>% 
  map_dfr(~ read_csv(.x))

# # previous not pt specific, check if dimension match (previously no level 4 geo decon)
# prev_combined_df <- read.csv(paste0(prev_deconresultpath, disease, "_batched_decon_long.csv"));
# all(dim(prev_combined_df) == dim(combined_df))

write.csv(combined_df, row.names = FALSE,
          file.path(deconresultpath, paste0(disease, "_batched_decon_long.csv")))

# # Plot ------------------------------------------------------------------
# p <- ggplot(combined_df, aes(x=CellType, y=Fraction)) +
#   geom_boxplot() +
#   geom_jitter(aes(color=sample_name), size=1, alpha=0.9) +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust=1),
#         panel.spacing=unit(1.5,"lines"),
#         panel.grid = element_blank(),
#         strip.text.x = element_text(size = 13.5, face = "bold"),
#         strip.background=element_rect(fill="#DEDEDE")) +
#   facet_wrap(~cell_fraction, ncol = 5) +
#   guides(col = guide_legend(override.aes = list(size = 2))) +
#   labs(color = "Section", x = "", y = "Cell type fraction") +
#   ggtitle("")
# 
# p


