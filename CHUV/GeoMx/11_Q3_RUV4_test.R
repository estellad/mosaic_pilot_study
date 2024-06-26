absolute_path_cur <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/"
source(paste0(absolute_path_cur, "env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/GeoMx_init.R"))
read_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_test_norm/")
# save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_Final_Qcd/")
library(ggspavis)
library(patchwork)

disease = "lung"
spe <- readRDS(file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/", paste0(disease, "_spe.rds")))


# Q3 ----------------------------------------------------------------------
assay = 7
which.assay = "upperquartile"

set.seed(100)
spe <- scater::runPCA(spe, assay.type = which.assay)
set.seed(100)
spe <- scater::runUMAP(spe, dimred = "PCA")
set.seed(100)
spe <- scater::runTSNE(spe, dimred = "PCA")


# RUV4 --------------------------------------------------------------------
# Batch correction --------------------------------------------------------
spe <- findNCGs(spe, batch_name = "slide_name", top_n = 300)
metadata(spe) |> names()

## Max biology cluster distinction, and minimize batch distinction
for(i in 10:20){
  spe_ruv <- geomxBatchCorrection(spe, factors = "cell_fraction", n_assay = assay,
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = "logcounts", n_dimension = 4, color = cell_fraction, title = paste0("k = ", i)))
  
}

spe_ruv <- geomxBatchCorrection(spe, factors = "cell_fraction", n_assay = assay,
                                NCGs = metadata(spe)$NCGs, k = 20)
set.seed(100)
spe_ruv <- scater::runPCA(spe_ruv, assay.type = "logcounts")


# Decon on Q3 + RUV4 ------------------------------------------------------
# NegProbe name
negProbeName <- rownames(spe)[which(grepl("Neg", rownames(spe)))]

## Normed Batched
# force the normed, batch corrected count matrix to be not a df, but just matrix
assay(spe_ruv, "logcounts") <- as.matrix(assay(spe_ruv, "logcounts"))
assay(spe_ruv, "q3_batched_linear_scale") <- as.matrix(exp(assay(spe_ruv, "logcounts")) - 1)
spd_normed_batched <- prepareSpatialDecon(spe_ruv, assay2use = "q3_batched_linear_scale", negProbeName = negProbeName)

new_ref_matrix <- readRDS(paste0(absolute_path_cur, "Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5/sigmat_breast_add_lung_healthy.rds"))

# -------------------------------------------------------------------------
library(SpatialDecon)
res_norm_batch <- spatialdecon(norm = spd_normed_batched$normCount,
                               bg = spd_normed_batched$backGround,
                               X = new_ref_matrix,
                               align_genes = TRUE)

subset_prop <- res_norm_batch$prop_of_all

# -------------------------------------------------------------------------
df <- data.frame(t(subset_prop))

df$sample <- rownames(df)

library(tidyr)
gathered_df <- gather(df, key = "column", value = "value", -sample)

CD <- as.data.frame(colData(spe)[, c("sample_id2", "section_id", "cell_fraction")])
CD <- CD %>% rename(sample = sample_id2)

library(dplyr)
gathered_df <- gathered_df %>%
  left_join(CD, by = "sample") %>%
  filter(cell_fraction != "PanCK-") %>%
  rename(CellType = column,
         Fraction = value)

# Plot
p <- ggplot(gathered_df, aes(x=CellType, y=Fraction)) +
  geom_boxplot() +
  geom_jitter(aes(color=section_id), size=1, alpha=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~cell_fraction, ncol = 5) + 
  ylim(0.0, 1.0) + 
  # ggtitle("Quantile Normed Batch Corrected Expr Decon Result - full genes")
  ggtitle("Quantile Normed Batch Corrected Expr Decon Result - Linear scale")
p


