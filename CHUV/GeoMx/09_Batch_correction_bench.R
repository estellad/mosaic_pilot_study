source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")
library(stringr)
library(tidyr)
library(dplyr)
datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"

disease = "lung"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/00_GeoMx_Paths.R")
spe <- readRDS(file.path(datapath, paste0(disease, "_spe.rds")))

assay = 8
which.assay = "quantile"

# Correction RUV4 --------------------------------------------------------
spe <- findNCGs(spe, batch_name = "slide_name", top_n = 300)
metadata(spe) |> names()

## Max biology cluster distinction, and minimize batch distinction
spe_ruv <- geomxBatchCorrection(spe, factors = "cell_fraction", n_assay = assay,
                                NCGs = metadata(spe)$NCGs, k = 9)


## Correction Limma ------------------------------------------------------
spe_lrb <- geomxBatchCorrection(spe,
                                batch = spe$slide_name, method = "Limma",
                                design = model.matrix(~cell_fraction, data = colData(spe)))


spe_list <- list(spe, spe_ruv, spe_lrb)

colors <- c("#05c88dff", "#063764ff", "#ffa64fff")
methods <- c("Raw","RUV4","Limma")
names(colors) <- methods

# plotClusterEvalStats(spe_list = spe_list,
#                      bio_feature_name = "cell_fraction",
#                      batch_feature_name = "slide_name",
#                      data_names = methods,
#                      colors = colors)

library("cluster")
library("mclustcomp")
stat_bio <- mutate(bind_rows(lapply(spe_list, function(x) {
  computeClusterEvalStats(x, "cell_fraction")
})), from = rep(methods, each = 6))
stat_batch <- mutate(bind_rows(lapply(spe_list, function(x) {
  computeClusterEvalStats(x, "slide_name")
})), from = rep(methods, each = 6))


# -------------------------------------------------------------------------
p_bio <- ggplot(mutate(mutate(stat_bio, from = factor(from, levels = methods)), 
                       scores = as.numeric(scores)), aes(from, scores, fill = from)) + 
  geom_bar(stat = "identity", col = "black", width = 0.4) + 
  facet_wrap(~types, scales = "free_y", ncol = 3) + 
  theme_bw() + theme(legend.position = "none") + 
  xlab("Count data") + ylab("Scores") + ggtitle("Biology")

p_batch <- ggplot(mutate(mutate(stat_batch, from = factor(from, levels = methods)), 
                         scores = as.numeric(scores)), aes(from, scores, fill = from)) + 
  geom_bar(stat = "identity", col = "black", width = 0.4) + 
  facet_wrap(~types, scales = "free_y", ncol = 3) + 
  theme_bw() + theme(legend.position = "none") + 
  xlab("Count data") + ylab("Scores") + ggtitle("Batch")

p_bio <- p_bio + scale_fill_manual(values = colors)
p_batch <- p_batch + scale_fill_manual(values = colors)
p_bio + p_batch + patchwork::plot_layout(1, 2)


