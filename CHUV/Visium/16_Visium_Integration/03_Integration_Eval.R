library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(SingleCellExperiment)
library(standR)
library(cluster)
library(mclustcomp)


disease = "lung"
## LogNorm ##############################################################
## Before SpotClean -----------------------------------------------------
datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
savepath.filtered = paste0(datapath.filtered, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-logNormnoSpotClean_filtered.rds")
savepath = savepath.filtered
pipeline = 1
seu <- readRDS(paste0(savepath, save_rds_name))
assign(paste0("seu_", pipeline), seu)

## After SpotClean ------------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-logNormpostSpotClean.rds")
savepath = savepath.spotclean
pipeline = 2

seu <- readRDS(paste0(savepath, save_rds_name))
assign(paste0("seu_", pipeline), seu)

## SCT #################################################################
# Before spotclean -----------------------------------------------------
datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
savepath.filtered = paste0(datapath.filtered, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTnoSpotClean_filtered.rds")
savepath = savepath.filtered
pipeline = 3

seu <- readRDS(paste0(savepath, save_rds_name))
assign(paste0("seu_", pipeline), seu)

# After spotclean ------------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
savepath = savepath.spotclean
pipeline = 4

seu <- readRDS(paste0(savepath, save_rds_name))
assign(paste0("seu_", pipeline), seu)


#######################################################################
# Preprocessing -----------------------------------------------------------
seu_1 <- JoinLayers(seu_1) # Active assay Spatial
spe1 <- as.SingleCellExperiment(seu_1)

seu_2 <- JoinLayers(seu_2) # Active assay Spatial
spe2 <- as.SingleCellExperiment(seu_2)

seu_3[['Spatial']] <- NULL # Active assay SCT
spe3 <- as.SingleCellExperiment(seu_3)
assay(spe3, "logcounts") <- as(GetAssayData(seu_3, assay = "SCT", layer = "counts"), "dgCMatrix") # counts() here should take 'Spatial' joint counts, but not needed for downstream

seu_4[['Spatial']] <- NULL # Active assay SCT
spe4 <- as.SingleCellExperiment(seu_4)
assay(spe4, "logcounts") <- as(GetAssayData(seu_4, assay = "SCT", layer = "counts"), "dgCMatrix") # counts() here should take 'Spatial' joint counts, but not needed for downstream


# Convert Seurat to SPE ---------------------------------------------------
object_names <- paste0("seu_", 1:4)
seu_list <- mget(object_names)

spe_list <- lapply(seu_list, as.SingleCellExperiment)


# Subset to healthy vs tumor ----------------------------------------------
spe_list_tu <- lapply(spe_list, function(x){
  x <- x[ , grepl("Tu_L", x$Level4_decon_max)]
})

spe_list_he <- lapply(spe_list, function(x){
  x <- x[ , !grepl("Tu_L", x$Level4_decon_max)]
})


# Plot --------------------------------------------------------------------
colors <- c("#daead3ff", "#6aa84fff", "#ead2dcff", "#a64d79ff")
methods <- c("LogNorm", "LogNorm\nSpotClean", "SCTransform", "SCTransform\nSpotClean")
names(colors) <- methods

# For tumor spots, want sample id to be separate
stat_tu_id <- mutate(bind_rows(lapply(spe_list_tu, function(x) {
  computeClusterEvalStats(x, "sample_id")
})), from = rep(methods, each = 6))

# For healthy spots, want sample id to be overlapped
stat_he_id <- mutate(bind_rows(lapply(spe_list_he, function(x) {
  computeClusterEvalStats(x, "sample_id")
})), from = rep(methods, each = 6))

# For healthy spots, want biology to be still distinct
stat_he_bio <- mutate(bind_rows(lapply(spe_list_he, function(x) {
  computeClusterEvalStats(x, "Level4_decon_max")
})), from = rep(methods, each = 6))


# -------------------------------------------------------------------------
p_tu_id <- ggplot(mutate(mutate(stat_tu_id, from = factor(from, levels = methods)),
                       scores = as.numeric(scores)), aes(from, scores, fill = from)) +
  geom_bar(stat = "identity", col = "black", width = 0.4) +
  facet_wrap(~types, scales = "free_y", ncol = 6) +
  theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(),
                     axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5)) +
  ylab("Scores") + ggtitle("Sample (Tumor spots)")

p_he_id <- ggplot(mutate(mutate(stat_he_id, from = factor(from, levels = methods)),
                       scores = as.numeric(scores)), aes(from, scores, fill = from)) +
  geom_bar(stat = "identity", col = "black", width = 0.4) +
  facet_wrap(~types, scales = "free_y", ncol = 6) +
  theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(),
                     axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5)) +
  ylab("Scores") + ggtitle("Sample (Healthy spots)")

p_he_bio <- ggplot(mutate(mutate(stat_he_bio, from = factor(from, levels = methods)),
                         scores = as.numeric(scores)), aes(from, scores, fill = from)) +
  geom_bar(stat = "identity", col = "black", width = 0.4) +
  facet_wrap(~types, scales = "free_y", ncol = 6) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5)) +
  xlab("Count data") + ylab("Scores") + ggtitle("Cell type deconvolution majority vote (Healthy spots)")


p_tu_id <- p_tu_id + scale_fill_manual(values = colors)
p_he_id <- p_he_id + scale_fill_manual(values = colors)
p_he_bio <- p_he_bio + scale_fill_manual(values = colors)


p_tu_id / p_he_id / p_he_bio



# # -------------------------------------------------------------------------
# stat_patho <- mutate(bind_rows(lapply(spe_list, function(x) {
#   computeClusterEvalStats(x, "Region")
# })), from = rep(methods, each = 6))
# stat_decon <- mutate(bind_rows(lapply(spe_list, function(x) {
#   computeClusterEvalStats(x, "Level4_decon_max")
# })), from = rep(methods, each = 6))
# stat_patient <- mutate(bind_rows(lapply(spe_list, function(x) {
#   computeClusterEvalStats(x, "sample_id")
# })), from = rep(methods, each = 6))
# 
# 
# # -------------------------------------------------------------------------
# p_patho <- ggplot(mutate(mutate(stat_patho, from = factor(from, levels = methods)), 
#                        scores = as.numeric(scores)), aes(from, scores, fill = from)) + 
#   geom_bar(stat = "identity", col = "black", width = 0.4) + 
#   facet_wrap(~types, scales = "free_y", ncol = 6) + 
#   theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(),
#                      axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5)) + 
#   ylab("Scores") + ggtitle("Pathology label")
# 
# p_decon <- ggplot(mutate(mutate(stat_decon, from = factor(from, levels = methods)), 
#                        scores = as.numeric(scores)), aes(from, scores, fill = from)) + 
#   geom_bar(stat = "identity", col = "black", width = 0.4) + 
#   facet_wrap(~types, scales = "free_y", ncol = 6) + 
#   theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(),
#                      axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5)) + 
#   ylab("Scores") + ggtitle("Cell type deconvolution majority vote")
# 
# p_patient <- ggplot(mutate(mutate(stat_patient, from = factor(from, levels = methods)), 
#                          scores = as.numeric(scores)), aes(from, scores, fill = from)) + 
#   geom_bar(stat = "identity", col = "black", width = 0.4) + 
#   facet_wrap(~types, scales = "free_y", ncol = 6) + 
#   theme_bw() + theme(legend.position = "none",
#                      axis.text.x = element_text(angle = 30, vjust = 0.7, hjust=0.5)) + 
#   xlab("Count data") + ylab("Scores") + ggtitle("Sample")
# 
# 
# p_patho <- p_patho + scale_fill_manual(values = colors)
# p_decon <- p_decon + scale_fill_manual(values = colors)
# p_patient <- p_patient + scale_fill_manual(values = colors)
# 
# 
# p_patho / p_decon / p_patient

