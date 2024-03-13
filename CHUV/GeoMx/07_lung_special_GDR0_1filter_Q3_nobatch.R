disease = "lung"
# disease = "breast"
# disease = "dlbcl"

save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_test_norm/")
spe <- readRDS(file.path(save_path, paste0(disease, ".rds")))

spe <- spe[, spe$gene_detection_rate >= 0.1] # 122 -> 47

## Upperquartile (Q3)
spe_upperquatile <- geomxNorm(spe, method = "upperquartile", log = FALSE)
assays(spe, withDimnames=FALSE)$upperquartile <- assay(spe_upperquatile, "logcounts")

assay = 7
which.assay = "upperquartile"

set.seed(100)
spe <- scater::runPCA(spe, assay.type = which.assay)
set.seed(100)
spe <- scater::runUMAP(spe, dimred = "PCA")

# PCA before batch
p1 <- plotDimRed(spe, type = "PCA", annotate = "section_id", pt.size = 1.5) + 
  theme(legend.position = "right",
        legend.title=element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) + 
  ggtitle(str_to_title(disease))

p2 <- plotDR(spe, dimred = "PCA", col = slide_name) + xlab("PC_1") + ylab("PC_2") + 
  theme(legend.position = "right", legend.title=element_blank(), plot.title = element_text(hjust = 0.5)) + ggtitle("Slide")
p3 <- plotDR(spe, dimred = "PCA", col = cell_fraction) + xlab("PC_1") + ylab("PC_2") + 
  theme(legend.position = "right", legend.title=element_blank(), plot.title = element_text(hjust = 0.5)) + ggtitle("Cell fraction")

p1 | p2 | p3

# UMAP before batch
plotDR(spe, dimred = "UMAP", col = section_id) + xlab("UMAP_1") + ylab("UMAP_2") | 
  plotDR(spe, dimred = "UMAP", col = slide_name) + xlab("UMAP_1") + ylab("UMAP_2") | 
  plotDR(spe, dimred = "UMAP", col = cell_fraction) + xlab("UMAP_1") + ylab("UMAP_2")


plotDR(spe, dimred = "UMAP", col = gene_detection_rate) + xlab("UMAP_1") + ylab("UMAP_2") 

##############################################################
#                          Decon                             #
##############################################################

# NegProbe name
negProbeName <- rownames(spe)[which(grepl("Neg", rownames(spe)))]

# ## Normed
spd_normed <- prepareSpatialDecon(spe, assay2use = which.assay, negProbeName = negProbeName)

# ## Normed Batched
# # force the normed, batch corrected count matrix to be not a df, but just matrix
# assay(spe_ruv, "logcounts") <- as.matrix(assay(spe_ruv, "logcounts"))
# spd_normed_batched <- prepareSpatialDecon(spe_ruv, assay2use = "logcounts", negProbeName = negProbeName)

# Signature matrix
new_ref_matrix <- as.matrix(read.csv(paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/Signature_matrices_CHUV_lung_breast.csv"), 
                                     row.names = 1))

library(SpatialDecon)
res_norm <- spatialdecon(norm = spd_normed$normCount,
                         bg = spd_normed$backGround,
                         X = new_ref_matrix,
                         # cellmerges = safeTME.matches,
                         align_genes = TRUE)

subset_prop <- res_norm$prop_of_all


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
  rename(CellType = column,
         Fraction = value)

# Plot
p <- ggplot(gathered_df, aes(x=CellType, y=Fraction)) +
  geom_boxplot() +
  geom_jitter(aes(color=section_id), size=1, alpha=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~cell_fraction, ncol = 5) + 
  # ggtitle("Quantile Normed Expr Decon Result")
  # ggtitle("Quantile Normed RUV4 Batch Corrected Expr Decon Result")
  ggtitle("GDR filtered Q3 Normed Expr Decon Result")










