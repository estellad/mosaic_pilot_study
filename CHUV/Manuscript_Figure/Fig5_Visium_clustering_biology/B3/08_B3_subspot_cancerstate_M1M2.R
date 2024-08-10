library(ggspavis)
library(dplyr)
library(tidyr)
library(tidyverse)
library(BayesSpace)
library(scater)
library(readxl)
library(janitor)

## M1 M2 aggregated

## T cells

## Cancer state


datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/"

# noSpotClean + logNorm + BayesSpace
sce <- readRDS(file.path(datapath, "B3_2_baye_clustered.rds")) # noSpotClean + logNorm + BayesSpace
sce_enhanced <- readRDS(file.path(datapath, "B3_2_baye_clustered_all_enhanced_expr.rds"))

sce$pxl_col_in_fullres <- NULL
sce$pxl_row_in_fullres <- NULL

# Subspot -----------------------------------------------------------------
# set.seed(100)
# sce_enhanced <- runUMAP(sce_enhanced, assay.type = "log1p")
# set.seed(100)
# sce_enhanced <- runTSNE(sce_enhanced, assay.type = "log1p")
# 
# sce_enhanced$spatial.cluster <- as.factor(sce_enhanced$spatial.cluster)
# saveRDS(sce_enhanced, file.path(datapath, "B3_2_baye_clustered_all_enhanced_expr.rds"))
plotSpots(sce_enhanced, annotate = "spatial.cluster", text_by = "spatial.cluster", 
          x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", in_tissue = NULL,
          y_reverse = FALSE, pt.size = 0.3)

plotDimRed(sce_enhanced, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster") 

# # Merge in pathology 
# CD_enhanced <- as.data.frame(colData(sce_enhanced))
# CD <- as.data.frame(colData(sce))
# 
# CD_enhanced
# 
# table(sce_enhanced$, sce_enhanced$spatial.cluster)

# -------------------------------------------------------------------------
disease = "breast"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
set.seed(100)
sce <- runUMAP(sce, assay.type = "logcounts")
set.seed(100)
sce <- runTSNE(sce, assay.type = "logcounts")

# Visualize immune marker spatially in 159 vs 14 and on UMAP ----------------

# Ind markers -------------------------------------------------------------
marker = "SDC4"  # "CD163"  #  "AKT1"  # "ERBB2" # "TP53" # "MKI67"  # "GSTP1" # "KRAS" # "CD86" # "TGFA"  #"CD68" #"CD19" # "KLRK1"  

# # HER2 positive genes 
# ERBB2, GRB7, TOP2A, PIK3CA, PTEN, and AKT1

# ER markers
# ESR1

# PR markers
# PGR, PGRMC1

plotSpots(sce, annotate = marker, y_reverse = FALSE, pt.size = 1, assay = "log1p") |
  plotSpots(sce_enhanced, annotate = marker, y_reverse = FALSE, pt.size = 0.3, 
            assay = "log1p", in_tissue = NULL,
            x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres")

plotDimRed(sce, type = "UMAP", annotate = marker, assay = "log1p") |
  plotDimRed(sce_enhanced, type = "UMAP", annotate = marker, assay = "log1p") 



## For aggregated markers ----------------------------------------------
# # Itai Yanai paper: 
# * Table S3 - Gene modules for diff cancer cell states
# * Table S4 - Regulon (transcription factor) for diff cancer cell states
# * Table S6 - Differentially expressed M1 and M2

M1_M2 <- readxl::read_excel("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Itai_Yanai_PanCancer_Supp_Table.xlsx", sheet = "TableS6", col_names = FALSE) %>%
  filter(!row_number() %in% c(1, 2, 3)) %>%
  `colnames<-`(c("M1", "M2"))

M1 <- na.omit(M1_M2$M1)
M2 <- na.omit(M1_M2$M2)

cancellstate <- readxl::read_excel("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Itai_Yanai_PanCancer_Supp_Table.xlsx", sheet = "TableS3", col_names = FALSE) %>%
  filter(!row_number() %in% c(1, 2)) %>%
  row_to_names(row_number = 1)

cancer_state <- colnames(cancellstate)

one_cancer_state_markers <- function(state = "Cycle", itaitable = cancellstate){
  string <- na.omit(itaitable[[state]])
  markers <- string[string %in% rownames(sce)]
  return(markers)
}

# TODO: here plot all cancer states spatial distribution

# Cycle <- na.omit(cancellstate$Cycle)
# Stress <- na.omit(cancellstate$Stress)
# Interferon <- na.omit(cancellstate$Interferon)
# Hypoxia <- na.omit(cancellstate$Hypoxia)
# Oxphos <- na.omit(cancellstate$Oxphos)
# Metal <- na.omit(cancellstate$Metal)
# Mesenchymal <- na.omit(cancellstate$Mesenchymal)
# pEMT <- na.omit(cancellstate$pEMT)
# Alveolar <- na.omit(cancellstate$Alveolar)
# Basal <- na.omit(cancellstate$Basal)
# Squamous <- na.omit(cancellstate$Squamous)
# Glandular <- na.omit(cancellstate$Glandular)
# Ciliated <- na.omit(cancellstate$Ciliated)
# AC <- na.omit(cancellstate$AC)
# OPC <- na.omit(cancellstate$OPC)
# NPC <- na.omit(cancellstate$NPC)

one_cancer_state_markers(state = cancer_state[i], itaitable = cancellstate)

# -----------------------------------------------------------------------
# marker <- M1[M1 %in% rownames(sce)]
# marker <- M2[M2 %in% rownames(sce)]

# marker <- Cycle[Cycle %in% rownames(sce)] # both
# marker <- Stress[Stress %in% rownames(sce)] # 14 a bit more
# marker <- Interferon[Interferon %in% rownames(sce)] # neither
# marker <- Hypoxia[Hypoxia %in% rownames(sce)] # both stroma
# marker <- Oxphos[Oxphos %in% rownames(sce)] # tiny bit both at the tip
# marker <- Metal[Metal %in% rownames(sce)] # 14 is more metal 
# marker <- Mesenchymal[Mesenchymal %in% rownames(sce)] # neither
# marker <- pEMT[pEMT %in% rownames(sce)] # both
# marker <- Alveolar[Alveolar %in% rownames(sce)] # 14 is more 
# marker <- Basal[Basal %in% rownames(sce)] # stroma region of both
# marker <- Squamous[Squamous %in% rownames(sce)] # both
# marker <- Glandular[Glandular %in% rownames(sce)] # both
# marker <- Ciliated[Ciliated %in% rownames(sce)]
# marker <- AC[AC %in% rownames(sce)]
# marker <- OPC[OPC %in% rownames(sce)] # neither
# marker <- NPC[NPC %in% rownames(sce)] # both stroma

# ----------------------------
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/15_Visium_Breast_Volcano/Tumor_Subtype_DE/plotSpots_helper.R")
marker <- M2[M2 %in% rownames(sce)]# rownames(sce)[grepl("TNFA", rownames(sce))]
marker <- M1[M1 %in% rownames(sce)] 

plotSpots_new(sce, annotate = marker, y_reverse = FALSE, pt.size = 1, assay = "log1p") + ggtitle("") |
  plotSpots_new(sce_enhanced, annotate = marker, y_reverse = FALSE, pt.size = 0.3, 
                assay = "log1p", in_tissue = NULL,
                x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres")

plotDimRed_new(sce, type = "UMAP", annotate = marker, assay = "log1p") |
  plotDimRed_new(sce_enhanced, type = "UMAP", annotate = marker, assay = "log1p") 

# Macrophage M1 
rownames(sce)[grepl("TNFA", rownames(sce))]
# "TNFAIP8L2" "TNFAIP6"   "TNFAIP8"   "TNFAIP3"   "TNFAIP2"   "TNFAIP8L3" "TNFAIP1"   "TNFAIP8L1"

# PanCancer paper pro-inflammatory: TNF, SPP1, ISG15

# Macrophage M2
rownames(sce)[grepl("TGFB", rownames(sce))]
# "TGFBR3"   "TGFB2"    "TGFBRAP1" "TGFBR2"   "TGFBI"    "TGFBR1"   "TGFB3"    "TGFB1I1"  "TGFBR3L"  "TGFB1" 

# PanCancer paper anti-inflammatory: HLA-DRA, C1QA, CD163

# Separate to sub SCE ------------------------------------------------------
# plot library size
plotSpots(sce, annotate = "sum", y_reverse = FALSE, pt.size = 1)
sce_159 <- sce[, sce$spatial.cluster %in% c("1", "5", "9")] # 17883  1036
avg_libsize_159 <- round(sum(sce_159$sum)/ncol(sce_159), 2) # 11969.74

sce_14 <- sce[, sce$spatial.cluster == "14"] # 17883   225
avg_libsize_14 <- round(sum(sce_14$sum)/ncol(sce_14), 2)  # 7951.14


############################## Submited a job #################################
# # Re-cluster 159 ---------------------------------------------------------
# set.seed(42)
# sce_159 <- spatialPreprocess(sce_159, n.PCs = 50, log.normalize = TRUE)
# set.seed(100)
# sce_159 <- runUMAP(sce_159, assay.type = "logcounts")
# set.seed(100)
# sce_159 <- runTSNE(sce_159, assay.type = "logcounts")
# 
# set.seed(123)
# sce_159 = spatialCluster(sce_159, use.dimred = "PCA", q = 2, nrep = 50000, d = 50, burn.in = 1000)

## Re-cluster 2 clusters
sce_159reclus <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered_159reclus_2clus.rds")

## Re-cluster 3 clusters
sce_159reclus <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered_159reclus_3clus.rds")

sce_159reclus$spatial.cluster <- as.factor(sce_159reclus$spatial.cluster)
plotSpots(sce_159reclus, annotate = "spatial.cluster", y_reverse = FALSE, pt.size = 1, assay = "log1p") |
  plotDimRed(sce_159reclus, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster") |
  plotDimRed(sce_159reclus, type = "TSNE", annotate = "spatial.cluster", text_by = "spatial.cluster")


## Previous cluster 1 5 9 ------------------
sce_159$spatial.cluster <- as.factor(sce_159$spatial.cluster)
plotSpots(sce_159, annotate = "spatial.cluster", y_reverse = FALSE, pt.size = 1, assay = "log1p") |
  plotDimRed(sce_159, type = "UMAP", annotate = "spatial.cluster", text_by = "spatial.cluster") |
  plotDimRed(sce_159, type = "TSNE", annotate = "spatial.cluster", text_by = "spatial.cluster") 

## Merge in reclus 2 clusters label back to original sce object "spatial.cluster"
CD_159reclus <- as.data.frame(colData(sce_159reclus)) %>%
  select(Barcode, spatial.cluster) %>%
  dplyr::rename(spatial.cluster159reclus = spatial.cluster)

CD <- data.frame(colData(sce)) %>%
  left_join(CD_159reclus, by = "Barcode") %>%
  mutate(spatial.cluster159reclus = ifelse(is.na(spatial.cluster159reclus), spatial.cluster, paste0("159_", spatial.cluster159reclus)))

rownames(CD) <- CD$Barcode

colData(sce) <- as(CD, "DFrame")


# Decon not so reliable, focus on marker
# ###########################################################################
# # Decon result of immune cell types fractions in these two regions --------
# decon <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_spot_level1_5_decon.csv") %>%
# # decon <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_spot_Level4_decon.csv") %>%
#   rename(Barcode = X)
# CD <- as.data.frame(colData(sce)) %>%
#   select(Barcode, spatial.cluster)
# 
# decon_sub <- decon %>%
#   left_join(CD) %>%
#   filter(spatial.cluster %in% c("1", "5", "9", "14")) %>%
#   mutate(spatial.cluster2 = ifelse(spatial.cluster == "14", "14", "1_5_9")) %>%
#   select(-spatial.cluster) %>%
#   column_to_rownames("Barcode")
# 
# decon_sub_long <- decon_sub %>%
#   gather(CellType, Fraction, -spatial.cluster2)
# 
# # decon_sub_long
# 
# p <- ggplot(decon_sub_long, aes(x=CellType, y=Fraction, fill= spatial.cluster2)) +
#   geom_boxplot(outlier.size = 0.1) + 
#   theme_bw() +
#   scale_fill_manual(values = c("#0000ff", "#ff00ff")) + 
#   theme(axis.text.x = element_text(size = 12.5, angle = 90, vjust = 0.5, hjust=1),
#         panel.spacing=unit(1.5,"lines"),
#         panel.grid = element_blank(), 
#         strip.text.x = element_text(size = 13.5, face = "bold"), 
#         strip.background=element_rect(fill="#DEDEDE"))

saveRDS(sce, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered.rds")
sce <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered.rds")
# class: SpatialExperiment 
# dim: 17883 2156 
# metadata(1): BayesSpace.data
# assays(3): counts logcounts log1p
# rownames(17883): SAMD11 NOC2L ... MT-ND6 MT-CYB
# rowData names(8): symbol gene_id ... low_abungene is.HVG
# colnames(2156): AACACTTGGCAAGGAA-1 AACAGGATTCATAGTT-1 ... TGTTGGATGGACTTCT-1 TGTTGGCCAGACCTAC-1
# colData names(22): in_tissue array_row ... Level4_decon_max Level4_decon_mixtypes
# reducedDimNames(3): PCA UMAP TSNE
# mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

sce$Tu_Consensus <- ifelse(sce$Level4_decon_max == "Tu_B3_CYP4F8" & sce$Region == "Tumor_pure", TRUE, FALSE)
plotSpotQC(sce, annotate = "Tu_Consensus", type = "spot", y_reverse = FALSE)

sce$Tu_consensus_14 <- ifelse(sce$spatial.cluster_tu_consensus == "14", TRUE, FALSE)
sce$Tu_consensus_159 <- ifelse(sce$spatial.cluster_tu_consensus == "159", TRUE, FALSE)
sce$TME_14 <- ifelse(sce$spatial.cluster_tu_consensus == "14_TME", TRUE, FALSE)
sce$TME_159 <- ifelse(sce$spatial.cluster_tu_consensus == "159_TME", TRUE, FALSE)
plotSpotQC(sce, annotate = "Tu_consensus_14", type = "spots", y_reverse = FALSE)
plotSpotQC(sce, annotate = "Tu_consensus_159", type = "spots", y_reverse = FALSE)
plotSpotQC(sce, annotate = "TME_14", type = "spots", y_reverse = FALSE)
plotSpotQC(sce, annotate = "TME_159", type = "spots", y_reverse = FALSE)

plotDimRed(sce, type = "UMAP", annotate = "Tu_consensus_14") | 
  plotDimRed(sce, type = "UMAP", annotate = "Tu_consensus_159") |
  plotDimRed(sce, type = "UMAP", annotate = "TME_14") |
  plotDimRed(sce, type = "UMAP", annotate = "TME_159")


# M1 M2 ratio -------------------------------------------------------------------------





