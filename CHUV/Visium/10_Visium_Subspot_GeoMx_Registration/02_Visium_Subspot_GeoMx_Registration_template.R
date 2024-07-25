
# # Left name is HE, right name is GeoMX
# he_to_fluo_matching = {
#   # TODO: No patho Visium annotation "B1-2": "B1-1",
#   "B1-4": "B1-3",
#   "B2-2": "B2-1",
#   "B3-2": "B3-1",
#   "B4-2": "B4-1",
#   "L1-2": "L1-1",
#   # L1-4": None,
#   "L2-2": "L2-1", # Could also be L2-3, don't know why it was replicated while there is only one HE
#   # TODO: No patho Visium annotaiton "L3-2": "L3-3", # /!\ should be L3-1 but fluo didnt work so they re-did it and replaced it w L3-3
#   "L4-2": "L4-1"
# }

# disease = "breast"
# source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
# samples_for_registration <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2") # Visium
# geo_reg_names <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1")


disease = "lung"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
samples_for_registration <- c("L1_2", "L2_2", "L3_2", "L4_2") # Visium
geo_reg_names <- c("L1_1", "L2_1", "L3_3", "L4_1")


i = 3
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
# Enhanced Expr -----------------------------------------------------------
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/", foldername, "/", samples_for_registration[i], "/")
sce_spot <- readRDS(file.path(save_bs_path, paste0(samples_for_registration[i], "_baye_clustered.rds"))) # spot
sce <- readRDS(file.path(save_bs_path, paste0(samples_for_registration[i], "_baye_clustered_all_enhanced_expr.rds"))) # subspot

# Map in barcode to subspot L3_2 --------------
CD_spot <- as.data.frame(colData(sce_spot))
CD_spot$barcode <- rownames(CD_spot)
CD <- as.data.frame(colData(sce))

CD_spot <- CD_spot %>%
  dplyr::rename(spot.row = array_row,
                spot.col = array_col) %>%
  select(spot.row, spot.col, barcode
         #, Region
         )

CD <- CD %>%
  left_join(CD_spot, by = c("spot.row", "spot.col"))
rownames(CD) <- colnames(sce)

colData(sce) <- as(CD, "DFrame")

library(SingleCellExperiment)
assays(sce)
# logcounts(sce)

# markers <- list()
# markers[["T-cell"]] <- c("CD2", "CD3D", "CD3E", "CD3G", "CD7")
# markers[["B-cell"]] <- c("CD38", "CD86", "IGKC", "CD27") # "CD19", 
# 
# all(unlist(markers) %in% rownames(sce))
# 
# expr_subset <- assay(sce, "log1p")[unlist(markers), ]

head(colData(sce))
# DataFrame with 6 rows and 11 columns
# spot.idx subspot.idx  spot.row  spot.col array_row array_col pxl_row_in_fullres pxl_col_in_fullres spatial.cluster            barcode              Region
# <numeric>   <integer> <integer> <integer> <numeric> <numeric>          <numeric>          <numeric>       <numeric>        <character>         <character>
#   subspot_1.1         1           1        47        71   47.3333  71.33333            1613.27           1067.603               1 AACACTTGGCAAGGAA-1 Intratumoral_Stroma
# subspot_2.1         2           1        49        43   49.3333  43.33333            1651.27            765.603              10 AACAGGATTCATAGTT-1    Tumor_Stroma_mix
# subspot_3.1         3           1        51        41   51.3333  41.33333            1688.27            743.603               5 AACAGGTTCACCGAAG-1    Tumor_Stroma_mix
# subspot_4.1         4           1        24         6   24.3333   6.33333            1180.27            365.603               4 AACAGTCAGGCTCCGC-1 Intratumoral_Stroma
# subspot_5.1         5           1        58        26   58.3333  26.33333            1820.27            581.603               1 AACATAGTCTATCTAC-1 Intratumoral_Stroma
# subspot_6.1         6           1        67        55   67.3333  55.33333            1989.27            894.603               2 AACATCTAATGACCGG-1          Tumor_pure

# read in mapped subspot with GeoMx 
library(readxl)
library(dplyr)
final_subspot <- read_xlsx("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Vis_Geo_Match/subspot_aoi_matched.xlsx")

head(final_subspot)
# # A tibble: 6 × 15
#Section_ID   ROI.ID     roi    Spot.array.location Spot.barcode  segment Total.intersection.a…¹ Subspot.area..px. Section_ID_Vis overlap_pct Sample_ID Cell_fraction Cell_fraction_location
#  <chr>      <chr>     <chr>   <chr>               <chr>         <chr>                    <dbl>             <dbl> <chr>                <dbl> <chr>     <chr>         <chr>                 
# 1 B1_1       ROI:0  A_Islet_2 93x35.2             CGTTGTGAGCAA… PanCK+                  1586.              5124. B1_2               0.310   DSP-1001… Malignant     Malignant             
# 2 B1_1       ROI:0  A_Islet_2 93x35.1             CGTTGTGAGCAA… PanCK+                  2801.              5174. B1_2               0.541   DSP-1001… Malignant     Malignant             
# 3 B1_1       ROI:0  A_Islet_2 95x35.6             ATTCCGCAACGC… PanCK+                    50.4             5131  B1_2               0.00982 DSP-1001… Malignant     Malignant             
# 4 B1_1       ROI:0  A_Islet_2 95x35.2             ATTCCGCAACGC… PanCK+                  2669.              5124. B1_2               0.521   DSP-1001… Malignant     Malignant             
# 5 B1_1       ROI:0  A_Islet_2 95x35.1             ATTCCGCAACGC… PanCK+                   159.              5174. B1_2               0.0307  DSP-1001… Malignant     Malignant             
# 6 B1_1       ROI:0  A_Islet_2 90x36.2             GCATTGACTACA… PanCK+                   847.              5124. B1_2               0.165   DSP-1001… Malignant     Malignant 


cc <- strsplit(final_subspot$Spot.array.location, "[x]")
spot_col <- unlist(cc)[2*(1:length(final_subspot$Spot.array.location))-1]
spot_row_subspot <- unlist(cc)[2*(1:length(final_subspot$Spot.array.location))]

cc2 <- strsplit(spot_row_subspot, "[.]")
spot_row <- unlist(cc2)[2*(1:length(spot_row_subspot))-1]
subspot_idx <- unlist(cc2)[2*(1:length(spot_row_subspot))]

final_subspot$spot.col <- as.numeric(spot_col)
final_subspot$spot.row <- as.numeric(spot_row)
final_subspot$subspot.idx <- as.numeric(subspot_idx)

df <- data.frame(colData(sce))

final_subspot_sample <- final_subspot %>%
  filter(Section_ID_Vis == samples_for_registration[i])

# Sanity check --------------------------
# length(unique(final_subspot$Spot.array.location))
# range(final_subspot$spot.col) #  0 125
# range(final_subspot$spot.row) #  6 69

dim(df) # 11994     9
range(df$spot.col) #  31 126
range(df$spot.row) #  4 77

dim(final_subspot_sample) # 3128   18
length(unique(final_subspot_sample$Spot.array.location)) # 1606 # still not unique, need to filter out overlap < 70%
range(final_subspot_sample$spot.col) #  37 96   # within range of df, good
range(final_subspot_sample$spot.row) #  6 69    # within range of df, good

# -----------------------------------------------------

final_subspot_sample_filtered <- final_subspot_sample %>%
  filter(overlap_pct >= 0.7) # for now to ensure uniqueness, later set to 70%
# dim(final_subspot_sample_filtered) # 1075   18
# length(unique(final_subspot_sample_filtered$Spot.array.location)) # 1075, so no more subspot assigned to two diff cell fractions

CD <- df %>%
  left_join(final_subspot_sample_filtered, by = c("spot.col", "spot.row", "subspot.idx")) # Okay CD = df number of rows/subspots

rownames(CD) <- rownames(df)

# # Merge in spot -------------------------------------------------------
# sce_spot <- readRDS(file.path(save_bs_path, "B1_4_baye_clustered.rds")) # 18405  1999
# df_spot <- data.frame(colData(sce_spot))
# range(df_spot$array_col) # 31 126
# range(df_spot$array_row) # 4 77
# 
# df_spot <- df_spot %>%
#   rename(spot.col = array_col,
#          spot.row = array_row) %>%
#   select(spot.col, spot.row, barcode)
# # -----------------------------------------------------------------------
# 
# CD <- CD %>%
#   left_join(df_spot, by = c("spot.col", "spot.row"))
# 
# 
# # Map in Pathology annotation --------------------------
# vis_anno <- read_excel("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Visium_combined_annotations_CS.xlsx")
# vis_anno_sample <- vis_anno %>%
#   filter(Section == "B1_4") %>%
#   select(Barcode, Region) %>%
#   rename(barcode = Barcode) # 1999   17
# 
# # -------------------------------------------------------------------------
# CD <- CD %>%
#   left_join(vis_anno_sample, by = "barcode")

# Adjust subspot to same color ----------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n) # "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"

colData(sce) <- as(CD, "DFrame")
sce$Cell_fraction_color <- case_when(sce$Cell_fraction == "Macro" ~ "#A3A500",
                                     sce$Cell_fraction == "Malignant" ~ "#F8766D",
                                     sce$Cell_fraction == "Other" ~ "#00BF7D",
                                     sce$Cell_fraction == "PanCK-" ~ "#00B0F6",
                                     sce$Cell_fraction == "T_cells" ~ "#E76BF3",
                                     is.na(sce$Cell_fraction) ~ "#757B82")

names(cols) <- c("Malignant", "Macro", "Other" ,"PanCK-", "T_cells")

library(BayesSpace)
library(patchwork)
library(ggspavis)
dim(sce)
sce <- sce[, !is.na(sce$barcode)]

plotSpots(sce, in_tissue = NULL, annotate = "Cell_fraction",
          x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres") + 
  scale_color_manual(values = c(cols[names(table(sce$Cell_fraction))], "#757B82"))
clusterPlot(sce, label = "Region")
table(sce$Region, sce$Cell_fraction)
table(sce$Cell_fraction)


# Save location -----------------------------------------------------------
regis_savepath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration/"
saveRDS(sce, file.path(paste0(regis_savepath, "/Visium_mapped_subspot_obj"),
                       paste0(samples_for_registration[i], "_mapped_subspot.rds")))

geomx_mapped_subspot_LW <- CD %>%
  filter(!is.na(segment))
write.csv(geomx_mapped_subspot_LW, 
          file.path(paste0(regis_savepath, "/GeoMx_mapped_subspot_LW"), 
                    paste0(geo_reg_names[i], "_mapped_LW.csv"))) # column Sample_ID to be filtered out


# Plot --------------------------------------------------------------------
for(i in 1:length(samples_for_registration)){
  sce <- readRDS(file.path(paste0(regis_savepath, "/Visium_mapped_subspot_obj"),
                           paste0(samples_for_registration[i], "_mapped_subspot.rds")))
  p_i <- plotSpots(sce, in_tissue = NULL, annotate = "Cell_fraction",
            x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres") + 
    scale_color_manual(values = c(cols[names(table(sce$Cell_fraction))], "#757B82"))
  
  assign(paste0("p", i), p_i)
}























