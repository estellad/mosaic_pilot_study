library(dplyr)
library(tibble)
library(tidyverse)
# Decon not so reliable, focus on marker
###########################################################################
# Decon result of immune cell types fractions in these two regions --------
sce <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/B3_2_baye_clustered.rds")

# decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/Breast/B3_2/"
decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/C2L"
# decon <- read.csv(file.path(decon_path, "B3_2_spot_level1_5_decon.csv")) %>%
# decon <- read.csv(file.path(decon_path, "B3_2_spot_level4_RCTD.csv")) %>% # RCTD
# decon <- read.csv(file.path(decon_path, "B3_2_spot_Level4_decon.csv")) %>% # CARD
decon <- read.csv(file.path(decon_path, "B3_2.csv")) %>% dplyr::rename(Fibroblast = Fibroblast_B3) %>% # C2L
  dplyr::rename(Barcode = X)
sce$x_coord <- spatialCoords(sce)[, 1]
sce$y_coord <- spatialCoords(sce)[, 2]
CD <- as.data.frame(colData(sce)) %>%
  select(Barcode, spatial.cluster, spatial.cluster_merge_final_consensus_TME, x_coord, y_coord)

# # -------------------------------------------------------------------------
# decon_sub <- decon %>%
#   left_join(CD)
# 
# library(ggplot2)
# p1 <- ggplot(decon_sub, aes(x_coord, y_coord)) + 
#   geom_point(aes(colour = Tu_B3_CYP4F8/Fibroblast_B3)) + 
#   theme_bw() + 
#   scale_color_viridis_c(option = "magma")
# 
# p2 <- ggplot(decon_sub, aes(x_coord, y_coord)) + 
#   geom_point(aes(colour = Fibroblast_B3)) + 
#   theme_bw() + 
#   scale_color_viridis_c(option = "magma")
# 
# p3 <- ggplot(decon_sub, aes(x_coord, y_coord)) + 
#   geom_point(aes(colour = Tu_B3_CYP4F8)) + 
#   theme_bw() + 
#   scale_color_viridis_c(option = "magma")
# 
# p1 | p2 | p3

# # Consensus -------------------------------------------------------------
decon_sub <- decon %>%
  left_join(CD) %>%
  filter(spatial.cluster_merge_final_consensus_TME %in% c("1_5_7_9_Consensus", "1_5_7_9_TME", "4_14_Consensus", "4_14_TME")) %>%
  column_to_rownames("Barcode")

decon_sub_long <- decon_sub %>%
  gather(CellType, Fraction, -spatial.cluster_merge_final_consensus_TME)
decon_sub_long$spatial.cluster_merge_final_consensus_TME <- factor(decon_sub_long$spatial.cluster_merge_final_consensus_TME,
                                                                   levels = c("4_14_Consensus", "1_5_7_9_Consensus", "4_14_TME", "1_5_7_9_TME"))

decon_sub_long_fibro_macro <- decon_sub_long %>%
  filter(CellType %in% c("Fibroblast_B3", "Macrophage", "Endothelia_vascular", "Pericyte", "Tu_B3_CYP4F8")) %>%
  mutate(Area = spatial.cluster_merge_final_consensus_TME) %>%
  mutate(Area = case_when(Area == "4_14_Consensus" ~ "Area A Tumor",
                          Area == "1_5_7_9_Consensus" ~ "Area B Tumor",
                          Area == "4_14_TME" ~ "Area A TME",
                          Area == "1_5_7_9_TME" ~ "Area B TME")) %>%
  filter(Area %in% c("Area A Tumor", "Area B Tumor"))

decon_sub_long_fibro_macro$Area <- factor(decon_sub_long_fibro_macro$Area,
                                          levels = c("Area A Tumor", "Area B Tumor")) #, "Area A TME", "Area B TME"))

decon_sub_long_fibro_macro$Fraction <- as.numeric(decon_sub_long_fibro_macro$Fraction)
# decon_sub_long
p <- ggplot(decon_sub_long_fibro_macro, aes(x=CellType, y=Fraction, fill= Area)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("#0000ff", "#ff00ff" )) + # , "#00008B", "#AA336A")) +
  theme(axis.text.x = element_text(size = 12.5, angle = 45, hjust=1),
        panel.spacing=unit(1.5,"lines"),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 13.5, face = "bold"),
        strip.background=element_rect(fill="#DEDEDE"))


# # Not restrict to Consensus -------------------------------------------------------------
decon_sub <- decon %>%
  left_join(CD) %>%
  filter(spatial.cluster %in% c("1", "5", "9", "14")) %>%
  column_to_rownames("Barcode") %>%
  mutate(spatial.cluster_merge_final_consensus_TME = ifelse(spatial.cluster == "14", spatial.cluster, "1_5_9"))

decon_sub_long <- decon_sub %>%
  gather(CellType, Fraction, -spatial.cluster_merge_final_consensus_TME)
decon_sub_long$spatial.cluster_merge_final_consensus_TME <- factor(decon_sub_long$spatial.cluster_merge_final_consensus_TME,
                                                                   levels = c("1_5_9", "14"))

decon_sub_long_fibro_macro <- decon_sub_long %>%
  filter(CellType %in% c("Fibroblast", "Macrophage", "Endothelia_vascular", "Pericyte", "Tu_B3_CYP4F8")) %>%
  mutate(Area = spatial.cluster_merge_final_consensus_TME) %>%
  mutate(Area = case_when(Area == "14" ~ "Area A",
                          Area == "1_5_9" ~ "Area B")) %>%
  filter(Area %in% c("Area A", "Area B"))

decon_sub_long_fibro_macro$Area <- factor(decon_sub_long_fibro_macro$Area,
                                          levels = c("Area A", "Area B")) #, "Area A TME", "Area B TME"))
decon_sub_long_fibro_macro$Fraction <- as.numeric(decon_sub_long_fibro_macro$Fraction)

# decon_sub_long
p <- ggplot(decon_sub_long_fibro_macro, aes(x=CellType, y=Fraction, fill= Area)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("#0000ff", "#ff00ff" )) + # , "#00008B", "#AA336A")) +
  theme(axis.text.x = element_text(size = 12.5, angle = 45, hjust=1),
        panel.spacing=unit(1.5,"lines"),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 13.5, face = "bold"),
        strip.background=element_rect(fill="#DEDEDE"))



# T-test ------------------------------------------------------------------
library(rstatix)
stat.test_Consensus <- decon_sub_long_fibro_macro %>%
  # filter(spatial.cluster_merge_final_consensus_TME %in% c("4_14_Consensus", "1_5_7_9_Consensus")) %>%
  filter(spatial.cluster_merge_final_consensus_TME %in% c("14", "1_5_9")) %>%
  group_by(CellType) %>%
  t_test(Fraction ~ spatial.cluster_merge_final_consensus_TME) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test_Consensus

# CellType      .y.      group1         group2        n1    n2 statistic    df       p   p.adj p.adj.signif
# <chr>         <chr>    <chr>          <chr>      <int> <int>     <dbl> <dbl>   <dbl>   <dbl> <chr>       
#   1 Fibroblast_B3 Fraction 4_14_Consensus 1_5_7_9_C…   197   674     -2.31  324. 2.14e-2 4.28e-2 *           
#   2 Macrophage    Fraction 4_14_Consensus 1_5_7_9_C…   197   674      4.65  257. 5.43e-6 1.09e-5 **** 

# CellType            .y.      group1         group2               n1    n2 statistic    df        p   p.adj p.adj.signif
# <chr>               <chr>    <chr>          <chr>             <int> <int>     <dbl> <dbl>    <dbl>   <dbl> <chr>       
#   1 Endothelia_vascular Fraction 4_14_Consensus 1_5_7_9_Consensus   197   674    -0.517  281.  6.06e-1 1   e+0 ns          
#   2 Fibroblast_B3       Fraction 4_14_Consensus 1_5_7_9_Consensus   197   674    -2.31   324.  2.14e-2 8.56e-2 ns          
#   3 Macrophage          Fraction 4_14_Consensus 1_5_7_9_Consensus   197   674     4.65   257.  5.43e-6 2.17e-5 ****        
#   4 Pericyte            Fraction 4_14_Consensus 1_5_7_9_Consensus   197   674    -6.12   377.  2.41e-9 9.64e-9 ****    

stat.test_TME <- decon_sub_long_fibro_macro %>%
  filter(spatial.cluster_merge_final_consensus_TME %in% c("4_14_TME", "1_5_7_9_TME")) %>%
  group_by(CellType) %>%
  t_test(Fraction ~ spatial.cluster_merge_final_consensus_TME) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test_TME

# CellType      .y.      group1   group2         n1    n2 statistic    df       p   p.adj p.adj.signif
# <chr>         <chr>    <chr>    <chr>       <int> <int>     <dbl> <dbl>   <dbl>   <dbl> <chr>       
#   1 Fibroblast_B3 Fraction 4_14_TME 1_5_7_9_TME   234   431     -3.82  499. 0.00015 0.0003  ***         
#   2 Macrophage    Fraction 4_14_TME 1_5_7_9_TME   234   431      3.38  442. 0.00079 0.00158 **     

# CellType            .y.      group1   group2         n1    n2 statistic    df           p      p.adj p.adj.signif
# <chr>               <chr>    <chr>    <chr>       <int> <int>     <dbl> <dbl>       <dbl>      <dbl> <chr>       
#   1 Endothelia_vascular Fraction 4_14_TME 1_5_7_9_TME   234   431     -1.89  467. 0.0588      0.235      ns          
#   2 Fibroblast_B3       Fraction 4_14_TME 1_5_7_9_TME   234   431     -3.82  499. 0.00015     0.0006     ***         
#   3 Macrophage          Fraction 4_14_TME 1_5_7_9_TME   234   431      3.38  442. 0.00079     0.00316    **          
#   4 Pericyte            Fraction 4_14_TME 1_5_7_9_TME   234   431     -5.05  516. 0.000000628 0.00000251 ****  

test <- decon_sub_long %>%
  filter(CellType != "Tu_B3_CYP4F8") %>%
  group_by(CellType) %>%
  summarise(mean_val = mean(Fraction))

plot(density(test$mean_val))
abline(v = mean(test$mean_val))
summary(test$mean_val)

test_ <- test %>%
  filter(mean_val > mean(test$mean_val))

# We investigated the average deconvolution fraction for healthy cell types in Area A and Area B,
# for which four cell types stand out with average deconvolution fractions that are above the mean value across all cell types. 





