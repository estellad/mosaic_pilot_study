# Only breast and lung ----------------------------------------------------
# Mapped spot and AOI IDs (spot)
mapped_geo_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration/GeoMx_mapped_spot_LW"
geo_samples <- c("B1_1", "B1_3", "B2_1", "B3_1", "B4_1", "L1_1", "L2_1", "L3_1", "L4_1")
vis_samples <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2", "L1_2", "L2_2", "L3_2", "L4_2")
disease_list <- c(rep("breast", 5), rep("lung", 4))

mapped_all <- NULL
for(i in 1:length(geo_samples)){
  mapped_i <- read.csv(file.path(mapped_geo_path, paste0(geo_samples[i], "_mapped_LW.csv")), row.names = 1) %>% 
    select(barcode, Section_ID, Sample_ID, Cell_fraction) %>%
    rename(cell_fraction = Cell_fraction) %>%
    mutate(Section_ID_Vis = vis_samples[i])
  
  # Merge in pathology ------------------------------------------------------
  sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease_list[i], "_qcd")
  sce <- readRDS(file.path(sce_path, paste0(vis_samples[i], "_qcd.rds")))
  
  mapped_i <- mapped_i %>%
    left_join(data.frame(barcode = colnames(sce),
                         Region = sce$Region))
  
  mapped_all <- rbind(mapped_all, mapped_i)
}

write.csv(mapped_all, 
          file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration", "spot_mapped_cf_region.csv"), row.names = FALSE)


# Geo decon results (spot)
geo_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results"
geo_decon <- rbind(read.csv(file.path(geo_decon_path, "breast_batched_decon_long.csv")),
                   read.csv(file.path(geo_decon_path, "lung_batched_decon_long.csv"))
                   # read.csv(file.path(geo_decon_path, "dlbcl_batched_decon_long.csv")) %>% rename(section_id = patient)
)
# geo_macro <- read.csv(file.path(geo_decon_path, "geo_batched_macro_long.csv"))
# geo_macro <- geo_macro %>%
#   filter(section_id == "L3_1") # 237 -> 32

geo_decon_mapped <- rbind(
  geo_decon[geo_decon$sample %in% unique(mapped_all$Sample_ID), ] #,
  # geo_macro
  )%>%
  rename(sample_id = sample) # 431 -> 668 by adding 237 macros 
geo_decon_mapped$cell_fraction = ifelse(geo_decon_mapped$cell_fraction == "Macro", "Macrophage", geo_decon_mapped$cell_fraction)  


# > head(geo_decon_mapped) # 431   5
#                     sample CellType   Fraction section_id cell_fraction
# 14 DSP-1001660018473-B-D03     T_NK 0.18231696       B4_1         Other
# 16 DSP-1001660018473-B-D05     T_NK 0.17330728       B4_1     Malignant
# 17 DSP-1001660018473-B-D06     T_NK 0.19302811       B4_1         Other
# 19 DSP-1001660018473-B-D08     T_NK 0.26285257       B4_1         Other
# 21 DSP-1001660018473-B-D10     T_NK 0.14996134       B4_1         Other
# 24 DSP-1001660018473-B-E01     T_NK 0.09200995       B4_1         Other

# Vis decon results (spot)
vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2c"
# ## CARD
vis_decon <- rbind(read.csv(file.path(vis_decon_path, "vis_breast_decon_long.csv")),
                   read.csv(file.path(vis_decon_path, "vis_lung_decon_long.csv"))
                   )
# ## RCTD
# vis_decon <- rbind(read.csv(file.path(vis_decon_path, "vis_breast_decon_long_RCTD.csv")),
#                    read.csv(file.path(vis_decon_path, "vis_lung_decon_long_RCTD.csv"))
# )
vis_macro <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2d/vis_spot_70pct_myeloid_breast_lung_decon_long.csv")
vis_macro <- vis_macro %>%
  mutate(Sample_ID = "-",
         Section_ID = "L3_1",
         Section = "L3_2",
         cell_fraction = "Macrophage")

vis_decon_mapped_ <- mapped_all %>%
  rename(Barcode = barcode,
         Section = Section_ID_Vis) %>%
  left_join(vis_decon, by = c("Barcode", "Section")) %>%
  filter(cell_fraction != "PanCK-")
  
vis_decon_mapped <- rbind(vis_decon_mapped_
                          #, vis_macro
                          ) %>% # 3385 -> 3425 # increase 40 rows, 5 spots for macro
  rename(section_id = Section_ID,
       sample_id = Sample_ID) %>%
  mutate(cell_fraction = ifelse(cell_fraction == "T_cells", "T cells", cell_fraction)) %>%
  mutate(CellType = ifelse(CellType == "Fibro_Muscle", "Fibro_muscle", CellType))

# >  head(vis_decon_mapped)   # no macrophage: 3385    8
#              Barcode Section_ID               Sample_ID cell_fraction Section Region     CellType     Fraction
# 1 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE         T_NK 0.0002814366
# 2 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE            B 0.0248407509
# 3 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE        Tumor 0.1323340298
# 4 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE       Vessel 0.0085140267
# 5 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE Fibro_Muscle 0.8316757981
# 6 AAGGCCGATTCTGAGC-1       B1_1 DSP-1001660018473-B-E09         Other    B1_2  FALSE      Myeloid 0.0004073472
  
vis_geo_spot_decon <- rbind(
  geo_decon_mapped %>% select(CellType, Fraction, cell_fraction, section_id, sample_id) %>% 
    mutate(Platform = "GeoMx"),
  vis_decon_mapped %>% select(CellType, Fraction, cell_fraction, section_id, sample_id) %>% 
    mutate(Platform = "Visium")
)

vis_geo_spot_decon_test <- vis_geo_spot_decon %>%
  filter(is.na(CellType))

library(rstatix)
stat.test <- vis_geo_spot_decon %>%
  filter(cell_fraction == "Malignant") %>%
  group_by(CellType) %>%
  t_test(Fraction ~ Platform) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
#     CellType     .y.      group1 group2    n1    n2 statistic    df        p    p.adj p.adj.signif
#     <chr>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>       
#   1 B            Fraction GeoMx  Visium    28   117    -8.53  125.  4.12e-14 3.30e-13 ****       #  
#   2 Epithelia    Fraction GeoMx  Visium    17    73    -4.30   30.2 1.63e- 4 1.30e- 3 **         # 
#   3 Fibro_muscle Fraction GeoMx  Visium    28   117    -6.45  136.  1.79e- 9 1.43e- 8 ****       #  
#   4 Granulocyte  Fraction GeoMx  Visium    28   117    -1.37   31.4 1.79e- 1 1   e+ 0 ns          
#   5 Myeloid      Fraction GeoMx  Visium    28   117     0.438  59.2 6.63e- 1 1   e+ 0 ns          
#   6 T_NK         Fraction GeoMx  Visium    28   117     4.40   33.3 1.04e- 4 8.32e- 4 ***        #  
#   7 Tumor        Fraction GeoMx  Visium    28   117     5.26   74.1 1.35e- 6 1.08e- 5 ****       # 
#   8 Vessel       Fraction GeoMx  Visium    28   117    -0.497  47.0 6.21e- 1 1   e+ 0 ns     

stat.test <- vis_geo_spot_decon %>%
  filter(cell_fraction == "Other") %>%
  filter(CellType != "Epithelia") %>% # just one sample per platform
  group_by(CellType) %>%
  t_test(Fraction ~ Platform) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
#     CellType     .y.      group1 group2    n1    n2 statistic    df          p     p.adj p.adj.signif
#    <chr>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>      <dbl>     <dbl> <chr>       
#   1 B            Fraction GeoMx  Visium    28   350    -2.50   43.9 0.0164     0.115     ns          
#   2 Fibro_muscle Fraction GeoMx  Visium    28   350     1.19   71.7 0.237      1         ns          
#   3 Granulocyte  Fraction GeoMx  Visium    28   350    -5.69   33.9 0.00000219 0.0000153 ****        
#   4 Myeloid      Fraction GeoMx  Visium    28   350    -0.574  33.5 0.57       1         ns          
#   5 T_NK         Fraction GeoMx  Visium    28   350    -4.46   34.8 0.0000824  0.000577  ***         
#   6 Tumor        Fraction GeoMx  Visium    28   350    -1.44   76.4 0.155      1         ns          
#   7 Vessel       Fraction GeoMx  Visium    28   350     3.32   60.4 0.00151    0.0106    *      # 

stat.test <- vis_geo_spot_decon %>%
  filter(cell_fraction == "T cells") %>%
  group_by(CellType) %>%
  t_test(Fraction ~ Platform) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test
# CellType     .y.      group1 group2    n1    n2 statistic    df        p   p.adj p.adj.signif
#  <chr>        <chr>    <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>   <dbl> <chr>       
# 1 B            Fraction GeoMx  Visium     3     6     0.362  4.04 0.735    1       ns          
# 2 Fibro_muscle Fraction GeoMx  Visium     3     6    -0.823  6.86 0.438    1       ns          
# 3 Granulocyte  Fraction GeoMx  Visium     3     6    -4.66   5.82 0.00376  0.0263  *           
# 4 Myeloid      Fraction GeoMx  Visium     3     6     0.260  3.32 0.81     1       ns          
# 5 T_NK         Fraction GeoMx  Visium     3     6     0.261  5.56 0.804    1       ns          
# 6 Tumor        Fraction GeoMx  Visium     3     6    -0.630  6.91 0.549    1       ns          
# 7 Vessel       Fraction GeoMx  Visium     3     6     6.32   7.00 0.000397 0.00278 **           #

annotation_df <- data.frame(
  color = c("E", "H"),
  start = c("Good", "Fair"),
  end = c("Very Good", "Good"),
  y = c(3.6, 4.7),
  label = c("Comp. 1", "Comp. 2")
)

p <- ggplot(vis_geo_spot_decon, aes(x=CellType, y=Fraction, fill= Platform)) +
  geom_boxplot(outlier.size = 0.1) + 
  theme_bw() +
  scale_fill_manual(values = c("#92d051ff", "#257cf2ff")) + 
  theme(axis.text.x = element_text(size = 12.5, angle = 90, vjust = 0.5, hjust=1),
        panel.spacing=unit(1.5,"lines"),
        panel.grid = element_blank(), 
        strip.text.x = element_text(size = 13.5, face = "bold"), 
        strip.background=element_rect(fill="#DEDEDE")) + 
  facet_wrap(~cell_fraction, ncol = 4)

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures"
pdf(file = file.path(paste0(fig_path, "/Fig2d"), "Fig2d_AOI_facet.pdf"),
    width = 10,
    height = 5.5)
print(p)
dev.off()


