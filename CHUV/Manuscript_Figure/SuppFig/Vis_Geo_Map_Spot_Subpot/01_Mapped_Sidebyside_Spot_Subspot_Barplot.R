library(dplyr)
# Reading mapping results -------------------------------------------------
mapped_all_spot <- read.csv(file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration", "spot_mapped_cf_region.csv"))
mapped_all_spot$indication <- ifelse(substr(mapped_all_spot$Section_ID, 1, 1) == "B", "Breast", "Lung")
mapped_all_subspot <- read.csv(file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration", "subspot_mapped_cf_region.csv"))
mapped_all_subspot$indication <- ifelse(substr(mapped_all_subspot$Section_ID, 1, 1) == "B", "Breast", "Lung")

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/SuppFig_regi_subspot"
###########################################################################
# Frequency Table spot subspot by indication and CF -----------------------
table_spot <- table(mapped_all_spot$indication, mapped_all_spot$cell_fraction)
table_subspot <- table(mapped_all_subspot$indication, mapped_all_subspot$cell_fraction)

plt_df_spot <- data.frame(table_spot[, c("Malignant", "Other", "T_cells")]) %>% 
  rename(Indication = Var1,
         cell_fraction = Var2,
         Count = Freq)%>%
  mutate(cell_fraction = factor(cell_fraction, levels = c("T_cells", "Malignant", "Other")),
         Resolution = "Spot")

plt_df_subspot <- data.frame(table_subspot[, c("Malignant", "Other", "T_cells")]) %>% 
  rename(Indication = Var1,
         cell_fraction = Var2,
         Count = Freq) %>%
  mutate(cell_fraction = factor(cell_fraction, levels = c("T_cells", "Malignant", "Other")),
         Resolution = "Subspot")

plt_df <- rbind(plt_df_spot, plt_df_subspot)


p <- ggplot(mapping = aes(y = Count, linetype = Resolution, fill = Indication)) + 
  geom_col(data = plt_df %>% filter(Resolution == "Spot"), aes(x = as.numeric(cell_fraction) + 0.15), width = 0.28, color = "black") +
  geom_text(data = plt_df %>% filter(Resolution == "Spot") %>% filter(Count > 0), aes(x = as.numeric(cell_fraction) + 0.15, label = Count), hjust=1, vjust=-3.5, size = 3.5, position = "stack") + 
  geom_col(data = plt_df %>% filter(Resolution == "Subspot"), aes(x = as.numeric(cell_fraction) - 0.15), width = 0.28, color = "black") +
  geom_text(data = plt_df %>% filter(Resolution == "Subspot") %>% filter(Count > 0), aes(x = as.numeric(cell_fraction) - 0.15, label = Count), hjust=1, vjust=4.5, size = 3.5, position = "stack") + 
  scale_fill_manual(name = "Indication", labels = c("Breast", "Lung"), values = c("#00c78cff", "#ffa54fff")) + 
  scale_linetype_manual(name = "Resolution", labels = c("Spot", "Subspot"), values = c("solid", "dashed")) +
  scale_x_continuous(breaks = c(3.0, 2.0, 1.0),
                   labels = c("Other", "Malignant", "T cells")) + 
  coord_flip() + 
  labs(x = "Cell fraction") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(), 
        legend.spacing.y = unit(1.0, 'cm'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"),
        legend.key.size = unit(2, "lines"), 
        axis.text.y = element_text(size = 13, angle = 90, vjust = 1, hjust=0.5),
        axis.title = element_text(size = 14)) + 
  guides(linetype = guide_legend(byrow = TRUE, override.aes = list(fill = c(NA, NA))),
         fill = guide_legend(override.aes = list(color = c(NA, NA))))


pdf(file = file.path(fig_path, "regi_count_stacked_indication_side_by_side_spot_subspot.pdf"),
    width = 12.5,
    height = 7)
print(p)
dev.off()


###########################################################################
# Side by side spot subspot  ----------------------------------------------
patho_order <- c("Tumor_pure", "Tumor_Stroma_mix", "Most_likely_tumor",
                 "Intratumoral_Stroma", "Intratumoral_Vessel", "Acellular mucin", 
                 "Tumor_Immune_mix", "Immune_Cell_mix", "Lymphocytes")

get_mapped_plt_df <- function(mapped_all_df = mapped_all_subspot, resolution = "Subspot"){
  mapped_all_plt <- mapped_all_df %>% # 4728
    select(cell_fraction, Region) %>%
    filter(!(Region %in% c("NA", "FALSE"))) %>% # 4164
    filter(cell_fraction != "PanCK-") %>% # 3192
    group_by(cell_fraction, Region) %>%
    summarise(Count = n()) %>% # 22 combinations
    ungroup() %>%
    mutate(prop = Count/sum(Count),
           Resolution = resolution) %>%
    filter(Region != "Necrosis_Debris") %>% # 21 combinations
    mutate(Region = factor(Region, 
                           levels = patho_order))
  
  return(mapped_all_plt)
}

mapped_all_plt <- rbind(get_mapped_plt_df(mapped_all_df = mapped_all_spot, resolution = "Spot"), 
                        get_mapped_plt_df(mapped_all_df = mapped_all_subspot, resolution = "Subspot"))

patho_anno <- c(
  "Tumor_pure",
  "Most_likely_tumor",
  "Tumor_Immune_mix",
  "Lymphocytes",
  "Immune_Cell_mix",
  "Tumor_Stroma_mix",
  "Intratumoral_Stroma",
  "Intratumoral_Vessel",
  "Acellular mucin"
)

patho_color <- c(
  "#FFA54F",
  "#c47f3d",
  "#00e04e",
  "#4169E1",
  "#b2ff00",
  "#808000",# "#00c3c3",
  "#388E8E",
  "#ba5440",
  "#FFC1C1"
)

names(patho_color) <- patho_anno

patho_color <- patho_color[patho_order]


mapped_all_plt <- rbind(mapped_all_plt, 
                        data.frame(cell_fraction = c("Malignant", rep("T_cells", 4)),
                                   Region = c("Immune_Cell_mix", "Tumor_Stroma_mix", "Intratumoral_Stroma", "Tumor_Immune_mix", "Immune_Cell_mix"),
                                   Count = rep(0, 5),
                                   prop = rep(0, 5),
                                   Resolution = rep("Spot", 5)))

bar_width <- 0.8
plot_each_region <- function(mapped_all_plt_subset){
  ggplot(data = mapped_all_plt_subset, aes(x = Region, y = prop, fill = Region, linetype = Resolution, group = Resolution)) +
    geom_bar(stat = "identity", position = position_dodge2(width = bar_width, preserve = "single"), color = "black", width = bar_width) +
    scale_fill_manual(values = patho_color) + 
    scale_linetype_manual(name = "Resolution", labels = c("Spot", "Subspot"), values = c("solid", "dashed")) +
    theme_bw() +
    theme(axis.text.x = element_blank(), # element_text(size = 13, angle = 90, vjust = 0.5, hjust=1),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(), 
          legend.position =  "none",
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          plot.title = element_text(size = 18),
          axis.title = element_text(size = 14)) + 
    labs(y = "Proportion", x = "") + 
    # facet_wrap(~cell_fraction) +
    guides(linetype = guide_legend(override.aes = list(fill = c(NA, NA))),
           fill = guide_legend(override.aes = list(color = rep(NA, length(unique(mapped_all_plt_subset$Region))))))
}

p1 <- plot_each_region(mapped_all_plt_subset = mapped_all_plt %>% filter(cell_fraction == "Malignant")) + ggtitle("Malignant") + theme(plot.title = element_text(hjust = 0.5, size = 17))
p2 <- plot_each_region(mapped_all_plt_subset = mapped_all_plt %>% filter(cell_fraction == "Other")) + ggtitle("Other") + theme(plot.title = element_text(hjust = 0.5, size = 17)) + ylab(" ")
p3 <- plot_each_region(mapped_all_plt_subset = mapped_all_plt %>% filter(cell_fraction == "T_cells")) + ggtitle("T cells") + theme(plot.title = element_text(hjust = 0.5, size = 17)) + ylab(" ")
p4 <- plot_each_region(mapped_all_plt_subset = mapped_all_plt) + 
  theme(legend.position = "bottom", # for the complete legend
        legend.spacing.y = unit(1.0, 'cm'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"),
        legend.key.size = unit(2, "lines")) + 
  guides(colour = guide_legend(nrow = 1))


p_combined <- p1 | p2 | p3
pdf(file = file.path(fig_path, "AOI_spot_subspot_minimal_bar.pdf"),
    width = 15,
    height = 7.5)
print(p_combined)
dev.off()

pdf(file = file.path(fig_path, "forlegend_all_spot_subspot.pdf"),
    width = 20,
    height = 5)
print(p4)
dev.off()



