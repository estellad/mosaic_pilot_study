# Heatmap of Region versus sample id ----------------------------------
seu_sub <- seu 
seu_sub$sample_id <- paste0(substr(seu$sample_id, 1, 1), substr(seu$sample_id, 7, 7))
df <- as.data.frame(table(seu_sub$Region, seu_sub$sample_id))
df2 <- df %>%
  group_by(Var1) %>%
  mutate(countT= sum(Freq)) %>%
  mutate(per=round(100*Freq/countT,0))

data_frame <- data.frame(Region = df2$Var1, sample = df2$Var2, per = df2$per)
data_frame$sample <- factor(data_frame$sample, levels = rev(levels(factor(data_frame$sample))))

# Create the heatmap pt per Region
p <- ggplot(data_frame, aes(x = Region, y = sample, fill = per)) + 
  geom_tile() + 
  scale_x_discrete(position = "top") +
  scale_fill_viridis_c(option = "rocket", direction = -1) + 
  theme_minimal() +
  theme(legend.position = "bottom") + 
  labs(fill = "Per-region Percentage", title = "", 
       x = "Pathology Annotation", y = "Sample") + 
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                               barwidth = unit(23, "cm"), barheight = unit(0.5, "cm")))

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "/DLBCL_Vis_Pt_Region_Heatmap.pdf"), width = 10, height = 5)
print(p)
dev.off()