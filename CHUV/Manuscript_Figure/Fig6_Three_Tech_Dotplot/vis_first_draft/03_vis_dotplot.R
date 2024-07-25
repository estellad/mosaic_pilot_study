# Dot plot ----------------------------------------------------------------
DLBCLnChromium_Marker_Gene_List <- c(
  "MS4A1", "TNFRSF13C", "CD79B", "CD37", "PSMB8", "CD19", "TYMS",
  "TUBB", "TOP2A", "POLD4", "CD47", "CD52", "BLK", "CD38", "MAP2K1",
  "CD40", "BCL2L1", "TNFRSF8", "SMO", "RARA", "TYK2", "TNFRSF10B"
)

vis <- vis[rownames(vis) %in% DLBCLnChromium_Marker_Gene_List, ]

vis_small <- vis[, vis$decon_max != "Mix"]

p <- DotPlot(
  vis_small,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5))

plot_title = "Geo_DLBLCnChromium_dotplot.pdf"
pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig6_drug_target_dlbcl_vis/", plot_title),
    width = 5,
    height = 20)
print(p)
dev.off()