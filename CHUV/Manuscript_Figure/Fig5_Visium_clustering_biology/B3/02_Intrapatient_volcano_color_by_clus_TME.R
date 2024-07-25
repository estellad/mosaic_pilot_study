contrast <- "1&5&9_14"
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2"
DE_result <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))

results <- as.data.frame(DE_result)

results$avg_log2FC <- ifelse(results$cluster == "14_TME", -1 * results$avg_log2FC, results$avg_log2FC) # upper tissue location one
# results$avg_log2FC <- ifelse(results$cluster == "9", -1 * results$avg_log2FC, results$avg_log2FC)
results = results %>% arrange(p_val_adj)
head(results)

# Categorize results based on P-value & FDR for plotting
results$Color[results$cluster == "159_TME" & results$p_val_adj < 0.05] <- "Clus159_05"
results$Color[results$cluster == "159_TME" & results$p_val_adj < 0.01] <- "Clus159_01"
results$Color[results$cluster == "14_TME" & results$p_val_adj < 0.05] <- "Clus14_05"
results$Color[results$cluster == "14_TME" & results$p_val_adj < 0.01] <- "Clus14_01"
# results$Color[abs(results$avg_log2FC) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                        levels = c("Clus159_05", "Clus159_01",
                                   "Clus14_05", "Clus14_01"))

top_n_genes = 12
# results$gene <- rownames(results)
# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$p_val)) * sign(results$avg_log2FC)
top_g = c(results[, 'gene'][
  order(results[, "invert_P"], decreasing = TRUE)[1:top_n_genes]],
  results[, 'gene'][
    order(results[, 'invert_P'], decreasing = FALSE)[1:top_n_genes]])
top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix

# Graph results
p <- ggplot(results,
            aes(x = avg_log2FC, y = -log10(p_val),
                color = Color, label = gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point(size = 5) +
  labs(x = expression("log"[2]*"FoldChange"),
       y = expression("Significance, -log"[10]*"p-value"),
       color = "Significance") +
  scale_color_manual(values = c("Clus159_05" = "#fc98fc",
                                "Clus159_01" = "#ff00ff",
                                "Clus14_05" = "#9393ff",
                                "Clus14_01" = "#0000ff"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  ggrepel::geom_text_repel(data = subset(results, gene %in% top_g & p_val_adj < 0.001),
                           size = 8, point.padding = 0.15, color = "black",
                           min.segment.length = .1, box.padding = .2, lwd = 2,
                           max.overlaps = 50) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5/B3"

plot_title = "Vis_B3_2_DE_Volcano_159_14_.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 10)
print(p)
dev.off()

# ggsave(file.path(file.path(out_dir, "volcanoplot.png")))