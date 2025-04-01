library(ggplot2)
library(ggrepel)

# Extract PCA data
pcaData <- plotPCA(rld, intgroup = "Size", returnData = TRUE)

# Calculate the percentage of variance explained by PC1 and PC2
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create PCA plot with only Size
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Size)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(
    aes(label = name),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  stat_ellipse(aes(fill = Size), geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE) +
  labs(
    title = "PCA Analysis of Gene Expression",
    subtitle = "Based on rlog-transformed counts",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.box = "vertical"
  )

# Display the plot
print(pca_plot)

# Save the plot as high-resolution PDF
ggsave("pca_plot_size_only.pdf", plot = pca_plot, width = 10, height = 8, dpi = 300)
