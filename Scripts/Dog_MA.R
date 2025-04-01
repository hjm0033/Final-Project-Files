library(ggplot2)
library(ggrepel)

# Extract the data from the DESeq2 results
ma_data <- data.frame(
  mean = res_characterized$baseMean,
  lfc = res_characterized$log2FoldChange,
  padj = res_characterized$padj,
  sig = ifelse(res_characterized$padj < 0.05, "Significant", "Not Significant"),
  gene = rownames(res_characterized)
)

# Find top genes to label (most significant or highest fold changes)
top_genes <- ma_data[order(ma_data$padj),][1:20,]

# Create the MA plot
ggplot(ma_data, aes(x = log10(mean), y = lfc, color = sig)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    title = "MA Plot: Large vs Small Dogs",
    x = "Log10 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Save as high-resolution PDF
ggsave("ma_plot_dog_size.pdf", width = 9, height = 7, dpi = 300)