library(ggplot2)

# Get top 10 up-regulated and down-regulated genes
topUpGenes <- rownames(res_characterized)[order(-res_characterized$log2FoldChange, res_characterized$padj)][1:10]
topDownGenes <- rownames(res_characterized)[order(res_characterized$log2FoldChange, res_characterized$padj)][1:10]

# Create a dataframe for plotting
top_genes_df <- rbind(
  data.frame(
    Gene = topUpGenes,
    log2FC = res_characterized[topUpGenes, "log2FoldChange"],
    padj = res_characterized[topUpGenes, "padj"],
    Direction = "Up-regulated"
  ),
  data.frame(
    Gene = topDownGenes,
    log2FC = res_characterized[topDownGenes, "log2FoldChange"],
    padj = res_characterized[topDownGenes, "padj"],
    Direction = "Down-regulated"
  )
)

# Reorder genes by fold change
top_genes_df$Gene <- factor(top_genes_df$Gene, levels = top_genes_df$Gene[order(top_genes_df$log2FC)])

# Create the dot plot
ggplot(top_genes_df, aes(x = log2FC, y = Gene, color = Direction, size = -log10(padj))) +
  geom_point() +
  scale_color_manual(values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8")) +
  scale_size_continuous(range = c(2, 8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  labs(
    title = "Top Differentially Expressed Genes",
    subtitle = "Large vs Small Dogs",
    x = "Log2 Fold Change",
    y = "",
    size = "-log10(padj)",
    color = "Regulation"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Save as high-resolution PDF
ggsave("top_genes_dotplot.pdf", width = 10, height = 8, dpi = 300)