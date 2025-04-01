library(pheatmap)
library(viridis)
library(RColorBrewer)
library(DESeq2)

# Get the top 50 most variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

# Extract data for these genes
mat <- assay(vsd)[topVarGenes, ]

# Center each row (gene) by subtracting the row mean
mat_centered <- mat - rowMeans(mat)

# Create annotation dataframe for samples
annotation_col <- data.frame(
  Size = colData(vsd)$Size,
  Sex = colData(vsd)$sex,
  row.names = colnames(vsd)
)

# Define custom colors for annotations
anno_colors <- list(
  Size = c(Small = "#1B9E77", Large = "#D95F02"),
  Sex = c(male = "#7570B3", female = "#E7298A", Nan = "#66A61E")
)

# Get gene names (assuming they're in the rownames)
gene_names <- rownames(mat)

# Create the enhanced heatmap
pheatmap(
  mat_centered,
  main = "Top 50 Most Variable Genes Across Samples",
  annotation_col = annotation_col,
  annotation_colors = anno_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = viridis(100),  # Using viridis color palette for better visualization
  border_color = NA,
  cellwidth = 15,
  cellheight = 10,
  treeheight_row = 30,
  treeheight_col = 30,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  filename = "highly_variable_genes_heatmap.pdf",  # Save directly as PDF
  width = 10, 
  height = 12
)

# If you want additional control or to save in different formats,
# you can also save using pdf() function
pdf("highly_variable_genes_heatmap_alt.pdf", width = 10, height = 12)
pheatmap(
  mat_centered,
  main = "Top 50 Most Variable Genes Across Samples",
  annotation_col = annotation_col,
  annotation_colors = anno_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = viridis(100),
  border_color = NA,
  cellwidth = 15,
  cellheight = 10,
  treeheight_row = 30,
  treeheight_col = 30,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2"
)
dev.off()

# To create a PNG version for presentations or quick viewing
png("highly_variable_genes_heatmap.png", width = 10, height = 12, units = "in", res = 300)
pheatmap(
  mat_centered,
  main = "Top 50 Most Variable Genes Across Samples",
  annotation_col = annotation_col,
  annotation_colors = anno_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  color = viridis(100),
  border_color = NA,
  cellwidth = 15,
  cellheight = 10,
  treeheight_row = 30,
  treeheight_col = 30,
  clustering_distance_rows = "euclidean", 
  clustering_method = "ward.D2"
)
dev.off()