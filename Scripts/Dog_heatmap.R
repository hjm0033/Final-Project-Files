library(pheatmap)
library(viridis)

# Calculate sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

# Create proper row names
rownames(sampleDistMatrix) <- paste0(
  colnames(rld), 
  " (", colData(rld)$Size, 
  ", ", colData(rld)$sex, ")"
)
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)

# Create annotation for the heatmap
annotation_row <- data.frame(
  Size = colData(rld)$Size,
  Sex = colData(rld)$sex,
  row.names = rownames(sampleDistMatrix)
)

# Define custom colors
anno_colors <- list(
  Size = c(Small = "#1B9E77", Large = "#D95F02"),
  Sex = c(male = "#7570B3", female = "#E7298A", Nan = "#66A61E")
)

# Create the enhanced heatmap
pheatmap(
  sampleDistMatrix,
  main = "Sample-to-Sample Distance",
  annotation_row = annotation_row,
  annotation_col = annotation_row,
  annotation_colors = anno_colors,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = viridis(100, direction = -1), # Using viridis color palette
  border_color = NA,
  cellwidth = 15,
  cellheight = 15,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  filename = "sample_distance_heatmap.pdf",
  width = 10, 
  height = 9
)