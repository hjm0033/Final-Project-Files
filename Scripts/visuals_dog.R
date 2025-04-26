## ──────────────────────────────── 0.  Setup ──────────────────────────────── ##
# Run this *once* per machine to install packages
# if (!requireNamespace("DESeq2", quietly = TRUE)) {
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("DESeq2")
# }
# if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
# if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
# if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
#     if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#     BiocManager::install("EnhancedVolcano")
# }


library(DESeq2)
library(pheatmap)        # heat‐maps
library(RColorBrewer)    # colour palettes
library(ggplot2)         # Enhanced plotting
library(ggrepel)         # For non-overlapping labels
library(EnhancedVolcano) # For Volcano plots

# Set working directory if needed
# setwd("~/Downloads") # Or use RStudio ▶ Session ▶ Set WD ▶ To Source File Dir

## ───────────────────── 5A. Publication Volcano Plot ─────────────────────── ##

# Prepare data frame for EnhancedVolcano
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df) # Add gene names as a column

# Define significance thresholds
p_cutoff <- 0.05
fc_cutoff <- 1.0 # Log2 fold change threshold (e.g., 1.0 means 2-fold change)

# --- FIX IS HERE ---
# Identify significant genes, explicitly handling potential NAs using which()
# which() returns the indices of rows where the condition is TRUE (and ignores NAs)
sig_indices <- which(res_df$padj < p_cutoff & abs(res_df$log2FoldChange) > fc_cutoff)

# Get the rownames for those significant genes
sig_genes_rownames <- rownames(res_df[sig_indices, ])

# Select the top 15 significant genes (e.g., by adjusted p-value) to label
# First, order the significant genes by p-value
ordered_sig_genes <- res_df[sig_indices, ]
ordered_sig_genes <- ordered_sig_genes[order(ordered_sig_genes$padj), ]
genes_to_label <- head(rownames(ordered_sig_genes), 15)
# --- END FIX ---


# Create the Volcano Plot
volcano_plot <- EnhancedVolcano(res_df,
                                lab = res_df$gene, # Use gene names for labels
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Differential Expression (Large vs Small)',
                                pCutoff = p_cutoff,
                                FCcutoff = fc_cutoff,
                                pointSize = 2.0,
                                labSize = 3.0,
                                colAlpha = 0.8, # Point transparency
                                legendPosition = 'right',
                                legendLabSize = 12,
                                legendIconSize = 4.0,
                                # Use the safely selected gene names
                                selectLab = genes_to_label,
                                drawConnectors = TRUE, # Draw lines to labels
                                widthConnectors = 0.5,
                                colConnectors = 'grey50',
                                boxedLabels = TRUE, # Put labels in boxes
                                subtitle = paste0('Cutoffs: adj. p < ', p_cutoff, ', |log2FC| > ', fc_cutoff)
                                # Add more customizations as needed:
                                # xlim = c(-8, 8),
                                # ylim = c(0, -log10(1e-50)), # Adjust y-limit if needed
                                # gridlines.major = FALSE,
                                # gridlines.minor = FALSE
)

# Print the plot to the viewer
print(volcano_plot)

# Save the plot to a file (e.g., PDF)
ggsave("19_Dog_Volcano_Plot.pdf", plot = volcano_plot, width = 10, height = 8, units = "in", dpi = 300)
ggsave("19_Dog_Volcano_Plot.png", plot = volcano_plot, width = 10, height = 8, units = "in", dpi = 300)

# Optional: Check for NAs if curious
# summary(res_df$padj)
# summary(res_df$log2FoldChange)

# Optional: Keep the original MA plot if desired
# pdf("19_Dog_MA_Plot.pdf", width = 7, height = 5)
# plotMA(res, ylim = c(-8, 8), main = "MA Plot")
# dev.off()

## ────────────────── 5B. Publication Sample Distance Heatmap ───────────────── ##

# Ensure rld is calculated (from original script section 5)
# rld <- rlog(dds)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

# Use Size information for row labels, keep original sample IDs for columns if preferred
# Or use Size for both if that makes sense for your comparison
rownames(sampleDistMatrix) <- paste(coldata$Size, rownames(coldata), sep="-") # More informative labels
colnames(sampleDistMatrix) <- NULL # Hide column names for cleaner look if row names are sufficient

# Define colors
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Create annotation data frame (optional, but good practice)
df_annot <- data.frame(Size = coldata$Size, row.names = rownames(sampleDistMatrix))

# Generate the heatmap
pheatmap_dist <- pheatmap(sampleDistMatrix,
                          clustering_distance_rows = sampleDists,
                          clustering_distance_cols = sampleDists,
                          col = colors,
                          annotation_row = df_annot, # Add annotation bar for rows
                          # annotation_col = df_annot, # Add annotation bar for columns if needed
                          main = "Sample-to-Sample Distances (based on rlog transformed data)",
                          fontsize = 10,
                          fontsize_row = 8, # Adjust row label size
                          # cellwidth = 15, # Uncomment and adjust for fixed cell size if needed
                          # cellheight = 15, # Uncomment and adjust for fixed cell size if needed
                          border_color = "grey60" # Add subtle borders
)

# Save the heatmap to a file
# Use the pheatmap object directly if saving via ggsave (requires converting)
# Or use standard R devices like pdf() or png() which often work better for pheatmap

# Option 1: Using pdf()/png() devices (Recommended for pheatmap)
pdf("19_Dog_Sample_Distance_Heatmap.pdf", width = 8, height = 7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         annotation_row = df_annot,
         main = "Sample-to-Sample Distances (rlog transformed)",
         fontsize = 10,
         fontsize_row = 8,
         border_color = "grey60")
dev.off()

png("19_Dog_Sample_Distance_Heatmap.png", width = 8, height = 7, units = "in", res = 300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         annotation_row = df_annot,
         main = "Sample-to-Sample Distances (rlog transformed)",
         fontsize = 10,
         fontsize_row = 8,
         border_color = "grey60")
dev.off()

# Note: The plot might look slightly different when saved vs viewed in RStudio plot pane.
# Adjust width/height/fontsize in the pdf()/png() calls as needed.

## ───────────────────── 5C. Publication PCA Plot ─────────────────────────── ##
# Ensure rld is calculated
# rld <- rlog(dds)

# Use DESeq2's plotPCA function to get the data
pcaData <- plotPCA(rld, intgroup = "Size", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar")) # Get variance explained

# Add other metadata if desired (e.g., sex) by merging
pcaData <- merge(pcaData, coldata[, c("sex"), drop = FALSE], by.x="name", by.y="row.names")

# Create the ggplot object
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Size, shape = sex)) + # Color by Size, shape by sex
  geom_point(size = 4, alpha = 0.8) + # Adjust point size and transparency
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  # Add labels using ggrepel to avoid overlap
  geom_text_repel(aes(label = name), size = 3, max.overlaps = Inf) +
  ggtitle("PCA Plot (rlog transformed data)") +
  # Use a clean theme
  theme_bw(base_size = 14) + # Adjust base font size
  theme(
    plot.title = element_text(hjust = 0.5), # Center title
    legend.position = "right" # Adjust legend position
  ) +
  # Optionally add ellipses around groups
  # stat_ellipse(aes(group = Size), type = "t", level = 0.95) +
  # Customize colors if needed
  scale_color_brewer(palette = "Set1") # Example using a Brewer palette

# Print the plot
print(pca_plot)

# Save the plot
ggsave("19_Dog_PCA_Plot.pdf", plot = pca_plot, width = 8, height = 6, units = "in", dpi = 300)
ggsave("19_Dog_PCA_Plot.png", plot = pca_plot, width = 8, height = 6, units = "in", dpi = 300)


## ───────────── 5D. Publication Top Variable Genes Heatmap (with Gene Names) ──── ##
# Ensure vsd and coldata are calculated/available

# Select top 50 variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat_vsd <- assay(vsd)[topVarGenes, ]
# Center the data row-wise
mat_vsd <- mat_vsd - rowMeans(mat_vsd)

# Create annotation data frame for columns
df_annot_col <- coldata[, c("Size", "sex"), drop = FALSE]

# Check unique sex values (as done previously)
print("Unique values in coldata$sex being used for heatmap annotation:")
actual_sex_levels <- unique(df_annot_col$sex)
print(actual_sex_levels)

# Define annotation colors using the actual values
# *** Make sure this matches your actual levels ("Nan", "female", "male", etc.) ***
ann_colors = list(
  Size = c(Small="lightblue", Large="darkblue"),
  sex = c(
    Nan = "grey",
    female = "purple",
    male = "orange"
    # Add/modify based on your unique(coldata$sex) output
  )
)
if (!all(actual_sex_levels %in% names(ann_colors$sex))) {
  warning("Not all levels in coldata$sex have colours defined in ann_colors! Please update ann_colors.")
}

# --- CHANGES FOR GENE NAMES ---
# Set show_rownames = TRUE
# Adjust fontsize_row as needed
# Increase height in pdf()/png() calls
# ---

# Generate the heatmap (optional interactive display)
# This might look crowded in the RStudio plot pane
pheatmap_topvar <- pheatmap(mat_vsd,
                            annotation_col = df_annot_col,
                            annotation_colors = ann_colors,
                            main = "Top 50 Most Variable Genes (VST transformed, centered)",
                            show_rownames = TRUE,   # <--- CHANGE: Show gene names
                            show_colnames = TRUE,
                            fontsize = 10,
                            fontsize_row = 8,    # <--- CHANGE: Adjust gene name font size (start with 8)
                            fontsize_col = 8,
                            border_color = "grey60",
                            color = colorRampPalette(c("blue", "white", "red"))(100)
)

# Save the heatmap (using pdf/png is recommended)
pdf("19_Dog_Top50_Var_Genes_Heatmap_WithNames.pdf", width = 9, height = 12) # <--- CHANGE: Increased height
pheatmap(mat_vsd,
         annotation_col = df_annot_col,
         annotation_colors = ann_colors,
         main = "Top 50 Most Variable Genes (VST centered)",
         show_rownames = TRUE,   # <--- CHANGE: Show gene names
         show_colnames = TRUE,
         fontsize = 10,
         fontsize_row = 8,    # <--- CHANGE: Adjust gene name font size
         fontsize_col = 8,
         border_color = "grey60",
         color = colorRampPalette(c("blue", "white", "red"))(100)
)
dev.off()

png("19_Dog_Top50_Var_Genes_Heatmap_WithNames.png", width = 9, height = 12, units = "in", res = 300) # <--- CHANGE: Increased height
pheatmap(mat_vsd,
         annotation_col = df_annot_col,
         annotation_colors = ann_colors,
         main = "Top 50 Most Variable Genes (VST centered)",
         show_rownames = TRUE,   # <--- CHANGE: Show gene names
         show_colnames = TRUE,
         fontsize = 10,
         fontsize_row = 8,    # <--- CHANGE: Adjust gene name font size
         fontsize_col = 8,
         border_color = "grey60",
         color = colorRampPalette(c("blue", "white", "red"))(100)
)
dev.off()


## ───────────────────── 5E. Publication MA Plot (using ggplot2) ──────────── ##

# Ensure the DESeq2 results object 'res' exists and ggplot2 is loaded
# library(ggplot2)

# --- Data Preparation ---
# Convert results to a data frame
res_df_ma <- as.data.frame(res)
# Add gene names as a column (if not already present from previous steps)
if (!"gene" %in% colnames(res_df_ma)) {
  res_df_ma$gene <- rownames(res_df_ma)
}

# Define significance threshold
padj_threshold <- 0.05

# Create a column for significance status, handling NAs in padj
res_df_ma$significant <- ifelse(!is.na(res_df_ma$padj) & res_df_ma$padj < padj_threshold,
                                "Significant",
                                "Not Significant")
# Make it a factor for consistent coloring/legend
res_df_ma$significant <- factor(res_df_ma$significant, levels = c("Significant", "Not Significant"))

# Optional: Add labels for top genes (e.g., by smallest padj)
# Filter out NAs first for reliable ordering
res_df_ma_filtered <- na.omit(res_df_ma)
res_df_ma_filtered <- res_df_ma_filtered[order(res_df_ma_filtered$padj), ]
# Select top N significant genes to label (adjust N as needed)
top_genes_to_label <- head(res_df_ma_filtered[res_df_ma_filtered$significant == "Significant", "gene"], 10)
res_df_ma$label <- ifelse(res_df_ma$gene %in% top_genes_to_label, res_df_ma$gene, "")

# --- Create ggplot ---
ma_plot_gg <- ggplot(res_df_ma, aes(x = baseMean, y = log2FoldChange)) +
  # Add points, coloring by significance, with transparency
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  
  # Use log10 scale for the x-axis (mean expression) - common for MA plots
  # Add a small pseudocount (e.g., 1) to baseMean to avoid log10(0) issues if any genes have 0 mean
  scale_x_log10(name = "Mean of Normalized Counts (log10 scale)") +
  # Add log tick marks
  annotation_logticks(sides = "b", colour = "grey50") +
  
  # Set y-axis label
  scale_y_continuous(name = "Log2 Fold Change (Large vs Small)") +
  
  # Define custom colors for significance
  scale_color_manual(name = "Significance",
                     values = c("Significant" = "red", "Not Significant" = "grey50"),
                     labels = c(paste0("padj < ", padj_threshold), paste0("padj >= ", padj_threshold))) +
  
  # Add a horizontal line at y=0 (no fold change)
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkblue") +
  
  # Optional: Add labels to top genes using ggrepel (ensure library(ggrepel) is loaded)
  # geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf,
  #                 box.padding = 0.5, point.padding = 0.5,
  #                 segment.color = 'grey50', segment.size = 0.5) +
  
  # Add Title
  ggtitle("MA Plot: Differential Expression (Large vs Small)") +
  
  # Apply a clean theme
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center title
    legend.position = "bottom",            # Adjust legend position
    panel.grid.major = element_line(colour = "grey90"), # Lighter grid lines
    panel.grid.minor = element_blank()
  ) +
  # Optional: Adjust y-axis limits if needed (use coord_cartesian to zoom without removing data)
  coord_cartesian(ylim = c(-8, 8)) # Adjust limits based on your data range

# Print the plot to the viewer
print(ma_plot_gg)

# --- Save the plot ---
ggsave("19_Dog_MA_Plot_ggplot.png", plot = ma_plot_gg, width = 8, height = 7, units = "in", dpi = 300)

# Optional: Clean up intermediate objects if needed
# rm(res_df_ma, res_df_ma_filtered, top_genes_to_label, ma_plot_gg)

message("✓ MA plot created and saved using ggplot2.")

# Note: You might see warnings like "Removed X rows containing missing values"
# or "Transformation introduced infinite values in continuous x-axis". This is often
# expected if some genes have baseMean of 0 (log10(0) = -Inf) or NA p-values.
# The plot should still render correctly for the valid points.