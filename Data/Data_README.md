# Dog Brain RNA-Seq Data Files

This directory contains data files generated during the dog brain RNA-seq analysis project. These files include raw and processed data, analysis results, and quality control reports.

## Data Categories

### Main Analysis Results

- **gene_count_matrix.csv**: Gene-level count matrix generated from StringTie analysis, containing raw counts for all samples.
- **transcript_count_matrix.csv**: Transcript-level count matrix with counts for individual transcripts.
- **Dog_NormTransExpIDs.txt**: Normalized expression data for downstream analysis, transformed to stabilize variance.
- **Dog_DGESeq_results.csv**: Complete differential expression results from DESeq2 analysis comparing large vs. small dog brains.

### Analysis of Specific Samples

- **19_Dog_DESeq2_rank.rnk**: Ranked gene list for GSEA analysis from sample 19, with genes ranked by their differential expression significance.
- **19_Dog_DESeq2_results.csv**: Differential expression results for sample 19, including log2 fold changes and adjusted p-values.
- **19_Dog_normalised_expression.txt**: Normalized expression data specific to sample 19.
- **20_Dog_DESeq2_rank.rnk**: Ranked gene list for GSEA analysis from sample 20.
- **20_Dog_DESeq2_results.csv**: Differential expression results for sample 20.
- **20_Dog_normalised_expression.txt**: Normalized expression data specific to sample 20.

### Differential Expression Results

- **all_significant_DEGs.txt**: Comprehensive list of all statistically significant differentially expressed genes (FDR < 0.05).
- **top25_upregulated.txt**: List of the 25 most significantly upregulated genes in large dog brains.
- **top25_downregulated.txt**: List of the 25 most significantly downregulated genes in large dog brains.
- **top50_DEGs.txt**: Top 50 differentially expressed genes ranked by significance.
- **top50_variable_genes.txt**: Top 50 genes with highest variance across all samples.

### Gene Set Enrichment Analysis (GSEA) Results

- **Dog_DGErankName.rnk**: Ranked gene list used as GSEA input, ranking genes based on their differential expression.
- **Dog_GSEA_KEGG.GseaPreranked.1743603144299.rpt**: GSEA report showing enrichment for KEGG pathways.
- **gsea_report_for_na_pos_1743603144299.html**: HTML report for positively enriched gene sets (upregulated in large dogs).
- **gsea_report_for_na_neg_1743603144299.html**: HTML report for negatively enriched gene sets (downregulated in large dogs).
- **gene_set_sizes.tsv**: Information about the sizes of analyzed gene sets.
- **ranked_gene_list_na_pos_versus_na_neg_1743603144299.tsv**: Full ranked gene list comparing positive vs. negative enrichment.
- **pos_snapshot.html** and **neg_snapshot.html**: HTML snapshots of positively and negatively enriched pathways.

### IIS-Specific Analysis

- **gene_set_sizesIIS.rpt** and **gene_set_sizesIIS.tsv**: Reports for IIS-specific gene sets with information about their sizes.
- **gsea_report_for_na_pos_IIS_1743697612355.tsv** and **gsea_report_for_na_neg_IIS_1743697612355.tsv**: GSEA results for IIS-specific pathways.
- **pos_snapshot_IIS.html** and **neg_snapshot_IIS.html**: HTML snapshots of IIS-specific pathway enrichment visualization.
- **ranked_gene_list_IIS_na_pos_versus_na_neg_1743697612355.tsv**: Ranked gene list specific to the IIS analysis.

### Quality Control Reports

- **multiqc_report_original.html**: Initial MultiQC report showing quality metrics of the raw sequencing data.
- **multiqc_report_reduced_headcrop.html**: Quality report after headcrop trimming, showing improved quality metrics.
- **multiqc_report_stringent_window.html**: Quality report after stringent window trimming, showing final read quality.

### Metadata and Statistics

- **Proposed_Samples_Meta_Data.txt**: Comprehensive metadata for all samples used in the analysis, including size, sex, and other relevant factors.
- **overall_statistics.txt**: Summary statistics from the analysis, including mapping rates and count distributions.

## Usage Notes

1. For differential gene expression analysis, start with the **Dog_DGESeq_results.csv** file.
2. For pathway analysis, refer to the GSEA results in HTML format.
3. For quality assessment, review the MultiQC reports to evaluate the sequencing quality before and after trimming.
4. When working with the count data, use **gene_count_matrix.csv** for gene-level analysis or **transcript_count_matrix.csv** for transcript-level analysis.
5. The IIS-specific results focus on immune and inflammatory pathways.

## Data Processing Information

These data files were generated using the scripts available in the Scripts directory:

1. Raw reads were processed using Trimmomatic to remove adapters and low-quality sequences.
2. Clean reads were mapped to the dog reference genome (Boxer Tasha) using HISAT2.
3. Gene and transcript quantification was performed using StringTie.
4. Differential expression analysis was conducted using DESeq2 in R.
5. Pathway analysis was performed using GSEA with KEGG pathways.
6. Visualizations were generated using custom R scripts with ggplot2, pheatmap, and other packages.

## References

For more information on the analysis methods, please refer to:
- Love et al. (2014) - DESeq2 
- Subramanian et al. (2005) - GSEA
- Kim et al. (2019) - HISAT2
- Pertea et al. (2015) - StringTie