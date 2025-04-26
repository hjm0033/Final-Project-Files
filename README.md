# Dog Brain RNA-Seq Analysis Scripts

This directory contains scripts for processing and analyzing RNA-seq data from a dog brain study. The workflow includes read trimming, quality control, genome mapping, differential gene expression analysis, and various visualization techniques.

## Workflow Overview

1. **Read Quality Control and Trimming**
2. **Reference Genome Indexing**
3. **Read Mapping to Reference Genome**
4. **Transcript Assembly and Quantification**
5. **Differential Gene Expression Analysis**
6. **Visualization and Interpretation**
7. **Gene Set Enrichment Analysis**

## Shell Scripts

### 1. `2_Trimmomatic_Paired.sh` and `2_Trimmomatic_Unpaired.sh`

**Purpose:** Trim sequencing adapters and low-quality regions from the raw sequence data.

**Functionality:**
- Uses Trimmomatic to remove adapters and low-quality reads
- Processes paired-end and unpaired reads (in respective scripts)
- Performs quality control using FastQC before and after trimming
- Generates MultiQC reports to evaluate overall sequence quality

**Input:** Raw FASTQ read files  
**Output:** Trimmed FASTQ files and quality reports

**ASC Queue Parameters:**
- Queue: medium
- Cores: 6
- Time limit: 02:00:00
- Memory: 12GB

### 2. `3_Index_Ref_Genome.sh`

**Purpose:** Prepare the reference genome for read mapping by creating an index.

**Functionality:**
- Converts annotation file from GFF3 to GTF format using gffread
- Extracts exons and splice sites for improved mapping
- Creates a HISAT2 index for the dog reference genome (Boxer Tasha)

**Input:** Reference genome FASTA and GFF3 annotation file  
**Output:** Indexed genome files for use with HISAT2

**ASC Queue Parameters:**
- Queue: class or medium
- Cores: 6
- Time limit: 04:00:00
- Memory: 12GB

### 3. `4_Map_Reads_To_Genome.sh`

**Purpose:** Map cleaned reads to the reference genome and count reads mapped to features.

**Functionality:**
- Maps reads to the indexed genome using HISAT2
- Converts SAM to BAM format and sorts using SAMtools
- Quantifies gene and transcript expression using StringTie
- Generates count matrices for differential expression analysis

**Input:** Cleaned FASTQ files and indexed reference genome  
**Output:** BAM alignment files, StringTie quantification files, and count matrices

**ASC Queue Parameters:**
- Queue: class or medium
- Cores: 6
- Time limit: 04:00:00
- Memory: 12GB

### 4. `IIS_mapping_consensus_sequence.sh`

**Purpose:** Additional mapping and variant calling pipeline component.

**Functionality:**
- Extends the mapping process with additional options
- Includes variant calling functionality for single nucleotide polymorphism identification 
- Sets up directory structure for variant analysis in IIS (Immune and Inflammatory System) focused analysis

## R Analysis Scripts

### 1. `Dog_DESeq2.R`

**Purpose:** Perform differential gene expression analysis.

**Functionality:**
- Imports count matrices and metadata
- Runs DESeq2 analysis to identify differentially expressed genes
- Creates MA plots to visualize expression differences
- Generates normalized counts for downstream analysis
- Prepares ranked gene lists for Gene Set Enrichment Analysis (GSEA)

**Input:** Count matrices and sample metadata  
**Output:** Differential expression results, visualizations, and prepared data for GSEA

### 2. `Dog_HVGs.R`

**Purpose:** Identify and analyze highly variable genes.

**Functionality:**
- Identifies genes with high variance across samples
- Creates heatmap visualizations using multiple color palettes (viridis)
- Outputs high-resolution PDF and PNG visualizations
- Generates the `highly_variable_genes_heatmap.pdf` in the project root directory

### 3. `Dog_MA.R`

**Purpose:** Create MA plots for visualizing differential expression.

**Functionality:**
- Plots log fold change against mean expression
- Visualizes differentially expressed genes
- Labels top genes with gene symbols for easy identification
- Uses ggrepel for optimal label placement
- Outputs high-resolution PDF visualizations (`ma_plot_dog_size.pdf`)

### 4. `Dog_PCA.R`

**Purpose:** Perform principal component analysis for sample clustering.

**Functionality:**
- Creates PCA plots to visualize sample relationships
- Helps identify patterns and clusters in the data
- Groups samples by dog size (Large vs. Small)

### 5. `Dog_heatmap.R`

**Purpose:** Generate heatmaps for expression data visualization.

**Functionality:**
- Creates heatmaps of gene expression patterns
- Visualizes clustering of samples and genes
- Uses viridis color palette for better visualization 
- Outputs high-resolution PDF visualizations (`sample_distance_heatmap.pdf`)

### 6. `dog_volcano.R`

**Purpose:** Create volcano plots to visualize significance and fold change.

**Functionality:**
- Plots statistical significance against fold change
- Highlights significantly differentially expressed genes
- Creates customized dotplot visualization showing top differentially expressed genes
- Outputs high-resolution PDF visualizations (`top_genes_dotplot.pdf`)

## Other Files

## Project Root Files

Several important output files and metadata files exist in the project's root directory:

### Sample and Run Information
- **SraRunTable_Dogs_RNAseq_All_Possible_Samples.txt**: Complete list of all available SRA samples
- **SraRunTable_Dogs_RNAseq_Best_Samples.txt**: Curated list of best quality samples
- **SraRunTable_Dogs_RNAseq_Proposed_Sample_List.txt**: Final sample list used for analysis

### Analysis Setup Files
- **Adaptor_Sequence_Selection_For_Trimming**: Documentation on adapters used in trimming process

### Visualization Outputs
- **highly_variable_genes_heatmap.pdf**: Heatmap showing the top 50 most variable genes
- **top_genes_dotplot.pdf**: Dotplot visualization of the most significant differentially expressed genes

## Related Data Directory Content

The project includes a Data directory containing various output files from the analysis workflow. Key files include:

### Count Matrices and Expression Data
- **gene_count_matrix.csv**: Gene-level count matrix from StringTie analysis
- **transcript_count_matrix.csv**: Transcript-level count matrix
- **Dog_NormTransExpIDs.txt**: Normalized expression data for downstream analysis

### Differential Expression Results
- **Dog_DGESeq_results.csv**: Complete differential expression results from DESeq2 analysis
- **all_significant_DEGs.txt**: List of all statistically significant differentially expressed genes
- **top25_upregulated.txt** and **top25_downregulated.txt**: Top differentially expressed genes
- **top50_DEGs.txt**: Top 50 differentially expressed genes based on significance
- **top50_variable_genes.txt**: Genes with highest variance across samples

### Gene Set Enrichment Analysis (GSEA) Results
- **Dog_DGErankName.rnk**: Ranked gene list for GSEA input
- **Dog_GSEA_KEGG.GseaPreranked.1743603144299.rpt**: GSEA report using KEGG pathways
- **gsea_report_for_na_pos_1743603144299.html**: HTML report for positively enriched gene sets
- **gsea_report_for_na_neg_1743603144299.html**: HTML report for negatively enriched gene sets
- **gene_set_sizes.tsv**: Information about analyzed gene set sizes
- **ranked_gene_list_na_pos_versus_na_neg_1743603144299.tsv**: Ranked gene list comparing positive vs. negative enrichment
- **pos_snapshot.html** and **neg_snapshot.html**: HTML snapshots of positively and negatively enriched pathways
- **pos_snapshot_IIS.html** and **neg_snapshot_IIS.html**: HTML snapshots of IIS-specific pathway enrichment

### Quality Control Reports
- **multiqc_report_original.html**: Initial MultiQC report of raw sequencing data
- **multiqc_report_reduced_headcrop.html**: Quality report after headcrop trimming
- **multiqc_report_stringent_window.html**: Quality report after stringent window trimming

### Metadata and Statistics
- **Proposed_Samples_Meta_Data.txt**: Metadata for samples used in the analysis
- **overall_statistics.txt**: Summary statistics from the analysis

## Graphs Directory Content

The Graphs directory contains visualization outputs from the GSEA analysis, showing enrichment plots for various KEGG pathways. Key files include:

### GSEA Enrichment Plots
- **enplot_KEGG_ALZHEIMERS_DISEASE_9.png**: Enrichment plot for Alzheimer's disease pathway
- **enplot_KEGG_PARKINSONS_DISEASE_3.png**: Enrichment plot for Parkinson's disease pathway
- **enplot_KEGG_HUNTINGTONS_DISEASE_6.png**: Enrichment plot for Huntington's disease pathway

### Neurological Pathway Visualizations
- **enplot_KEGG_LONG_TERM_POTENTIATION_11.png**: Enrichment plot for long-term potentiation
- **enplot_KEGG_AXON_GUIDANCE_31.png**: Enrichment plot for axon guidance

### Metabolic and Cellular Process Visualizations
- **enplot_KEGG_OXIDATIVE_PHOSPHORYLATION_2.png**: Enrichment plot for oxidative phosphorylation
- **enplot_KEGG_CITRATE_CYCLE_TCA_CYCLE_29.png**: Enrichment plot for citric acid cycle
- **enplot_KEGG_RIBOSOME_1.png**: Enrichment plot for ribosomal functions
- **enplot_KEGG_SPLICEOSOME_4.png**: Enrichment plot for splicing machinery

### Immune-Related Pathway Visualizations
- **enplot_KEGG_COMPLEMENT_AND_COAGULATION_CASCADES_59.png**: Enrichment plot for complement cascade
- **enplot_KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_77.png**: Enrichment plot for cytokine interactions
- **enplot_from_IIS_*.png**: Specific enrichment plots for immune and inflammatory system pathways

## Usage Instructions

1. Run scripts in numerical order (2_*, 3_*, 4_*)
2. After generating count matrices, proceed with R scripts for analysis:
   - `Dog_DESeq2.R` - Start with differential expression analysis
   - `Dog_HVGs.R`, `Dog_MA.R`, `Dog_PCA.R`, etc. - Generate visualizations
3. Review GSEA results in the Data and Graphs directories

## Dependencies

- **Bioinformatics Tools:** Trimmomatic, FastQC, MultiQC, HISAT2, StringTie, SAMtools, BCFtools
- **R Packages:** 
  - **Core Analysis:** DESeq2, edgeR
  - **Visualization:** ggplot2, pheatmap, RColorBrewer, viridis, ggrepel
  - **Statistics:** stats, BiocManager
- **Other:** Python (for prepDE.py)

## Notes

- These scripts were designed to run on the Alabama Super Computer (ASC) and may need adjustments for other systems
- Parameters and file paths should be modified according to your specific project setup
- Please ensure all required modules and dependencies are installed before running these scripts
- Results files in the Data directory are produced by executing the scripts in this directory
- The Graphs directory contains visualization outputs from GSEA and other analyses
- This project focuses on comparing gene expression differences between large and small dog brains
- Key findings include enrichment of neurological disease pathways and metabolic functions in the dataset

