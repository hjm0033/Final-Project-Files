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

### 2. `visuals_dog.R`

**Purpose:** Comprehensive visualization of RNA-seq analysis results.

**Functionality:**
- Creates multiple visualization types for differential gene expression data
- Generates PCA plots to visualize sample clustering
- Creates heatmaps of highly variable genes
- Produces MA plots showing differential expression patterns
- Creates volcano plots to visualize significant genes
- Uses advanced visualization libraries including ggplot2, pheatmap, and viridis
- Outputs high-resolution PDF and PNG visualizations

**Input:** DESeq2 results and normalized expression data  
**Output:** Various visualization files including PCA plots, heatmaps, MA plots, and volcano plots

## Usage Instructions

1. Run scripts in numerical order (2_*, 3_*, 4_*)
2. After generating count matrices, proceed with R scripts for analysis:
   - `Dog_DESeq2.R` - Start with differential expression analysis
   - `visuals_dog.R` - Generate various visualizations
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

