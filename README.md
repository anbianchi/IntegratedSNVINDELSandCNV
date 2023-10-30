# Integrated Analysis of Germline SNVs, INDELs, and CNVs in Breast Cancer Whole Exome Sequencing Data

## Table of Contents
- [Description](#description)
  - [Motivation](#motivation)
  - [Results](#results)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Dataset Retrieval](#dataset-retrieval)
  - [Tools & Setup](#tools--setup)
  - [Pipeline Directory Setup](#pipeline-directory-setup)
  - [Pipeline Configuration](#pipeline-configuration)
  - [Installation Options](#installation-options)
- [Dataset Information](#dataset-information)
- [Directory Structure](#directory-structure)

## Description

### Motivation
The rapidly evolving field of genomics emphasizes the need for a holistic understanding of the genetic structures associated with diseases like breast cancer. Current methods often analyze genomic variants such as Single Nucleotide Variants (SNVs), insertions and deletions (INDELs), and Copy Number Variants (CNVs) separately, leaving a gap for integrated analysis.

### Results
We introduce an adaptable method for analyzing SNVs, INDELs, and CNVs from Whole Exome Sequencing (WES) data, emphasizing germline variants. Our approach rigorously validates and filters variants for accuracy. Through the reanalysis of two public WES datasets, our tool highlights its versatility, uncovering both novel and known variations crucial for breast cancer predisposition.

## Installation

### Prerequisites
- Ensure ~200GB of disk space.
- Clone or download this repository.

### Dataset Retrieval
1. [Download the first dataset from ENA browser](https://www.ebi.ac.uk/ena/browser/view/PRJEB3235).
2. For the second dataset, request access.

### Tools & Setup
1. **Annovar**: [Download and configure Annovar](https://annovar.openbioinformatics.org/en/latest/user-guide/download/), following its official documentation. 
2. **GATK Bundle**: [Acquire the GATK bundle (for hg19 genome)](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references).
3. **AnnotSV**: [Install AnnotSV](https://github.com/lgmgeo/AnnotSV), following its official documentation.

### Pipeline Directory Setup
- Use VS Code to navigate to the pipeline directory.
- Follow the provided structure for placing data, configuration files, and tools.

### Pipeline Configuration
1. Choose the appropriate snakefile from the Snakefiles folder.
2. Open Snakefile and set the number of threads for Snakemake.
3. Create config files for datasets based on provided examples.

### Installation Options
1. **Using Docker**:
   - [Download and install Docker](https://www.docker.com/products/docker-desktop/).
   - Configure Docker via the app settings.
   - Build and run the pipeline using the provided commands.
   - Monitor progress and retrieve results within the Docker app.
   - Optionally, execute specific parts of the pipeline as directed.

2. **Using Snakemake (No Docker)**:
   - [Download and install Anaconda](https://docs.anaconda.com/free/anaconda/install/index.html).
   - Set up the environment and install necessary tools using the provided commands.
     ```bash
      conda create -c conda-forge -c bioconda -n #environmentname snakemake -y 
      conda activate approach
      conda install -y c conda-forge perl
      conda install -y -c conda-forge r-base
      conda install -y -c bioconda samtools
      conda install -y -c bioconda trimmomatic
      conda install -y -c conda-forge -c bioconda gatk4
      conda install -y -c bioconda annotsv
      conda install -y -c bioconda bcftools
      R -e "install.packages('BiocManager', repos='http://cran.us.r-project.org')"
      R -e 'BiocManager::install("ExomeDepth")'
      R -e 'BiocManager::install("DNAcopy")'
      R -e 'BiocManager::install("cn.mops")'
      ```
   - Navigate to the appropriate folder inside "replication_cnv_snakemake/Snakefiles" and execute the pipeline with Snakemake.

## Dataset Information
Datasets are publicly sourced. The second dataset includes 7 Illumina samples from [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6813924/). Requests for data can be made to the corresponding author.

### Key Features
- **Accessibility**: The dataset promotes open and collaborative research.

### Accessing the Dataset
- Visit [PRJEB3235](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB3235) and [PRJEB31704](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB31704).
- Download using methods from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/).

> **Note**: Adherence to usage guidelines ensures ethical data usage.

## Directory Structure

### 1. Using Docker
Ensure your directory structure matches this repository to avoid errors.

1. **100bp_exon.bed**: File .bed containing information about genomic regions of interest 
2. **alignedFiles/**: Folder cointaining intermediate file for determining SNV and Indels
3. **AnnotSV/**: 
4. **annovar/**
5. **clean_and_merge.py**: Python file for merging results from different CNV callers
6. **config_paired.csv**: Configuration file for paired_end data
7. **config_single.csv**: Configuration file for single_end data
8. **dockerfile**: docker file containing the instruction to create the docker image
9. **exomeDepth_paired.r**: R script that manages ExomeDepth method for paired_end reads
10. **exomeDepth_single.r**: R script that manages ExomeDepth method for paired_end reads
11. **final**: Folder containing final aligned, deduplicated and recalibrated bam files
12. **gatkbundle/**: Folder containing all the resources (variant files, genome files, ...) required from GATK 
13. **GENDB/**: 
14. **index/**: Folder containing index files for the reference genome
15. **logs/**: Folder containing logs from the execution
16. **mapped/**: Folder containing partial aligned and sorted bam files
17. **reads/**: Folder containing the raw sequencing reads that are to be analyzed in /single/ and /paired/ directories
18. **results/**: Folder containing results from cnv calling
19. **scriptmops.r**: cn.mops script for CNV detection and analysis
20. **Snakefiles/**: Folder containing snakemake's files to run the method
21. **trimmed**: Folder containing trimmed files
22. **tmpgenomicsdb**
23. **tmpPicard**
24. **tables**

### 2. Using Snakemake (No Docker)
Ensure your directory structure matches this repository to avoid errors.

1. **100bp_exon.bed**: File .bed containing information about genomic regions of interest 
2. **alignedFiles/**: Folder cointaining intermediate file for determining SNV and Indels
3. **clean_and_merge.py**: Python file for merging results from different CNV callers
4. **config_paired.csv**: Configuration file for paired_end data
5. **config_single.csv**: Configuration file for single_end data
6. **exomeDepth_paired.r**: R script that manages ExomeDepth method for paired_end reads
7. **exomeDepth_single.r**: R script that manages ExomeDepth method for paired_end reads
8. **final**: Folder containing final aligned, deduplicated and recalibrated bam files
9. **gatkbundle/**: Folder containing all the resources (variant files, genome files, ...) required from GATK 
10. **GENDB/**: 
11. **index/**: Folder containing index files for the reference genome
12. **logs/**: Folder containing logs from the execution
13. **mapped/**: Folder containing partial aligned and sorted bam files
14. **reads/**: Folder containing the raw sequencing reads that are to be analyzed in /single/ and /paired/ directories
15. **results/**: Folder containing results from cnv calling
16. **scriptmops.r**: cn.mops script for CNV detection and analysis
17. **Snakefiles/**: Folder containing snakemake's files to run the method
18. **trimmed**: Folder containing trimmed files
19. **tmpgenomicsdb**
20. **tmpPicard**
21. **tables**
