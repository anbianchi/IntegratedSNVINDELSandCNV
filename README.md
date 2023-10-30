# Replication Package - An integrated approach to discover germline SNVs, INDELs and CNVs from whole exome sequencing data of patients affected with breast cancer

## Description


### Motivation:
In the swiftly progressing domain of genomics, a profound understanding of the genetic frameworks underpinning complex diseases like breast cancer (BC) is imperative. This necessitates a comprehensive analysis encompassing various genomic variants such as Single Nucleotide Variants (SNVs), insertions and deletions (indels), and Copy Number Variants (CNVs). The existing methodologies primarily tackle these variants in isolation, thereby lacking a unified platform for an integrative analysis
### Results: 
We present a powerful, reproducible, and adaptable method for a comprehensive analysis of SNVs, indels, and CNVs produced from Whole Exome Sequencing (WES) data, with an emphasis on germinal variants, in order to fill this crucial gap. Our creative approach combines a rigorous validation process with a variant filtering strategy to accurately find variants. Re-analyzing two publicly available WES datasets allowed us to show the tool's adaptability by revealing a wide range of novel and well-known variations that are important for BC susceptibility. This method demonstrated accuracy as well as the capacity to include a broad range of genomic variations, demonstrating its capacity for a single, integrated study. Our tool is a ground-breaking solution in this field because to the diligent handling of technological issues. 
### Availability: 
The code is freely available for non-commercial users and can be accessed on the web at https://github.com/alessandrodimatteo97/BioInf-Replication- \textcolor{red}{prima di sottomettere, version finale e renderlo free il progetto github}.
The pipeline can be run using docker or without docker, just using anaconda.

## Installation

### Prerequisites

Ensure you have approximately 200GB of disk space available before proceeding with the following steps to run the pipeline:

1. **Clone or Download the Pipeline Code:**
   - Clone this repository to your local machine or download the entire code.

2. **Dataset Retrieval:**
   - Download the first dataset from [ENA browser](https://www.ebi.ac.uk/ena/browser/view/PRJEB3235).
   - For the second dataset, request access.

3. **Reference File Acquisition:**
   - Obtain the `hg19.fa`,`hg19.dict` and `hg19.fa.fai` reference file from ....

4. **Annovar Installation:**
   - Download **Annovar** from the official website (https://annovar.openbioinformatics.org/en/latest/user-guide/download/).
   - Download hg19_ALL.sites.2015_08.txt from http://www.openbioinformatics.org/annovar/download/hg19_1000g2015aug.zip and extract   
    **hg19_ALL.sites.2015_08.txt.idx**, **hg19_ALL.sites.2015_08.txt**
   - Download **hg19_clinvar_20190305.txt** (http://www.openbioinformatics.org/annovar/download/hg19_clinvar_20190305.txt.gz) 
   - Download **hg19_clinvar_20190305.txt.idx** (https://www.openbioinformatics.org/annovar/download/hg19_clinvar_20190305.txt.idx.gz)
   - Download **hg19_exac03.txt** (https://www.openbioinformatics.org/annovar/download/hg19_exac03.txt.gz)
   - Download **hg19_exac03.txt.idx** (https://www.openbioinformatics.org/annovar/download/hg19_exac03.txt.idx.gz)   
   - Download **hg19_gnomad211_exome.txt** (https://www.openbioinformatics.org/annovar/download/hg19_gnomad211_exome.txt.gz)
   - Download **hg19_gnomad211_exome.txt.idx** (https://www.openbioinformatics.org/annovar/download/hg19_gnomad211_exome.txt.idx.gz)
   - Download **hg19_ljb26_all.txt** (http://www.openbioinformatics.org/annovar/download/hg19_ljb26_all.txt.gz)
   - Download **hg19_ljb26_all.txt.idx** (http://www.openbioinformatics.org/annovar/download/hg19_ljb26_all.txt.idx.gz)
   - Download **hg19_popfreq_max_20150413.txt** (https://www.openbioinformatics.org/annovar/download/hg19_popfreq_max_20150413.txt.gz)
   - Downlaod **hg19_popfreq_max_20150413.txt.idx** (https://www.openbioinformatics.org/annovar/download/hg19_popfreq_max_20150413.txt.idx.gz)
   - Put all of these file inside Annovar/humandb


5. **GATK Bundle:**
   - Acquire the GATK bundle from the specified source (https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references).

6. **IDE Selection:**
   - [Download](https://code.visualstudio.com/download) and install Visual Studio Code (VS Code) or use an IDE of your choice.


7. **Pipeline Directory Setup:**
   - Using VS Code, navigate to the pipeline directory.
   - Place unpaired data in `reads/single/` and paired data in `reads/paired/` within the pipeline directory.
   - Place the `hg19.fa` and `hg19.dict` reference file in `index/` directory (Update rule `bwa index`, `bwa mem`, and `bwa mem paired` in the Snakefile if reference file name changes).
   - Place the `gatkbundle` directory in `gatkbundle/hg19/`, ensuring it contains the necessary files as listed in the pipeline documentation.
   - Place the Annovar directory within the pipeline directory.

8. **Pipeline Configuration:**
    - Select from Snakefiles folder the right snakefile to run and move to the pipeline folder
    - Open Snakefile and specify the number of threads for Snakemake to use.
    - Create config files for your paired and unpaired datasets, following the examples provided in `config_paired.csv` and `config_single_csv` by specifing the complete name of your reads files.


### 1. Installation using Docker


7. **Docker Installation:**
   - [Download](https://www.docker.com/products/docker-desktop/) and install Docker on your machine.


8. **Docker Configuration:**
   - Launch Docker App.
   - Navigate to Settings → Resources to allocate the desired amount of CPU, Memory, Swap, and Virtual Disk Limit for Docker containers.


11. **Pipeline Execution:**
    - Open terminal in VS Code (View → Terminal).
    - Execute the following commands to build and run the pipeline:
        ```bash
        docker build -t bioinfpipeline .
        docker run --name pipelinerun -it bioinfpipeline
        ```
    - Monitor the pipeline's progress within the Docker app by selecting containers → `pipelinerun`. View file outputs in the file section and running commands in the logs section.
    - Wait for the execution to complete.

12. **Result Retrieval:**
    - Once execution is complete, download and analyze the necessary files from the file window in the Docker app.

13. **Partial Pipeline Execution (Optional):**
    - If needed, execute specific parts of the pipeline changing the current snakefile with the respective Snakefiles located in `SNKFL/SNKFR`, `SNKFRPaired`, `SNKFRUnpaired` directories by moving the file to the main directory.



### 2. Installation withouth Docker


7. **Anaconda Installation:**
   - [Download](https://docs.anaconda.com/free/anaconda/install/index.html) and install Anaconda on your machine (depending on your OS following the reported instruction)

  
8. **Environment Setup:**
   - Open terminal in VS Code (View → Terminal).
   - Execute the following commands to set the conda environment "approach"
     ```bash
      conda create -c conda-forge -c bioconda -n approach snakemake -y 
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


9.  **Pipeline Execution:**
    - Execute the following command to build and run the pipeline:
        ```bash
        snakemake --use-conda --cores
        ```
    - Monitor the pipeline's progress within the terminal

10. **Result Retrieval:**
    - Once execution is complete, download and analyze the necessary files from the directory

11. **Partial Pipeline Execution (Optional):**
    - If needed, execute specific parts of the pipeline using the respective Snakefiles located in `SNKFL/SNKFR`, `SNKFRPaired`, `SNKFRUnpaired` directories.



### Dataset Information

The datasets were derived from sources in the public domain: [(https://www.ncbi.nlm.nih.gov/bioproject/PRJEB3235) and https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB31704]. Regarding the second dataset, we employed 7 Illumina samples. Those were provided by the authors of [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6813924/). Data will be shared on request to the corresponding author.

#### Key Features:

- **Accessibility**: The dataset is publicly available to researchers across the world, fostering a collaborative and open research environment.


#### Accessing the Dataset:
To access and use the datasets:

1. **SRA Access**: Visit the [PRJEB3235](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB3235) and [PRJEB31704](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB31704).

2. **Downloading the Data**: Download the datasets using one of the methods in https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/. 

> Note: Even though the dataset is publicly available, we strictly adhere to the usage guidelines, ensuring ethical use of the data in our research.


### Code Structure

#### <span style='color: red;'>Please be certain that your directory structure mirrors that of this repository; failure to do so may lead to potential errors.
</span>

1. **100bp_exon.bed**: File .bed containing information about genomic regions of interest 
2. **alignedFiles/**:
3. **alignedFiles/sortedBam/**:
4. **alignedFiles/sortedBam/gvcf**:
5. **AnnotSV/**: 
6. **annovar/**
7. **clean_and_merge.py**:
8. **config_paired.csv**:
9. **config_single.csv**:
10. **dockerfile**: docker file containing the instruction to create the docker image
11. **exomeDepth_paired.r**:
12. **exomeDepth_single.r**:
13. **final**:
14. **gatkbundle/hg19**: Folder containing all the resources (variant files, genome files, ...) required from GATK 
15. **GENDB/**: 
16. **GENDB/annovar/**:
17. **index/**: Folder containing index files for the reference genome
18. **logs/**
19. **mapped/**: it stores mapped files, post-alignment
20. **reads/**: Folder containing the sequencing reads that are to be analyzed in /single/ and /paired/ directories
21. **reads/single/**: Folder containing the sequencing single reads 
22. **reads/paired/**: Folder containing the sequencing paired reads
23. **results/**:
24. **scriptmops.r**: cn.mops script for CNV detection and analysis
25. **Snakefiles/**: directory containig other snakemake file to run the right pipeline
26. **exomeDepth.r**: ExomeDepth script for analyzing and identifying CNV
27. **trimmed**:
28. 


