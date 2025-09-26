# APIARC: Integrated Analysis Workflow for RNA-seq and ChIP-seq Data

## Overview

APIARC (Analysis Pipeline for Integrated RNA-seq and ChIP-seq) is a comprehensive computational workflow designed to identify and analyze common genes from both RNA-seq and ChIP-seq data, enabling the discovery of novel transcriptional regulatory mechanisms.

## Installation

### 1. Install Mamba and Snakemake

Follow these steps to set up the environment:

1. **Access your Linux terminal**
2. **Initialize the base environment:**
   ```bash
   source ~/.bashrc
   ```
3. **Install Mamba:**
   ```bash
   conda install -n base -c conda-forge mamba
   ```
4. **Create a new environment and install Snakemake:**
   ```bash
   mamba create -c conda-forge -c bioconda -n APIARC snakemake
   ```

### 2. Install Git and Clone Repository

1. **Check Git installation:**
   ```bash
   git --version
   ```
   *Note: If Git is not installed, download from [https://git-scm.com/downloads](https://git-scm.com/downloads)*

2. **Clone the APIARC repository:**
   ```bash
   git clone https://github.com/yxiaobo/APIARC.git
   ```

### 3. Activate and Update Environment

1. **Activate the environment:**
   ```bash
   conda activate APIARC;
   cd APIARC
   ```

2. **Install required dependencies:**
   ```bash
   mamba env update -n APIARC -f env/environment.yml
   ```

**Note:** For future sessions, activate the environment with:
```bash
source ~/.bashrc;
conda activate APIARC
```

---

### 2. Deploy and Configure Workflows  

APIARC consists of the following workflows:  

#### Reference  
[reference](./reference): Downloads reference data required for all workflows.  
   ```bash
   cd ~/PROJECT/APIARC/reference
   wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz -P ./genome
   wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz -P ./genome
   gunzip ./genome/GRCh38.p14.genome.fa.gz
   gunzip ./genome/GRCm38.p6.genome.fa.gz
   bowtie2-build --threads 5 ./genome/GRCh38.p14.genome.fa ./bowtie2_index/humanbowtie2Index/GRCh38
   bowtie2-build --threads 5 ./genome/GRCm38.p6.genome.fa ./bowtie2_index/mousebowtie2Index/GRCm38
   hisat2-build -p 5 ./genome/GRCh38.p14.genome.fa ./hisat2_index/humanhisat2Index/GRCh38
   hisat2-build -p 5 ./genome/GRCm38.p6.genome.fa ./hisat2_index/mousehisat2Index/GRCm38
   wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz -P ./annotation
   wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz -P ./annotation
   gunzip ./annotation/gencode.v44.annotation.gtf.gz
   gunzip ./annotation/gencode.vM25.annotation.gtf.gz
   ```
#### RNA-seq Workflow  
[RNAPROJ](./RNAPROJ): Processes RNA sequencing data.  

#### ChIP-seq Workflow  
[CHIPPROJ](./CHIPPROJ): Processes Chromatin Immunoprecipitation Sequencing data.  

#### Integrated Analysis Workflow  
[INTEGRATEDPROJ](./INTEGRATEDPROJ): Analyzes common genes identified by `RNAPROJ` and `CHIPPROJ`.  

**Important:**  
- Complete `RNAPROJ` and `CHIPPROJ` before running `INTEGRATEDPROJ`.  
- The workflow concludes after `INTEGRATEDPROJ` finishes execution.  

---

## System Requirements  

### Recommended Specifications  
- **Operating System:** Linux (tested on Ubuntu 20.04)  
- **Memory:** ≥ 16 GB RAM  
- **CPU:** ≥ 8 cores (10 cores tested)  
- **Storage:** ≥ 500 GB available space  

### Software Versions (Tested)  
- `mamba=2.0.8`  
- `snakemake=7.32.4`  

---

## Reproducibility  
APIARC is designed for Linux systems. Ensure all dependencies are installed as specified to guarantee reproducibility.  

For optimal performance, adhere to the recommended system specifications.
