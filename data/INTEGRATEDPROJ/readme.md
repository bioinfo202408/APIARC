# INTEGRATEDPROJ Analysis Pipeline


## Overview

INTEGRATEDPROJ is a highly automated integrated analysis pipeline for RNA-seq and ChIP-seq data, featuring two automated workflows: **Promoter** and **Enhancer** analysis. 

The **Promoter analysis workflow** focuses on studying the relationship between transcription factor (TF) binding and gene expression regulation, implementing promoter region analysis, motif and transcription factor enrichment analysis, multi-omics data integration and visualization, and functional annotation and pathway analysis.

The **Enhancer analysis workflow** investigates enhancer regulatory mechanisms, performing enhancer-gene association analysis, motif and transcription factor enrichment analysis, multi-omics data integration and visualization, and functional annotation and pathway analysis.

Utilizing the Snakemake workflow management system, this pipeline ensures reproducibility and efficiency throughout the entire analysis process.

> **Important**: Ensure you are in the APIARC directory with the appropriate environment (e.g., APIARC) activated. All analyses should be performed from this directory without changing paths. Complete the RNAPROJ and INTEGRATEDPROJ analysis workflows first to ensure input files are available.

## Promoter Analysis


## Quick Start Guide

### 1. Workflow Configuration

Modify parameters in [`data/INTEGRATEDPROJ/script/Promoter/config.yaml`](data/INTEGRATEDPROJ/script/Promoter/config.yaml) as follows:

#### System Configuration
- **`conda_install_path`**: Path to your Conda installation (e.g., `/home/user/miniconda3`)

#### File Group Parameters
*(These parameter settings should be consistent with CHIPPROJ's [`config.yaml`](data/CHIPPROJ/config.yaml))*

> **Naming Convention**: Biological replicates must follow the `groupname_repX` pattern, where:
> - `groupname`: Your custom identifier for the sample group
> - `X`: Replicate number (must be sequential starting from 1)
> - Example: For three replicates: `groupname_rep1`, `groupname_rep2`, `groupname_rep3`

**Required parameters for each sample group:**
- **`treat_bam`**: Treatment sample SRR ID
- **`control_bam`**: Control sample SRR ID  
- **`genome`**: Reference genome (`mm` for mouse, `hg` for human)
- **`qval`**: FDR cutoff for peak calling (default: 0.05)

**Configuration Example:**
```yaml
Macs2:
  mESC_18h_DOX_rep1:
    treat_bam: SRR5297329
    control_bam: SRR5297331
    genome: mm
    qval: 0.05

  mESC_18h_DOX_rep2:
    treat_bam: SRR5297330
    control_bam: SRR5297331
    genome: mm
    qval: 0.05

  mESC_18h_DOX_rep3:
    treat_bam: SRR5297332
    control_bam: SRR5297331
    genome: mm
    qval: 0.05
```

#### Experimental Group Assignment

> **Important**: Group names must match those defined in the Macs2 section (without `_repX` suffix) and correctly assign samples to experimental conditions.

**Group designation:**
- **`P`**: Treatment group
- **`T`**: Control group

**Configuration Example:**
```yaml
groupnames:
  mESC_18h_DOX: P    # Treatment group
  mESC_12h_DOX: T    # Control group

# Group list (note: space required after dash)
groups:  
  - mESC_18h_DOX
  - mESC_12h_DOX
```

---

### 2. Workflow Execution

1. **Test configuration (dry run):**
   ```bash
   snakemake --snakefile data/INTEGRATEDPROJ/script/Promoter/snakefile --configfile data/INTEGRATEDPROJ/script/Promoter/config.yaml -n -p
   ```

2. **Execute pipeline:**
   ```bash
   nohup snakemake --snakefile data/INTEGRATEDPROJ/script/Promoter/snakefile --configfile data/INTEGRATEDPROJ/script/Promoter/config.yaml -j <N> > log_Promoter.txt 2>&1 &
   ```
   Where `<N>` = number of CPU threads to use

---

## Configuration Template

### [`config.yaml`](data/INTEGRATEDPROJ/script/Promoter/config.yaml)
```yaml
# Your conda.sh path: {conda_install_path}/etc/profile.d/conda.sh
conda_install_path: /home/yxiaobo/miniconda3    ## Enter the prefix path of your Miniconda installation {conda_install_path}

## Note: Ensure {groupname_repx}, {treat_bam}, {control_bam}, {genome}, {qval}, {groupname} 
## match the corresponding parameters in CHIPPROJ's config.yaml (data/CHIPPROJ/config.yaml)
file:
  mESC_18h_DOX_rep1:        ## Name each of your experimental and control sample groups {groupname_repx}
    treat_bam: SRR5297329   ## Enter your treatment group SRR ID
    control_bam: SRR5297331 ## Enter your control group SRR ID
    genome: mm              ## Select your species: mm (mouse) or hg (human)
    qval: 0.05              ## Select the p-value threshold you want to use

  mESC_18h_DOX_rep2:    
    treat_bam: SRR5297330
    control_bam: SRR5297331
    genome: mm
    qval: 0.05
  
  mESC_12h_DOX_rep1:
    treat_bam: SRR5297327
    control_bam: SRR5297331
    genome: mm
    qval: 0.05

  mESC_12h_DOX_rep2:
    treat_bam: SRR5297328
    control_bam: SRR5297331
    genome: mm
    qval: 0.05

# Sample grouping - Remove repx and assign samples to experimental and control groups correctly
groupnames:      
  mESC_18h_DOX: P    ## P for treatment group
  mESC_12h_DOX: T    ## T for control group

# Group list corresponding to groupnames configuration (Note: a space is required after the dash)
groups:  
  - mESC_18h_DOX 
  - mESC_12h_DOX
```

---

### Key Output Files

#### GO Functional Enrichment Results
- **`GO_enrich_common_gene.csv`**:
  ```
  "ONTOLOGY","ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalue","rank","leading_edge","core_enrichment"
  "MF","GO:0019783","ubiquitin-like protein peptidase activity",11,0.869409068594305,1.49004920675373,0.00100275040686698,0.37212098155533,0.37212098155533,25,"tags=45%, list=3%, signal=45%","Usp17lc/Usp17ld/Usp17lb/Usp17le/Usp1
  ```

#### KEGG Pathway Enrichment Results
- **`KEGG_GSEA_common_gene.csv`**:
  ```
  "ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalue","rank","leading_edge","core_enrichment"
  "mmu05417","Lipid and atherosclerosis - Mus musculus (house mouse)",15,0.786218588344221,1.40881509861291,0.0108990966102584,0.27247741525646,0.263872865300993,64,"tags=27%, list=8%, signal=25%","Rxrg/Hspa1a/Hsp
  ```

#### Visualization Outputs
- **`common_up_genes_venn_diagram.pdf` / `common_down_genes_venn_diagram.pdf`**:  
  Overlap between RNA-seq (up/down-regulated) and promoter-bound ChIP-seq genes.
- **`common_genes_volcano_plot.pdf`**:  
  Volcano plots display log2FC and -log10(padj) for common genes.
- **`TF_up_gene_heatmap.pdf` / `TF_down_gene_heatmap.pdf`**:  
  Heatmap showing expression patterns of upregulated/downregulated transcription factors with hierarchical clustering.

---

| File Name | Directory | Description | Tool/Source |
| :--- | :--- | :--- | :--- |
| common_up_genes_venn_diagram.pdf | data/INTEGRATEDPROJ/Promoter | Venn diagram of common up-regulated genes | VennDiagram (R) |
| common_down_genes_venn_diagram.pdf | data/INTEGRATEDPROJ/Promoter | Venn diagram of common down-regulated genes | VennDiagram (R) |
| common_genes_volcano_plot.pdf | data/INTEGRATEDPROJ/Promoter | Volcano plot of common DEGs | ggplot2 (R) + ggrepel |
| RNA_ChIP_col_plot_same_direction.pdf | data/INTEGRATEDPROJ/Promoter | Correlation plot between RNA-seq and ChIP-seq log2FC | ggplot2 (R) |
| {groupname}_up_peaks.bed | data/INTEGRATEDPROJ/Promoter | BED file of up-regulated peaks | Promoter_plot.R |
| {groupname}_down_peaks.bed | data/INTEGRATEDPROJ/Promoter | BED file of down-regulated peaks | Promoter_plot.R |
| {groupname}_all_peaks.bed | data/INTEGRATEDPROJ/Promoter | BED file of all peaks | Promoter_plot.R |
| GO_enrich_common_gene.csv | data/INTEGRATEDPROJ/Promoter | GO enrichment results for common genes | clusterProfiler (R) |
| KEGG_GSEA_common_gene.csv | data/INTEGRATEDPROJ/Promoter | KEGG GSEA results for common genes | clusterProfiler (R) |
| Promoter_both_up_genes_target_seqs.fasta | data/INTEGRATEDPROJ/Promoter | Extracted promoter sequences for up-regulated genes | Custom Perl script |
| Promoter_both_down_genes_target_seqs.fasta | data/INTEGRATEDPROJ/Promoter | Extracted promoter sequences for down-regulated genes | Custom Perl script |
| Motif_result_up | data/INTEGRATEDPROJ/Promoter | Directory containing AME motif enrichment results for up-regulated genes | AME (MEME Suite) |
| Motif_result_down | data/INTEGRATEDPROJ/Promoter | Directory containing AME motif enrichment results for down-regulated genes | AME (MEME Suite) |
| TF_motif_network_up.txt | data/INTEGRATEDPROJ/Promoter | TF-motif network with significance for up-regulated genes | Custom Perl script |
| TF_motif_network_down.txt | data/INTEGRATEDPROJ/Promoter | TF-motif network with significance for down-regulated genes | Custom Perl script |
| TF_gene_mapping_up.txt | data/INTEGRATEDPROJ/Promoter | TF-to-gene mapping table for up-regulated genes | Custom Perl script |
| TF_gene_mapping_down.txt | data/INTEGRATEDPROJ/Promoter | TF-to-gene mapping table for down-regulated genes | Custom Perl script |
| Promoter_both_up_genes_matrix.gz | data/INTEGRATEDPROJ/Promoter | Compressed matrix for profile plotting of up-regulated genes | computeMatrix (deepTools) |
| Promoter_both_up_genes_signal.tab | data/INTEGRATEDPROJ/Promoter | Signal values for up-regulated regions | computeMatrix (deepTools) |
| Promoter_both_up_genes_plotProfile.pdf | data/INTEGRATEDPROJ/Promoter | Profile plot of signal around up-regulated peaks | plotProfile (deepTools) |
| Promoter_both_down_genes_matrix.gz | data/INTEGRATEDPROJ/Promoter | Compressed matrix for profile plotting of down-regulated genes | computeMatrix (deepTools) |
| Promoter_both_down_genes_signal.tab | data/INTEGRATEDPROJ/Promoter | Signal values for down-regulated regions | computeMatrix (deepTools) |
| Promoter_both_down_genes_plotProfile.pdf | data/INTEGRATEDPROJ/Promoter | Profile plot of signal around down-regulated peaks | plotProfile (deepTools) |
| Motif_result_up_list.txt | data/INTEGRATEDPROJ/Promoter | List of enriched motifs with TF mappings for up-regulated genes | extr_enriched_motifs.pl |
| Motif_result_down_list.txt | data/INTEGRATEDPROJ/Promoter | List of enriched motifs with TF mappings for down-regulated genes | extr_enriched_motifs.pl |
| TF_up_gene_heatmap.pdf | data/INTEGRATEDPROJ/Promoter | Heatmap of upregulated TFs | ComplexHeatmap (R) |
| TF_down_gene_heatmap.pdf | data/INTEGRATEDPROJ/Promoter | Heatmap of downregulated TFs | ComplexHeatmap (R) |

---

| Tool | Purpose | Documentation |
| :--- | :--- | :--- |
| **VennDiagram (R)** | Visualize overlaps between RNA-seq and ChIP-seq gene sets | [Link](https://cran.r-project.org/web/packages/VennDiagram/index.html) |
| **ggplot2 (R)** | Create publication-quality plots (volcano plots, correlation plots) | [Link](https://ggplot2.tidyverse.org/) |
| **ggrepel** | Prevent label overlapping in ggplot2 visualizations | [Link](https://ggrepel.slowkow.com/) |
| **Promoter_plot.R** | Integrate RNA-seq and ChIP-seq data, perform comparative analysis and visualization | [Link](data/INTEGRATEDPROJ/script/Promoter/Promoter_plot.R) |
| **MEME Suite (AME)** | Motif enrichment analysis in promoter regions | [Link](https://meme-suite.org/meme/) |
| **deeptools** | Process BigWig files and generate signal profiles around promoter peaks | [Link](https://deeptools.readthedocs.io/) |
| **TF_motif_enrichment_pipeline.pl** | Pipeline for TF motif discovery, network building, and TF-gene mapping | [Link](data/INTEGRATEDPROJ/script/Promoter/TF_motif_enrichment_pipeline.pl) |
| **extr_enriched_motifs.pl** | Extract significantly enriched motifs from AME results | [Link](data/INTEGRATEDPROJ/script/Promoter/extr_enriched_motifs.pl) |
| **ComplexHeatmap** | Create heatmaps for TF expression patterns | [Link](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html) |
| **TF_motif_plot.R** | Generate heatmaps for transcription factor expression analysis | [Link](data/INTEGRATEDPROJ/script/Promoter/TF_motif_plot.R) |
| **snakefile** | Promoter analysis workflow automation | [Link](data/INTEGRATEDPROJ/script/Promoter/snakefile) |

---

## Enhancer Analysis


## Quick Start Guide

### 1. Workflow Configuration

Modify parameters in [`data/INTEGRATEDPROJ/script/Enhancer/config.yaml`](data/INTEGRATEDPROJ/script/Enhancer/config.yaml) as follows:

#### System Configuration
- **`conda_install_path`**: Path to your Conda installation (e.g., `/home/user/miniconda3`)

#### File Group Parameters
*(These parameter settings should be consistent with CHIPPROJ's [`config.yaml`](data/CHIPPROJ/config.yaml))*

> **Naming Convention**: Biological replicates must follow the `groupname_repX` pattern, where:
> - `groupname`: Your custom identifier for the sample group
> - `X`: Replicate number (must be sequential starting from 1)
> - Example: For three replicates: `groupname_rep1`, `groupname_rep2`, `groupname_rep3`

**Required parameters for each sample group:**
- **`treat_bam`**: Treatment sample SRR ID
- **`control_bam`**: Control sample SRR ID  
- **`genome`**: Reference genome (`mm` for mouse, `hg` for human)
- **`qval`**: FDR cutoff for peak calling (default: 0.05)

**Configuration Example:**
```yaml
Macs2:
  mESC_18h_DOX_rep1:
    treat_bam: SRR5297329
    control_bam: SRR5297331
    genome: mm
    qval: 0.05

  mESC_18h_DOX_rep2:
    treat_bam: SRR5297330
    control_bam: SRR5297331
    genome: mm
    qval: 0.05

  mESC_18h_DOX_rep3:
    treat_bam: SRR5297332
    control_bam: SRR5297331
    genome: mm
    qval: 0.05
```

#### Experimental Group Assignment

> **Important**: Group names must match those defined in the Macs2 section (without `_repX` suffix) and correctly assign samples to experimental conditions.

**Group designation:**
- **`P`**: Treatment group
- **`T`**: Control group

**Configuration Example:**
```yaml
groupnames:
  mESC_18h_DOX: P    # Treatment group
  mESC_12h_DOX: T    # Control group

# Group list (note: space required after dash)
groups:  
  - mESC_18h_DOX
  - mESC_12h_DOX
```

---

### 2. Workflow Execution

1. **Test configuration (dry run):**
   ```bash
   snakemake --snakefile data/INTEGRATEDPROJ/script/Enhancer/snakefile --configfile data/INTEGRATEDPROJ/script/Enhancer/config.yaml -n -p
   ```

2. **Execute pipeline:**
   ```bash
   nohup snakemake --snakefile data/INTEGRATEDPROJ/script/Enhancer/snakefile --configfile data/INTEGRATEDPROJ/script/Enhancer/config.yaml -j <N> > log_Enhancer.txt 2>&1 &
   ```
   Where `<N>` = number of CPU threads to use

---

## Configuration Template

### [`config.yaml`](data/INTEGRATEDPROJ/script/Enhancer/config.yaml)
```yaml
# Your conda.sh path: {conda_install_path}/etc/profile.d/conda.sh
conda_install_path: /home/yxiaobo/miniconda3    ## Enter the prefix path of your Miniconda installation {conda_install_path}

## Note: Ensure {groupname_repx}, {treat_bam}, {control_bam}, {genome}, {qval}, {groupname} 
## match the corresponding parameters in CHIPPROJ's config.yaml (data/CHIPPROJ/config.yaml)
file:
  mESC_18h_DOX_rep1:        ## Name each of your experimental and control sample groups {groupname_repx}
    treat_bam: SRR5297329   ## Enter your treatment group SRR ID
    control_bam: SRR5297331 ## Enter your control group SRR ID
    genome: mm              ## Select your species: mm (mouse) or hg (human)
    qval: 0.05              ## Select the p-value threshold you want to use

  mESC_18h_DOX_rep2:    
    treat_bam: SRR5297330
    control_bam: SRR5297331
    genome: mm
    qval: 0.05
  
  mESC_12h_DOX_rep1:
    treat_bam: SRR5297327
    control_bam: SRR5297331
    genome: mm
    qval: 0.05

  mESC_12h_DOX_rep2:
    treat_bam: SRR5297328
    control_bam: SRR5297331
    genome: mm
    qval: 0.05

# Sample grouping - Remove repx and assign samples to experimental and control groups correctly
groupnames:      
  mESC_18h_DOX: P    ## P for treatment group
  mESC_12h_DOX: T    ## T for control group

# Group list corresponding to groupnames configuration (Note: a space is required after the dash)
groups:  
  - mESC_18h_DOX 
  - mESC_12h_DOX
```

---

### Key Output Files

#### GO Functional Enrichment Results
- **`GO_enrich_common_gene.csv`**:
  ```
  "ONTOLOGY","ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalue","rank","leading_edge","core_enrichment"
  "MF","GO:0031625","ubiquitin protein ligase binding",49,0.774426547805831,1.58771070580483,2.98242724643566e-06,0.00694010820245578,0.00694010820245578,155,"tags=27%, list=10%, signal=25%","Spopfm2/Tdpoz4/Tdpoz1/Tdpoz5/Tdpoz2/Hspa1a/Tdpoz8/Tdpoz3/Hspa1b/Tdpoz9/Sh3kbp1/Jun/Magea10"
  ```

#### Functional Enhancer-Gene Pairs
- **`Ture_enh_gene_within_1mb.bed`**:
  ```
  chrX	74865201	74865400	ENH504903	Ctag2l1
  chrX	74865201	74865400	ENH504903	Or13ae2
  chrX	74865201	74865400	ENH504903	Rab39b
  chr5	140556801	140557200	ENH168484	Elfn1
  ```

#### Visualization Outputs
- **`RNA_ChIP_col_plot_same_direction.pdf`**:  
  Correlation scatter plot demonstrating coordinated changes between RNA expression and ChIP-seq binding signals.
- **`Enhancer_both_up/down_genes_plotProfile.pdf`**:  
  Signal Profile - Histone modification patterns around enhancer regions.
- **`TF_up/down_gene_heatmap.pdf`**:  
  Heatmap showing expression patterns of upregulated/downregulated transcription factors with hierarchical clustering.

---

| File Name | Directory | Description | Tool/Source |
| :--- | :--- | :--- | :--- |
| Enhancer_both_up_genes_target_seqs.fasta | data/INTEGRATEDPROJ/Enhancer | Extracted enhancer sequences for up-regulated genes | TF_motif_enrichment_pipeline.pl |
| Enhancer_both_down_genes_target_seqs.fasta | data/INTEGRATEDPROJ/Enhancer | Extracted enhancer sequences for down-regulated genes | TF_motif_enrichment_pipeline.pl |
| Motif_result_up | data/INTEGRATEDPROJ/Enhancer | Directory containing AME motif enrichment results for up-regulated genes | AME (MEME Suite) |
| Motif_result_down | data/INTEGRATEDPROJ/Enhancer | Directory containing AME motif enrichment results for down-regulated genes | AME (MEME Suite) |
| TF_motif_network_up.txt | data/INTEGRATEDPROJ/Enhancer | TF-motif network with significance for up-regulated genes | TF_motif_enrichment_pipeline.pl |
| TF_motif_network_down.txt | data/INTEGRATEDPROJ/Enhancer | TF-motif network with significance for down-regulated genes | TF_motif_enrichment_pipeline.pl |
| TF_gene_mapping_up.txt | data/INTEGRATEDPROJ/Enhancer | TF-to-gene mapping table for up-regulated genes | TF_motif_enrichment_pipeline.pl |
| TF_gene_mapping_down.txt | data/INTEGRATEDPROJ/Enhancer | TF-to-gene mapping table for down-regulated genes | TF_motif_enrichment_pipeline.pl |
| Enhancer_both_up_genes_matrix.gz | data/INTEGRATEDPROJ/Enhancer | Compressed matrix for profile plotting of up-regulated genes | computeMatrix (deepTools) |
| Enhancer_both_up_genes_signal.tab | data/INTEGRATEDPROJ/Enhancer | Signal values for up-regulated regions | computeMatrix (deepTools) |
| Enhancer_both_up_genes_plotProfile.pdf | data/INTEGRATEDPROJ/Enhancer | Profile plot of signal around up-regulated peaks | plotProfile (deepTools) |
| Enhancer_both_down_genes_matrix.gz | data/INTEGRATEDPROJ/Enhancer | Compressed matrix for profile plotting of down-regulated genes | computeMatrix (deepTools) |
| Enhancer_both_down_genes_signal.tab | data/INTEGRATEDPROJ/Enhancer | Signal values for down-regulated regions | computeMatrix (deepTools) |
| Enhancer_both_down_genes_plotProfile.pdf | data/INTEGRATEDPROJ/Enhancer | Profile plot of signal around down-regulated peaks | plotProfile (deepTools) |
| enh_gene_within_1mb.bed | data/INTEGRATEDPROJ/Enhancer | BED file of enhancer-gene pairs within 1Mb | calculate_eRNA_gene_distance.pl |
| random_enh_gene_within_1mb.bed | data/INTEGRATEDPROJ/Enhancer | BED file of random enhancer regions | calculate_eRNA_gene_distance.pl |
| random_enh_gene_within_1mb_matrix.gz | data/INTEGRATEDPROJ/Enhancer | Compressed matrix for random enhancer regions profile plotting | computeMatrix (deepTools) |
| random_enh_gene_within_1mb_signal.tab | data/INTEGRATEDPROJ/Enhancer | Signal values for random enhancer regions | computeMatrix (deepTools) |
| enh_gene_within_1mb_signal.tab | data/INTEGRATEDPROJ/Enhancer | Signal values for enhancer-gene pairs within 1Mb | computeMatrix (deepTools) |
| enh_gene_within_1mb_matrix.gz | data/INTEGRATEDPROJ/Enhancer | Compressed matrix for enhancer-gene pairs profile plotting | computeMatrix (deepTools) |
| TF_up_gene_heatmap.pdf | data/INTEGRATEDPROJ/Enhancer | Heatmap of upregulated TFs expression | ComplexHeatmap (R) |
| TF_down_gene_heatmap.pdf | data/INTEGRATEDPROJ/Enhancer | Heatmap of downregulated TFs expression | ComplexHeatmap (R) |
| Motif_result_up_list.txt | data/INTEGRATEDPROJ/Enhancer | List of enriched motifs with TF mappings for up-regulated genes | extr_enriched_motifs.pl |
| Motif_result_down_list.txt | data/INTEGRATEDPROJ/Enhancer | List of enriched motifs with TF mappings for down-regulated genes | extr_enriched_motifs.pl |
| Ture_enh_gene_within_1mb.bed | data/INTEGRATEDPROJ/Enhancer | Filtered enhancer-gene pairs based on signal threshold | Calculate_signal.R |
| Ture_enh_gene_within_1mb_signal.tab | data/INTEGRATEDPROJ/Enhancer | Signal values for filtered enhancer-gene pairs | computeMatrix (deepTools) |
| Ture_enh_gene_within_1mb_profile.pdf | data/INTEGRATEDPROJ/Enhancer | Profile plot of signal around filtered enhancer peaks | plotProfile (deepTools) |
| Ture_enh_gene_within_1mb_matrix.gz | data/INTEGRATEDPROJ/Enhancer | Compressed matrix for filtered enhancer-gene pairs profile plotting | computeMatrix (deepTools) |
| common_up_genes_venn_diagram.pdf | data/INTEGRATEDPROJ/Enhancer | Venn diagram of common up-regulated genes | VennDiagram (R) |
| common_down_genes_venn_diagram.pdf | data/INTEGRATEDPROJ/Enhancer | Venn diagram of common down-regulated genes | VennDiagram (R) |
| common_genes_volcano_plot.pdf | data/INTEGRATEDPROJ/Enhancer | Volcano plot of common DEGs | ggplot2 (R) + ggrepel |
| RNA_ChIP_col_plot_same_direction.pdf | data/INTEGRATEDPROJ/Enhancer | Correlation plot between RNA-seq and ChIP-seq log2FC | ggplot2 (R) |
| {groupname}_up_peaks.bed | data/INTEGRATEDPROJ/Enhancer | BED file of up-regulated peaks | enhancer_cor_plot.R |
| {groupname}_down_peaks.bed | data/INTEGRATEDPROJ/Enhancer | BED file of down-regulated peaks | enhancer_cor_plot.R |
| {groupname}_all_peaks.bed | data/INTEGRATEDPROJ/Enhancer | BED file of all peaks | enhancer_cor_plot.R |
| GO_enrich_common_gene.csv | data/INTEGRATEDPROJ/Enhancer | GO enrichment results for common genes | clusterProfiler (R) |
| KEGG_GSEA_common_gene.csv | data/INTEGRATEDPROJ/Enhancer | KEGG GSEA results for common genes | clusterProfiler (R) |

---

| Tool | Purpose | Documentation |
| :--- | :--- | :--- |
| **VennDiagram (R)** | Visualize overlaps between RNA-seq and ChIP-seq gene sets | [Link](https://cran.r-project.org/web/packages/VennDiagram/index.html) |
| **ggplot2 (R)** | Create publication-quality plots (volcano plots, correlation plots) | [Link](https://ggplot2.tidyverse.org/) |
| **ggrepel** | Prevent label overlapping in ggplot2 visualizations | [Link](https://ggrepel.slowkow.com/) |
| **clusterProfiler (R)** | Perform GO and KEGG enrichment analysis on common genes | [Link](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) |
| **AME (MEME Suite)** | Motif enrichment analysis in Enhancer regions | [Link](https://meme-suite.org/meme/tools/ame) |
| **deeptools (computeMatrix/plotProfile)** | Process BigWig files and generate signal profiles around Enhancer peaks | [Link](https://deeptools.readthedocs.io/) |
| **ComplexHeatmap (R)** | Create heatmaps for TF expression patterns | [Link](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html) |
| **TF_motif_enrichment_pipeline.pl** | Pipeline for TF motif discovery, network building, and TF-gene mapping | [Link](data/INTEGRATEDPROJ/script/Enhancer/TF_motif_enrichment_pipeline.pl) |
| **extr_enriched_motifs.pl** | Extract significantly enriched motifs from AME results | [Link](data/INTEGRATEDPROJ/script/Enhancer/extr_enriched_motifs.pl) |
| **TF_motif_plot.R** | Generate heatmaps for transcription factor expression analysis | [Link](data/INTEGRATEDPROJ/script/Enhancer/TF_motif_plot.R) |
| **calculate_eRNA_gene_distance.pl** | Calculate distances between enhancers and genes, generate random control regions | [Link](data/INTEGRATEDPROJ/script/Enhancer/calculate_eRNA_gene_distance.pl) |
| **Calculate_signal.R** | Filter enhancer-gene pairs based on ChIP-seq signal thresholds | [Link](data/INTEGRATEDPROJ/script/Enhancer/Calculate_signal.R) |
| **enhancer_cor_plot.R** | Integrate RNA-seq and ChIP-seq data, perform comparative analysis and visualization | [Link](data/INTEGRATEDPROJ/script/Enhancer/enhancer_cor_plot.R) |
| **snakefile** | Enhancer analysis workflow automation | [Link](data/INTEGRATEDPROJ/script/Enhancer/snakefile) |

---

## Reference Data
All required reference files are documented in the [reference](/reference) directory.