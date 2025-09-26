# reference
The `reference` directory provides reference and annotation files for the entire APIARC analysis pipeline. The file names, species, location, usage, and a brief description are listed in the table below.

| File Name | Species | Location | Usage | Description |
| :--- | :--- | :--- | :--- | :--- |
| gencode.v44.annotation.gtf | Human | reference/annotation/ | [data/RNAPROJ/config.yaml](data/RNAPROJ/config.yaml) | Annotation file containing detailed information of genes, exons, etc. |
| gencode.vM25.annotation.gtf | Mouse | reference/annotation/ | [data/RNAPROJ/config.yaml](data/RNAPROJ/config.yaml) | Annotation file containing detailed information of genes, exons, etc. |
| gencode.v44.annotation.bed | Human | reference/annotation/ | [data/RNAPROJ/config.yaml](data/RNAPROJ/config.yaml) | File containing only genomic coordinates of gene regions |
| gencode.vM25.annotation.bed | Mouse | reference/annotation/ | [data/RNAPROJ/config.yaml](data/RNAPROJ/config.yaml) | File containing only genomic coordinates of gene regions |
| GRCh38 | Human | reference/hisat2_index/humanhisat2Index/ | [data/RNAPROJ/config.yaml](data/RNAPROJ/config.yaml) | Genome sequence index file for HISAT2 alignment software |
| GRCm38 | Mouse | reference/hisat2_index/mousehisat2Index/ | [data/RNAPROJ/config.yaml](data/RNAPROJ/config.yaml) | Genome sequence index file for HISAT2 alignment software |
| GRCh38 | Human | reference/bowtie2_index/humanbowtie2Index/ | [data/CHIPPROJ/config.yaml](data/CHIPPROJ/config.yaml) | Genome sequence index file for Bowtie2 alignment software |
| GRCm38 | Mouse | reference/bowtie2_index/mousebowtie2Index/ | [data/CHIPPROJ/config.yaml](data/CHIPPROJ/config.yaml) | Genome sequence index file for Bowtie2 alignment software |
| hsa_KEGG_anno.rds | Human | reference/genome/ | enhancer_cor_plot.R, Promoter_plot_1.R (* Built-in, no modification needed) | Human KEGG pathway annotation file (R format) |
| mmu_KEGG_anno.rds | Mouse | reference/genome/ | enhancer_cor_plot.R, Promoter_plot_1.R (* Built-in, no modification needed) | Mouse KEGG pathway annotation file (R format) |
| GRCh38.p12.genome.fa.gz | Human | reference/genome/ | TF_motif_enrichment_pipeline.pl (* Built-in, no modification needed) | Human reference genome sequence file (compressed) |
| GRCm38.p6.genome.fa.gz | Mouse | reference/genome/ | TF_motif_enrichment_pipeline.pl (* Built-in, no modification needed) | Mouse reference genome sequence file (compressed) |
| Homo_sapiens_TF.txt | Human | reference/Motif_tf_anno/ | TF_motif_enrichment_pipeline.pl, extr_enriched_motifs.pl (* Built-in, no modification needed) | Human transcription factor gene list file |
| Mus_musculus_TF.txt | Mouse | reference/Motif_tf_anno/ | TF_motif_enrichment_pipeline.pl, extr_enriched_motifs.pl (* Built-in, no modification needed) | Mouse transcription factor gene list file |
| JASPAR_CORE | Multiple | reference/Motif_tf_anno/ | snakefile both in Promoter and Enhancer (* Required, no modification needed) | JASPAR transcription factor motif database file |
| JASPAR_CORE_ID_to_NAME.txt | Multiple | reference/Motif_tf_anno/ | snakefile both in Promoter and Enhancer (* Required, no modification needed) | JASPAR database ID to name mapping file |
| refGene.ncbiRefSeq.genecode.filtered.bed | Multiple Species | reference/Enhancer_anno/ | snakefile in Enhancer (* Required, no modification needed) | Filtered and merged gene annotation BED file |
| picard-2.18.2 | - | reference/ | config.yaml both in RNAPROJ and CHIPPROJ (* Required, no modification needed) | Picard tools software package |