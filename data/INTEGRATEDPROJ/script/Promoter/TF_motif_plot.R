#!/usr/bin/env Rscript

rm(list = ls())

install_and_load <- function(package, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package %in% installed.packages()) {
      suppressWarnings(library(package, character.only = TRUE))
    } else {
      if (package %in% c("org.Mm.eg.db")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = repo)
        }
        BiocManager::install(package, ask = FALSE)
      } else {
        install.packages(package, repos = repo)
      }
    }
  }
  suppressWarnings(suppressPackageStartupMessages(library(package, character.only = TRUE)))
}

required_packages <- c(
  "dplyr",
  "ComplexHeatmap",
  "circlize",
  "optparse"
)
lapply(required_packages, install_and_load)

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(optparse)

option_list <- list(
  make_option(c("-r", "--RNAseq"), type = "character"),
  make_option(c("-c", "--TF_gene"), type = "character"),
  make_option(c("-o", "--output"), type = "character"),
  make_option(c("-p", "--picture"), type = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$RNAseq) | is.null(opt$TF_gene) | is.null(opt$output) | is.null(opt$picture)) quit(status=1)

RNAseq_gene_file <- opt$RNAseq
TF_gene_file <- opt$TF_gene
output_folder <- opt$output
heatmap_pdf <- opt$picture

TF_motif_data <- read.table(TF_gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(TF_motif_data) <- c("Ensembl", "TF", "pvalue")
TF_unique <- TF_motif_data[!duplicated(TF_motif_data$TF), ]
TF_filtered <- TF_unique[TF_unique$pvalue < 0.05, ]

RNAseq_gene <- read.csv(RNAseq_gene_file, stringsAsFactors = FALSE)
RNAseq_gene <- RNAseq_gene[!is.na(RNAseq_gene$padj), ]
RNAseq_gene$gene_id_noversion <- sub("\\..*", "", RNAseq_gene$gene_id)

non_sample_columns <- c("SYMBOL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_id", "gene_id_noversion")
samples <- setdiff(colnames(RNAseq_gene), non_sample_columns)

TF_gene_data_merged <- merge(
  TF_filtered,
  RNAseq_gene[, c("gene_id_noversion", "log2FoldChange", "padj", samples)],
  by.x = "Ensembl",
  by.y = "gene_id_noversion",
  all.x = TRUE
)

if (nrow(TF_gene_data_merged) == 0) quit(status = 0)
TF_gene_data_merged <- na.omit(TF_gene_data_merged)
if (nrow(TF_gene_data_merged) == 0) quit(status = 0)

expData <- TF_gene_data_merged[, samples, drop = FALSE]
rownames(expData) <- TF_gene_data_merged$TF
expData <- na.omit(expData)
if (nrow(expData) < 2 || ncol(expData) == 0) quit(status = 0)

valscaling <- function(x) { x / sum(x) }
zscore <- function(x) { (x - mean(x)) / sd(x) }

expData_scaled <- apply(expData, 2, valscaling)
if (is.null(dim(expData_scaled))) expData_scaled <- t(as.matrix(expData_scaled))
expData_scaled <- t(apply(expData_scaled, 1, zscore))
if (is.null(dim(expData_scaled))) expData_scaled <- t(as.matrix(expData_scaled))
expData_scaled <- as.matrix(expData_scaled)
rownames(expData_scaled) <- rownames(expData)

if (nrow(expData_scaled) < 2) quit(status = 0)

TF_gene_data_merged <- TF_gene_data_merged[TF_gene_data_merged$TF %in% rownames(expData_scaled), ]
UpTFs <- rownames(expData_scaled)[which(TF_gene_data_merged$log2FoldChange > 1)]
DownTFs <- rownames(expData_scaled)[which(TF_gene_data_merged$log2FoldChange < -1)]
matchIndexes <- match(c(UpTFs, DownTFs), rownames(expData_scaled))

f1 <- colorRamp2(seq(min(expData_scaled), max(expData_scaled), length = 4),
                 c("#CCFFFF", "#99CCFF", "#FF9966", "#CC3333"))

pdf(heatmap_pdf, width = 4, height = 10)
ha <- rowAnnotation(foo = anno_mark(at = matchIndexes, labels = rownames(expData_scaled)[matchIndexes]))
if (nrow(expData_scaled) >= 2) {
  km_num <- min(2, nrow(expData_scaled))
  if (km_num >= 2 && nrow(expData_scaled) > 2) {
    Heatmap(expData_scaled, name = "z-score",
            cluster_columns = FALSE,
            show_row_names = FALSE,
            right_annotation = ha,
            col = f1,
            row_names_gp = gpar(fontsize = 8),
            row_km = km_num)
  } else {
    Heatmap(expData_scaled, name = "z-score",
            cluster_columns = FALSE,
            show_row_names = FALSE,
            right_annotation = ha,
            col = f1,
            row_names_gp = gpar(fontsize = 8))
  }
}
dev.off()

