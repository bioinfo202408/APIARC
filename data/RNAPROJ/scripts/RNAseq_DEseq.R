#!/usr/bin/env Rscript

if(!require("yaml")) install.packages("yaml", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("optparse")) install.packages("optparse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
library(DESeq2)
library(yaml)
library(optparse)

.libPaths(c("/usr/local/R/lib64/R/library", .libPaths()))

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-i", "--input"), type = "character"),
  make_option(c("-y", "--yaml"), type = "character"),
  make_option(c("-o", "--output"), type = "character")
)))

database_all <- read.table(file = opt$input, sep = ",", header = TRUE)
config <- yaml.load_file(opt$yaml)

expdata <- database_all[, -1]
rownames(expdata) <- database_all[, 1]
sample_names <- names(config$group$sample_IDs)
conditions <- factor(unlist(config$group$sample_IDs), levels = c("T", "P"))

coldata <- data.frame(conditions = conditions, row.names = sample_names)
expdata <- expdata[, sample_names]
expdata <- expdata[complete.cases(expdata), ]

dds <- DESeqDataSetFromMatrix(countData = expdata, colData = coldata, design = ~ conditions)
dds <- DESeq(dds)

res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
row.names(resdata) <- resdata[, 1]
names(resdata)[names(resdata) == "Row.names"] <- "SYMBOL"
resdata$gene_id <- sub("\\|[^|]*$", "", resdata$SYMBOL)
resdata$SYMBOL <- sub("^.*\\|", "", resdata$SYMBOL)
resdata <- resdata[!is.na(resdata$padj), ]

write.csv(resdata, opt$output, row.names = FALSE)


# Rscript RNAseq_DEseq.R --input /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632_pipline/mRNA/gene_count_matrix.csv --yaml /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632_pipline/config.yaml --output /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632_pipline/mRNA/DEG_result.csv.csv

