#!/usr/bin/env Rscript

install_and_load <- function(package, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") {
    if (!requireNamespace(package, quietly = TRUE)) {
        if (package %in% installed.packages()) {
            suppressWarnings(library(package, character.only = TRUE))
        } else {
            if (package %in% c("yaml", "GenomicRanges", "TxDb.Mmusculus.UCSC.mm10.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Mm.eg.db", "org.Hs.eg.db", "rtracklayer", "ChIPseeker")) {
                BiocManager::install(package)
            } else {
                install.packages(package, repos = repo)
            }
        }
    } else {
        library(package, character.only = TRUE)
    }
}

install_and_load("BiocManager")
required_packages <- c(
    "optparse", "GenomicRanges", "TxDb.Mmusculus.UCSC.mm10.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Mm.eg.db", "org.Hs.eg.db", "rtracklayer", "ChIPseeker", "RColorBrewer", "ggplot2", "stringr", "tidyr", "dplyr", "yaml"
)
lapply(required_packages, install_and_load)

library(optparse)
library(yaml)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)

option_list <- list(
    make_option(c("-d", "--directory"), type = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL),
    make_option(c("-p", "--pdf"), type = "character", default = NULL),
    make_option(c("-y", "--yaml"), type = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$yaml)) quit(status = 1)

yaml_data <- yaml.load_file(opt$yaml)
first_sample <- names(yaml_data)[1]
genome <- yaml_data[[first_sample]]$genome

if (genome == "mm") {
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    annoDb <- "org.Mm.eg.db"
} else if (genome == "hg") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    annoDb <- "org.Hs.eg.db"
} else {
    quit(status = 1)
}

if (is.null(opt$output)) opt$output <- file.path(opt$directory, "annotations")
if (is.null(opt$pdf)) opt$pdf <- file.path(opt$directory, "pdf_plots")
dir.create(opt$output, showWarnings = FALSE)
dir.create(opt$pdf, showWarnings = FALSE)

merge_bed_files <- function(bed_list, min_overlap = 0.5) {
    common_seqlevels <- Reduce(intersect, lapply(bed_list, seqlevels))
    bed_list <- lapply(bed_list, function(bed) keepSeqlevels(bed, common_seqlevels, pruning.mode = "coarse"))
    merged_bed <- bed_list[[1]]
    for (i in seq_along(bed_list[-1])) {
        current_bed <- bed_list[[i + 1]]
        current_min_overlap <- as.integer(min(width(current_bed)) * min_overlap)
        merged_bed <- subsetByOverlaps(merged_bed, current_bed, type = "any", minoverlap = current_min_overlap)
    }
    return(merged_bed)
}

extract_bowtie2_stats <- function(file) {
    lines <- readLines(file, warn = FALSE)
    end_lines <- grep("======== SRR[0-9]+ Bowtie2 end ========", lines)
    result <- list()
    for (i in end_lines) {
        sample <- str_extract(lines[i], "SRR[0-9]+")
        stat_block <- lines[(i - 10):(i - 1)]
        unique <- as.numeric(str_extract(grep("aligned concordantly exactly 1 time", stat_block, value = TRUE), "\\d+\\.\\d+(?=%)"))
        multiple <- as.numeric(str_extract(grep("aligned concordantly >1 times", stat_block, value = TRUE), "\\d+\\.\\d+(?=%)"))
        unmapped <- as.numeric(str_extract(grep("aligned concordantly 0 times", stat_block, value = TRUE), "\\d+\\.\\d+(?=%)"))
        result[[length(result) + 1]] <- data.frame(sample = sample, Unique = unique, Multiple = multiple, Unmapped = unmapped)
    }
    if (length(result) == 0) return(NULL)
    return(do.call(rbind, result))
}

log_files <- list.files(opt$directory, pattern = "bowtie.log", recursive = TRUE, full.names = TRUE)
stats_list <- lapply(log_files, extract_bowtie2_stats)
stats_list <- stats_list[!sapply(stats_list, is.null)]
stats_df <- bind_rows(stats_list)

result <- pivot_longer(stats_df, cols = c("Unique", "Multiple", "Unmapped"), names_to = "type", values_to = "percent")
result$type <- factor(result$type, levels = c("Unique", "Multiple", "Unmapped"))
sample_count <- length(unique(result$sample))
final_width <- max(8, min(12 + 0.5 * (sample_count - 5), 30))
pdf(file.path(opt$pdf, "Bowtie2_result_bar_plot.pdf"), width = final_width, height = 5)
p <- ggplot(result, aes(x = sample, y = percent, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    geom_text(aes(label = sprintf("%.2f%%", percent)), position = position_dodge(width = 0.6), vjust = -0.3, size = 5) +
    scale_fill_manual(values = c("Unique" = "#4daf4a", "Multiple" = "#377eb8", "Unmapped" = "#e41a1c")) +
    ylim(0, 100) +
    labs(title = "Bowtie2 Mapping Statistics", y = "Percentage (%)", fill = "Type") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
print(p)
dev.off()

macs2_dir <- file.path(opt$directory, "Macs2")
sample_dirs <- list.dirs(macs2_dir, recursive = FALSE)
sample_groups <- list()
for (dir in sample_dirs) {
    sample_name <- basename(dir)
    group_name <- sub("_rep\\d+$", "", sample_name)
    sample_groups[[group_name]] <- c(sample_groups[[group_name]], dir)
}
for (group_name in names(sample_groups)) {
    group_output_dir <- file.path(opt$output, group_name)
    dir.create(group_output_dir, showWarnings = FALSE)
    bed_files <- unlist(lapply(sample_groups[[group_name]], function(dir) list.files(dir, pattern = "_peaks_tab\\.bed$", full.names = TRUE)))
    bed_list <- lapply(bed_files, function(file) import(file, format = "BED"))
    merged_bed <- merge_bed_files(bed_list, min_overlap = 0.5)
    merged_bed_file <- file.path(group_output_dir, paste0(group_name, "_merged_peaks_tab.bed"))
    export(merged_bed, merged_bed_file, format = "BED")
    peakAnno <- annotatePeak(readPeakFile(merged_bed_file), tssRegion = c(-2000, 2000), TxDb = txdb, annoDb = annoDb)
    annotationData <- as.data.frame(peakAnno)
    annotationData$annotation <- ifelse(grepl("Promoter", annotationData$annotation), "Promoter",
                                        ifelse(grepl("Intron", annotationData$annotation), "Intron",
                                               ifelse(grepl("Exon", annotationData$annotation), "Exon",
                                                      ifelse(grepl("3' UTR", annotationData$annotation), "3' UTR", "Intergenic"))))
    annotationCounts <- as.data.frame(table(annotationData$annotation))
    annotationCounts$Percentage <- annotationCounts$Freq / sum(annotationCounts$Freq) * 100
    total_peaks <- sum(annotationCounts$Freq)
    pdf(file.path(opt$pdf, paste0(group_name, "_annotation_plot.pdf")), width = 7.5, height = 7.5)
    pie(annotationCounts$Freq, labels = paste(annotationCounts$Var1, "\n(", annotationCounts$Freq, ", ", round(annotationCounts$Percentage, 2), "%)", sep = ""), col = brewer.pal(length(annotationCounts$Var1), "Set3"), border = "white", main = paste("Total peaks:", total_peaks))
    dev.off()
    peak_anno_file <- file.path(group_output_dir, "peak_anno.csv")
    write.csv(annotationData, peak_anno_file, row.names = FALSE)
    promoter_peaks <- annotationData[grepl("Promoter", annotationData$annotation), ]
    promoter_file <- file.path(group_output_dir, "promoter_peaks.csv")
    write.csv(promoter_peaks, promoter_file, row.names = FALSE)
}


# Rscript annoChIPPeaks.R --directory /home/yxiaobo/RCProj/ChIPseq/database/TF-binding/GSE85632_text --yaml /home/yxiaobo/RCProj/ChIPseq/database/TF-binding/GSE85632_text/bambw_list.yaml --output /home/yxiaobo/RCProj/ChIPseq/database/TF-binding/GSE85632_text --pdf /home/yxiaobo/RCProj/picture/ChIPseq/TF-binding/GSE85632_text
