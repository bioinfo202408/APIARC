#!/usr/bin/env Rscript
rm(list = ls())

install_and_load <- function(package, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package %in% c("org.Mm.eg.db", "org.Hs.eg.db", "clusterProfiler", "enrichplot")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = repo)
      }
      BiocManager::install(package, update = FALSE)
    } else {
      install.packages(package, repos = repo, dependencies = TRUE)
    }
  }
  library(package, character.only = TRUE)
}

required_packages <- c(
  "dplyr", "org.Mm.eg.db", "org.Hs.eg.db", "VennDiagram", "grid", "ggplot2",
  "clusterProfiler", "enrichplot", "ggrepel", "yaml", "optparse", "futile.logger"
)
lapply(required_packages, install_and_load)

option_list <- list(
  make_option(c("-r", "--RNAseq"), type = "character"),
  make_option(c("-c", "--ChIPseq"), type = "character"),
  make_option(c("-o", "--output"), type = "character"),
  make_option(c("-p", "--picture"), type = "character"),
  make_option(c("-y", "--yaml"), type = "character"),
  make_option(c("-m", "--matrix"), type = "character"),
  make_option(c("-e", "--enh"), type = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

RNAseq_gene_file <- opt$RNAseq
ChIPseq_gene_file <- opt$ChIPseq
output_folder <- opt$output
outpic <- opt$picture
yaml_file <- opt$yaml
matrix_file <- opt$matrix
enhancer_bed_file <- opt$enh

if (!file.exists(yaml_file)) quit(status=1)
if (!file.exists(matrix_file)) quit(status=1)
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
if (!dir.exists(outpic)) dir.create(outpic, recursive = TRUE)
rna_file <- file.path(RNAseq_gene_file)
if (!file.exists(rna_file)) quit(status=1)
RNAseq_gene <- read.csv(rna_file, stringsAsFactors = FALSE)
RNAseq_gene <- RNAseq_gene[!is.na(RNAseq_gene$padj), ]
yaml_info <- yaml.load_file(yaml_file)
if (is.null(yaml_info$groups)) quit(status=1)
chip_csv_files <- file.path(ChIPseq_gene_file, yaml_info$groups, "peak_anno.csv")
names(chip_csv_files) <- yaml_info$groups

chip_list <- lapply(names(chip_csv_files), function(sample_name) {
  f <- chip_csv_files[sample_name]
  if (file.exists(f)) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    return(df)
  } else {
    return(NULL)
  }
})
names(chip_list) <- yaml_info$groups
chip_list <- Filter(Negate(is.null), chip_list)
if (length(chip_list) == 0) quit(status=1)
if (!"file" %in% names(yaml_info)) quit(status=1)
first_macs_sample <- names(yaml_info$file)[1]
genome_type <- yaml_info$file[[first_macs_sample]]$genome
if (length(genome_type) != 1) quit(status=1)
anno_file <- switch(
  genome_type,
  mm = "reference/annotation/gencode.vM25.annotation.bed",
  hg = "reference/annotation/gencode.v44.annotation.bed",
  quit(status=1)
)
if (!file.exists(anno_file)) quit(status=1)
anno <- read.table(anno_file, header = FALSE, stringsAsFactors = FALSE)
matrix_data <- read.table(matrix_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(matrix_data) < 2) quit(status=1)
file_paths <- matrix_data[, 2]
sample_names <- matrix_data[, 1]
mean_matrix_list <- list()
for (i in seq_along(file_paths)) {
  f <- file_paths[i]
  if (!file.exists(f)) next
  mat <- read.table(
    f, header = TRUE, stringsAsFactors = FALSE, 
    skip = 2, fill = TRUE, check.names = FALSE
  )
  mat[is.na(mat)] <- 0
  unique_names <- unique(colnames(mat))
  mean_mat <- sapply(unique_names, function(name) {
    cols <- colnames(mat) == name
    rowMeans(mat[, cols, drop = FALSE])
  })
  mean_matrix_list[[sample_names[i]]] <- as.data.frame(mean_mat)
}
if (length(mean_matrix_list) == 0) quit(status=1)
ChIP_signal <- lapply(mean_matrix_list, function(x) cbind(anno, x))
enhancer_bed <- read.table(enhancer_bed_file, header = FALSE, stringsAsFactors = FALSE)
unique_symbols <- unique(enhancer_bed[, 5])
RNAseq_gene <- RNAseq_gene %>% filter(SYMBOL %in% unique_symbols)
RNAseq_gene_ids <- sub("\\..*", "", RNAseq_gene$gene_id)
chip_gene_ids_list <- lapply(chip_list, function(df) sub("\\..*", "", df$ENSEMBL))
chip_all_genes <- unique(unlist(chip_gene_ids_list))
gene_intersect <- intersect(RNAseq_gene_ids, chip_all_genes)
if (length(gene_intersect) == 0) quit(status=1)
ChIP_signal <- lapply(ChIP_signal, function(df) {
  gene_id_simple <- sub("\\..*", "", df$V4)
  df[gene_id_simple %in% gene_intersect, ]
})
group_map <- yaml_info$groupnames
if (is.null(group_map)) quit(status=1)

get_group_means <- function(df, group_name, input_pattern = "Input") {
  ip_cols <- grep(group_name, colnames(df), value = TRUE)
  ip_cols <- setdiff(ip_cols, grep(input_pattern, colnames(df), value = TRUE))
  input_cols <- grep(input_pattern, colnames(df), value = TRUE)
  if (length(ip_cols) == 0) quit(status=1)
  if (length(input_cols) == 0) quit(status=1)
  IP_mean <- rowMeans(df[, ip_cols, drop = FALSE])
  Input_mean <- rowMeans(df[, input_cols, drop = FALSE])
  tmp <- data.frame(
    gene = df$V4,
    IP_mean, Input_mean,
    check.names = FALSE
  )
  colnames(tmp)[2:3] <- paste0(group_name, c("_IP_mean", "_Input_mean"))
  return(tmp)
}
valid_groups <- intersect(names(group_map), names(ChIP_signal))
if (length(valid_groups) == 0) quit(status=1)
group_results <- lapply(valid_groups, function(grp) {
  get_group_means(ChIP_signal[[grp]], grp)
})
names(group_results) <- valid_groups
if (length(group_results) < 2) quit(status=1)
Result <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), group_results)
P_group_name <- names(group_map)[group_map == "P"]
T_group_name <- names(group_map)[group_map == "T"]
if (length(P_group_name) != 1) quit(status=1)
if (length(T_group_name) != 1) quit(status=1)
P_ip_col <- paste0(P_group_name, "_IP_mean")
T_ip_col <- paste0(T_group_name, "_IP_mean")
if (!P_ip_col %in% colnames(Result)) quit(status=1)
if (!T_ip_col %in% colnames(Result)) quit(status=1)
Result <- Result[Result[[T_ip_col]] > 0, ]
if (nrow(Result) == 0) quit(status=1)
Result$log2FC <- log2(Result[[P_ip_col]] / Result[[T_ip_col]])
Result$gene_id_simple <- sub("\\..*", "", Result$gene)
RNAseq_gene$gene_id_simple <- sub("\\..*", "", RNAseq_gene$gene_id)
Result <- Result %>%
  left_join(RNAseq_gene, by = "gene_id_simple") %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))
if (nrow(Result) == 0) quit(status=1)
RNAseq_up   <- Result %>% filter(padj < 0.05 & log2FoldChange > 1) %>% pull(SYMBOL)  
RNAseq_down <- Result %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(SYMBOL)
ChIPseq_up   <- Result %>% filter(log2FC > 0.5) %>% pull(SYMBOL)
ChIPseq_down <- Result %>% filter(log2FC < -0.5) %>% pull(SYMBOL)
futile.logger::flog.threshold(futile.logger::ERROR)
venn_up <- venn.diagram(
  x = list(
    RNAseq  = RNAseq_up,
    ChIPseq = ChIPseq_up
  ),
  fill = c("#fc8d62", "#66c2a5"),
  alpha = c(0.5, 0.5),
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  scaled = FALSE,
  main = paste("Common Up-regulated Genes in Promoter Region: ", length(intersect(RNAseq_up, ChIPseq_up))),
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans",
  margin = 0.1,
  filename = NULL
)
venn_pdf_up <- file.path(outpic, "common_up_genes_venn_diagram.pdf")
pdf(file = venn_pdf_up, width = 9, height = 6)
grid.draw(venn_up)
dev.off()
venn_down <- venn.diagram(
  x = list(
    RNAseq = RNAseq_down,
    ChIPseq = ChIPseq_down
  ),
  fill = c("#66c2a5", "#fc8d62"),
  alpha = c(0.5, 0.5),
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  scaled = FALSE,
  main = paste("Common Down-regulated Genes: ", length(intersect(RNAseq_down, ChIPseq_down))),
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans",
  margin = 0.1,
  filename = NULL
)
venn_pdf_down <- file.path(outpic, "common_down_genes_venn_diagram.pdf")
pdf(file = venn_pdf_down, width = 9, height = 6)
grid.draw(venn_down)
dev.off()
plot_df <- Result %>%
  filter(!is.na(padj) & padj > 0 & -log10(padj) <= 200) %>%
  mutate(group = case_when(
    padj < 0.05 & log2FoldChange > 1  ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "NS"
  )) %>%
  mutate(group = factor(group, levels = c("Up", "NS", "Down")))
topn <- 10
logfc <- 1
padj_cutoff <- 0.05
plot_df_nonzero <- plot_df %>% filter(padj != 0)
top_up <- plot_df_nonzero %>% filter(group == "Up") %>% arrange(padj) %>% slice_head(n = topn)
top_down <- plot_df_nonzero %>% filter(group == "Down") %>% arrange(padj) %>% slice_head(n = topn)
top_genes <- bind_rows(top_up, top_down)
volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = log2FoldChange, size = -log10(padj)), alpha = 0.7) +
  geom_point(
    data = top_genes,
    aes(color = log2FoldChange, size = -log10(padj)),
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = SYMBOL),
    size = 3, box.padding = 0.5, point.padding = 0.8,
    segment.color = "black", show.legend = FALSE, max.overlaps = Inf
  ) +
  scale_color_gradientn(
    colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
    values = scales::rescale(c(min(plot_df$log2FoldChange, na.rm = T), -1, 0, 1, max(plot_df$log2FoldChange, na.rm = T)))
  ) +
  scale_size(range = c(1, 7)) +
  geom_vline(xintercept = c(-logfc, logfc), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  xlim(c(-max(abs(plot_df$log2FoldChange), na.rm = T), max(plot_df$log2FoldChange, na.rm = T))) +
  ylim(c(0, max(-log10(plot_df$padj), na.rm = T) + 2)) +
  labs(
    x = "log2 (fold change)",
    y = "-log10 (padj)",
    title = "Volcano Plot (Common Genes)",
    subtitle = paste("Top", topn, "up/down regulated genes")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
volcano_pdf <- file.path(outpic, "common_genes_volcano_plot.pdf")
pdf(file = volcano_pdf, width = 10, height = 7.5)
print(volcano_plot)
dev.off()
result_df_same_direction <- Result[
  (Result$log2FoldChange > 0 & Result$log2FC > 0) |
    (Result$log2FoldChange < 0 & Result$log2FC < 0),
  ,]
result_df_same_direction <- result_df_same_direction[ 
  is.finite(result_df_same_direction$log2FoldChange) & 
    is.finite(result_df_same_direction$log2FC), 
  ,]
result_df_same_direction$color_group <- "Non-significant"
result_df_same_direction$color_group[
  result_df_same_direction$log2FoldChange > 1 & result_df_same_direction$log2FC > 0.5
] <- "Up"
result_df_same_direction$color_group[
  result_df_same_direction$log2FoldChange < -1 & result_df_same_direction$log2FC < -0.5
] <- "Down"
color_map <- c("Up" = "#D55E00", "Down" = "#0072B2", "Non-significant" = "grey70")
fit <- lm(log2FC ~ 0 + log2FoldChange, data=result_df_same_direction)
k <- coef(fit)[1]
cor_value <- cor(result_df_same_direction$log2FoldChange, result_df_same_direction$log2FC, method = "pearson")
cor_label <- sprintf("Pearson~italic(R) == %.2f", cor_value)
cor_pdf <- file.path(outpic, "RNA_ChIP_col_plot_same_direction.pdf")
pdf(file = cor_pdf, width=5, height=4)
ggplot(result_df_same_direction, aes(x=log2FoldChange, y=log2FC, color=color_group)) +
  geom_point(alpha=0.7, size=2) +
  scale_color_manual(values=color_map) +
  geom_abline(intercept=0, slope=k, color="black", linewidth=1) +
  geom_hline(yintercept=0, linetype="dashed", color="grey") +
  geom_vline(xintercept=0, linetype="dashed", color="grey") +
  theme_classic() +
  labs(x="RNA-seq log2FC", y="ChIP-seq log2FC", title="Same Direction (Up/Down)") +
  annotate("text", 
           x = min(result_df_same_direction$log2FoldChange, na.rm=TRUE), 
           y = max(result_df_same_direction$log2FC, na.rm=TRUE), 
           label = cor_label, 
           parse=TRUE, 
           hjust=0, vjust=1, size=5, color="black")
dev.off()
id2symbol <- dplyr::select(Result, gene_id_simple, SYMBOL) %>% distinct()
RNAseq_up_id   <- Result %>% filter(padj < 0.05 & log2FoldChange > 1) %>% pull(SYMBOL)
RNAseq_down_id <- Result %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(SYMBOL)
groupnames <- yaml_info$groupnames
p_sample <- names(groupnames)[groupnames == "P"]
t_sample <- names(groupnames)[groupnames == "T"]
p_peak <- chip_list[[p_sample]]
t_peak <- chip_list[[t_sample]]

process_peak_data <- function(peak_data, up_ids, down_ids, sample_name) {
  bed_up <- peak_data %>%
    filter(SYMBOL %in% up_ids) %>%
    distinct(seqnames, start, end, SYMBOL, .keep_all = TRUE) %>%
    dplyr::select(seqnames, start, end, SYMBOL)
  colnames(bed_up) <- c("chrom", "chromStart", "chromEnd", "name")
  write.table(bed_up,
              file = file.path(output_folder, paste0(sample_name, "_up_peaks.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  bed_down <- peak_data %>%
    filter(SYMBOL %in% down_ids) %>%
    distinct(seqnames, start, end, SYMBOL, .keep_all = TRUE) %>%
    dplyr::select(seqnames, start, end, SYMBOL)
  colnames(bed_down) <- c("chrom", "chromStart", "chromEnd", "name")
  write.table(bed_down,
              file = file.path(output_folder, paste0(sample_name, "_down_peaks.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

process_peak_data_all <- function(peak_data, sample_name) {
  bed_all <- peak_data %>%
    distinct(seqnames, start, end, SYMBOL, .keep_all = TRUE) %>%
    dplyr::select(seqnames, start, end, SYMBOL)
  colnames(bed_all) <- c("chrom", "chromStart", "chromEnd", "name")
  write.table(bed_all,
              file = file.path(output_folder, paste0(sample_name, "_all_peaks.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

process_peak_data(p_peak, RNAseq_up_id, RNAseq_down_id, p_sample)
process_peak_data(t_peak, RNAseq_up_id, RNAseq_down_id, t_sample)
process_peak_data_all(p_peak, p_sample)
process_peak_data_all(t_peak, t_sample)

geneList <- Result$log2FoldChange
names(geneList) <- Result$SYMBOL
geneList <- sort(geneList, decreasing = TRUE)

if (genome_type == "hg") {
  OrgDb <- org.Hs.eg.db
} else if (genome_type == "mm") {
  OrgDb <- org.Mm.eg.db
} else {
  quit(status=1)
}

gene_df <- tryCatch(
  bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb),
  error = function(e) quit(status=1)
)

geneList2 <- geneList[gene_df$SYMBOL]
names(geneList2) <- gene_df$ENTREZID

GO_kk <- gseGO(
  geneList      = geneList2,
  OrgDb         = OrgDb,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1
)

GO_kk_readable <- setReadable(GO_kk, OrgDb = OrgDb, keyType = 'ENTREZID')
GO_kk_df <- as.data.frame(GO_kk_readable)
write.csv(GO_kk_df, file = file.path(output_folder, "GO_enrich_common_gene.csv"), row.names = FALSE)

kegg_file <- switch(
  genome_type,
  mm = "reference/genome/mmu_KEGG_anno.rds",
  hg = "reference/genome/hsa_KEGG_anno.rds",
  quit(status=1)
)
if (!file.exists(kegg_file)) quit(status=1)
kegg_anno <- readRDS(kegg_file)

KEGG_kk_entrez <- GSEA(
  geneList      = geneList2,
  TERM2GENE     = kegg_anno$TERM2GENE,
  TERM2NAME     = kegg_anno$TERM2NAME,
  pvalueCutoff  = 1
)
KEGG_kk <- DOSE::setReadable(
  KEGG_kk_entrez,
  OrgDb = OrgDb,
  keyType = 'ENTREZID'
)
KEGG_kk_df <- as.data.frame(KEGG_kk)
write.csv(KEGG_kk_df, file = file.path(output_folder, "KEGG_GSEA_common_gene.csv"), row.names = FALSE)

KEGG_kk_df_top10 <- KEGG_kk_df %>%
  dplyr::arrange(pvalue) %>%
  dplyr::slice(1:min(10, nrow(.)))

outfile_pdf <- file.path(outpic, "KEGG_GSEA_Top10_Multipage.pdf")
pdf(outfile_pdf, width = 10, height = 5)
for (i in seq_len(nrow(KEGG_kk_df_top10))) {
  geneSetID <- KEGG_kk_df_top10$ID[i]
  plot_title <- KEGG_kk_df_top10$Description[i]
  p1 <- gseaplot2(
    KEGG_kk_entrez,
    geneSetID = geneSetID,
    color = 'red',
    rel_heights = c(1.5, 0.5, 1),
    subplots = 1:3,
    pvalue_table = TRUE,
    title = plot_title,
    ES_geom = 'line'
  )
  print(p1)
}
dev.off()