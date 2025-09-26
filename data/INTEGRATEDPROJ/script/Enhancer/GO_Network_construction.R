#!/usr/bin/env Rscript
rm(list = ls())

# 自动安装并加载包
install_and_load <- function(package, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package %in% c("org.Mm.eg.db", "org.Hs.eg.db", "clusterProfiler")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = repo)
      }
      BiocManager::install(package, update = FALSE)
    } else {
      install.packages(package, repos = repo, dependencies = TRUE)
    }
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

required_packages <- c("dplyr", "tidyr", "igraph", "optparse")
lapply(required_packages, install_and_load)

# 命令行参数
option_list <- list(
  make_option(c("-g", "--go_csv"), type = "character", help = "GO富集分析结果csv文件"),
  make_option(c("-o", "--output"), type = "character", default = ".", help = "输出文件夹，默认为当前目录")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$go_csv)) {
  cat("请使用 -g 参数指定GO富集结果csv文件！\n")
  quit(status = 1)
}
input_file <- opt$go_csv
output_folder <- opt$output

# 读取GO富集结果
GO_kk_df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!("core_enrichment" %in% colnames(GO_kk_df) && "ID" %in% colnames(GO_kk_df) && "pvalue" %in% colnames(GO_kk_df))) {
  cat("输入csv文件缺少 core_enrichment、ID 或 pvalue 列！\n")
  quit(status = 1)
}

## 只选pvalue最小的前10个GO term
GO_kk_df_top10 <- GO_kk_df %>%
  dplyr::arrange(pvalue) %>%
  dplyr::slice(1:10)

# 构建边表
edges <- GO_kk_df_top10 %>%
  dplyr::select(ID, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = core_enrichment)

# 构建节点表
GO_nodes <- data.frame(name = unique(edges$GO_ID), type = "GO")
Gene_nodes <- data.frame(name = unique(edges$Gene), type = "Gene")
nodes <- dplyr::bind_rows(GO_nodes, Gene_nodes)

# 构建igraph对象
g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# 输出graphml文件
output_graphml <- file.path(output_folder, "GO_network.graphml")
igraph::write_graph(g, file = output_graphml, format = "graphml")

cat("graphml文件已输出至：", output_graphml, "\n")

# Rscript your_script.R -g your_GO_result.csv -o ./output_folder/