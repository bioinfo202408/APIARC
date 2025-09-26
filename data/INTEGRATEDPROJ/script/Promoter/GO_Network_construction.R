#!/usr/bin/env Rscript
rm(list = ls())
install_and_load <- function(package, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = repo, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}
pkgs <- c("optparse", "igraph")
lapply(pkgs, install_and_load)

library(optparse)
option_list <- list(
  make_option(c("-g", "--go_csv"), type = "character"),
  make_option(c("-o", "--output"), type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$go_csv)) quit(status = 1)

go_df <- read.csv(opt$go_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("ID", "Description", "pvalue", "p.adjust", "qvalue", "geneID")
if (!all(required_cols %in% colnames(go_df))) quit(status=1)

go_sub <- go_df[order(go_df$pvalue), ][1:min(10, nrow(go_df)), ]

edges <- data.frame()
for (i in seq_len(nrow(go_sub))) {
  term <- go_sub$ID[i]
  genes <- strsplit(go_sub$geneID[i], "/")[[1]]
  edges <- rbind(edges, data.frame(from=term, to=genes, stringsAsFactors=FALSE))
}
if(nrow(edges)==0) quit(status=0)

terms <- unique(edges$from)
genes <- unique(edges$to)
nodes <- data.frame(
  name = c(terms, genes),
  type = c(rep("GO_term", length(terms)), rep("Gene", length(genes))),
  stringsAsFactors = FALSE
)
go_ann <- go_sub[, c("ID", "Description", "pvalue", "p.adjust", "qvalue")]
go_ann <- go_ann[match(terms, go_ann$ID), ]
nodes <- merge(nodes, go_ann, by.x="name", by.y="ID", all.x=TRUE)
nodes$negLog10Pvalue <- ifelse(!is.na(nodes$pvalue), -log10(nodes$pvalue), NA)

g <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=nodes)

output_folder <- opt$output
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
output_graphml <- file.path(output_folder, "GO_network_top10.graphml")
igraph::write_graph(g, output_graphml, format = "graphml")