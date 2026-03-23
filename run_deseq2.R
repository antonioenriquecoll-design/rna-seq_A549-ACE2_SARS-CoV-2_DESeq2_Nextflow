#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 0) return(character(0))
  args[i + 1]
}

infected_arg <- get_flag("--infected")
mock_arg     <- get_flag("--mock")
outfile      <- get_flag("--out")
plots        <- get_flag("--plots")

# IMPORTANTE: dividir por espacios las listas que llegan desde Nextflow
infected <- unlist(strsplit(infected_arg, " +"))
mock     <- unlist(strsplit(mock_arg,     " +"))

files <- c(infected, mock)

# Leer TSV (gene_id, NumReads)
count_list <- lapply(files, function(f) {
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
})

# merge iterativo por gene_id
counts <- Reduce(function(x, y) merge(x, y, by = "gene_id"), count_list)
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

# coldata según el orden de las columnas (infected primero, mock después)
coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(c(rep("infected", length(infected)),
                       rep("mock",     length(mock))))
)

dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                              colData   = coldata,
                              design    = ~ condition)

dds <- DESeq(dds)
res <- results(dds)

write.table(as.data.frame(res), outfile, sep = "\t", quote = FALSE, col.names = NA)

pdf(plots)
plotMA(res, main = "MA-plot")
plotPCA(vst(dds), intgroup = "condition")
dev.off()