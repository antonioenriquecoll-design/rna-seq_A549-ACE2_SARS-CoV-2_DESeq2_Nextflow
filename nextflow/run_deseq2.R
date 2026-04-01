#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 0) stop(paste("Falta el parámetro", flag))
  args[i + 1]
}

infected_arg <- get_flag("--infected")
mock_arg     <- get_flag("--mock")
outfile      <- get_flag("--out")
plots        <- get_flag("--plots")

infected <- unlist(strsplit(infected_arg, " +"))
mock     <- unlist(strsplit(mock_arg,     " +"))
files    <- c(infected, mock)

# Leer cada tabla y renombrar NumReads con el nombre del archivo
count_list <- lapply(files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
  sample_name <- sub("\\.tsv$", "", basename(f))
  colnames(df)[colnames(df) == "NumReads"] <- sample_name
  df
})

# Unir todas las tablas por gene_id
counts <- Reduce(function(x, y) merge(x, y, by = "gene_id"), count_list)
rownames(counts) <- counts$gene_id
counts$gene_id   <- NULL

# Metadatos de condición
coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(
    c(rep("infected", length(infected)),
      rep("mock",     length(mock))),
    levels = c("mock", "infected")
  )
)

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds, fitType = "parametric")
res <- results(dds, alpha = 0.1)

# Tabla de resultados
res_df         <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange",
                     "lfcSE", "stat", "pvalue", "padj")]

write.table(res_df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

# Plots
vsd <- vst(dds, blind = TRUE)

pdf(plots)

plotMA(res, main = "MA-plot: infected vs mock", alpha = 0.1)

plotPCA(vsd, intgroup = "condition")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
heatmap(sampleDistMatrix,
        main = "Sample distance heatmap",
        col  = colorRampPalette(c("steelblue", "white"))(100))

plotDispEsts(dds, main = "Dispersion estimates")

hist(res$pvalue[res$baseMean > 1],
     breaks = 50,
     col    = "steelblue",
     main   = "Histogram of p-values (baseMean > 1)",
     xlab   = "p-value")

dev.off()
