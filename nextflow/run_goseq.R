#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(goseq)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag) {
  i <- which(args == flag)
  if (length(i) == 0) stop(paste("Falta el parámetro", flag))
  args[i + 1]
}

input_file  <- get_flag("--input")
results_out <- get_flag("--results")
genes_out   <- get_flag("--genes")
pwf_plot    <- get_flag("--pwfplot")
go_plot     <- get_flag("--goplot")
gtf_file    <- get_flag("--gtf")

# Leer resultados de DESeq2
res <- read.table(input_file, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
res <- res[!is.na(res$padj), ]

# Limpiar versión del identificador ENST
res$transcript_id <- sub("\\..*$", "", res$gene_id)
res$enst_original <- res$gene_id
res$gene_id       <- NULL

# Leer GTF y extraer tanto longitudes de gen como mapeo ENST->ENSG
message("Leyendo GTF...")
con <- gzcon(file(gtf_file, "rb"))
gtf_lines <- readLines(con)
close(con)

gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]

# Extraer atributo genérico
get_attr <- function(attr, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regexpr(pattern, attr, perl = TRUE)
  if (m == -1) return(NA_character_)
  val <- regmatches(attr, m)
  gsub(paste0(key, ' "|"'), '', val)
}

# 1. Longitudes génicas desde líneas de tipo 'gene'
message("Extrayendo longitudes génicas...")
gene_lines <- gtf_lines[grepl("\tgene\t", gtf_lines)]
message(paste("Líneas de gen:", length(gene_lines)))

length_df <- do.call(rbind, lapply(gene_lines, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  if (length(fields) < 9) return(NULL)
  gene_id <- get_attr(fields[9], "gene_id")
  if (is.na(gene_id)) return(NULL)
  gene_id <- sub("\\..*$", "", gene_id)
  data.frame(gene_id = gene_id,
             length  = as.numeric(fields[5]) - as.numeric(fields[4]) + 1,
             stringsAsFactors = FALSE)
}))
length_df <- length_df[!duplicated(length_df$gene_id), ]
message(paste("Genes con longitud:", nrow(length_df)))

# 2. Mapeo ENST->ENSG desde líneas de tipo 'transcript'
message("Extrayendo mapeo ENST->ENSG...")
tx_lines <- gtf_lines[grepl("\ttranscript\t", gtf_lines)]
message(paste("Líneas de transcrito:", length(tx_lines)))

map_df <- do.call(rbind, lapply(tx_lines, function(line) {
  fields <- strsplit(line, "\t")[[1]]
  if (length(fields) < 9) return(NULL)
  gene_id <- get_attr(fields[9], "gene_id")
  tx_id   <- get_attr(fields[9], "transcript_id")
  if (is.na(gene_id) || is.na(tx_id)) return(NULL)
  data.frame(
    transcript_id = sub("\\..*$", "", tx_id),
    gene_id       = sub("\\..*$", "", gene_id),
    stringsAsFactors = FALSE
  )
}))
map_df <- map_df[!duplicated(map_df$transcript_id), ]
message(paste("Transcritos mapeados:", nrow(map_df)))

# Unir DESeq2 con mapeo ENST->ENSG
merged <- merge(res, map_df, by = "transcript_id")
message(paste("Genes tras mapeo ENST->ENSG:", nrow(merged)))

# Convertir a data.frame plano por seguridad
merged <- as.data.frame(merged)
merged$gene_id <- as.character(merged$gene_id)
merged$padj    <- as.numeric(merged$padj)

# Si varios transcritos mapean al mismo gen, conservar el de menor padj
merged <- merged[order(merged$gene_id, is.na(merged$padj), merged$padj), ]
merged <- merged[!duplicated(merged$gene_id), ]

# Vector binario 0/1
merged$de_flag <- ifelse(merged$padj < 0.05, 1, 0)
gene_vector <- merged$de_flag
names(gene_vector) <- merged$gene_id
message(paste("Genes DE (padj < 0.05):", sum(gene_vector)))

# Alinear longitudes
length_vec        <- length_df$length
names(length_vec) <- length_df$gene_id
length_vec        <- length_vec[names(gene_vector)]

valid        <- !is.na(length_vec)
gene_vector  <- gene_vector[valid]
length_vec   <- length_vec[valid]
message(paste("Genes con longitud disponible:", sum(valid)))

if (sum(valid) == 0) stop("No se pudieron alinear genes con longitudes. Revisar identificadores.")

# PWF
pwf <- nullp(gene_vector, bias.data = length_vec, plot.fit = FALSE)

pdf(pwf_plot)
plotPWF(pwf)
dev.off()

# GOseq
go_res <- goseq(
  pwf,
  genome    = "hg38",
  id        = "ensGene",
  test.cats = c("GO:BP"),
  method    = "Wallenius"
)

go_res$p_adjust_over_represented  <- p.adjust(go_res$over_represented_pvalue,  method = "BH")
go_res$p_adjust_under_represented <- p.adjust(go_res$under_represented_pvalue, method = "BH")

write.table(go_res, results_out, sep = "\t", quote = FALSE, row.names = FALSE)

# Tabla genes DE por categoría
sig_go      <- go_res[!is.na(go_res$p_adjust_over_represented) &
                       go_res$p_adjust_over_represented < 0.05, ]
de_gene_ids <- names(gene_vector)[as.logical(gene_vector)]
de_genes_per_cat <- getgo(de_gene_ids, "hg38", "ensGene")

if (!is.null(de_genes_per_cat) && length(de_genes_per_cat) > 0) {
  genes_table <- data.frame(
    category = names(de_genes_per_cat),
    de_genes = vapply(de_genes_per_cat,
                      function(x) paste(x, collapse = ","),
                      character(1)),
    stringsAsFactors = FALSE
  )
  genes_table <- genes_table[genes_table$category %in% sig_go$category, ]
} else {
  genes_table <- data.frame(category = character(0), de_genes = character(0))
}

write.table(genes_table, genes_out, sep = "\t", quote = FALSE, row.names = FALSE)

# Gráfico de burbujas
if (nrow(sig_go) > 0) {
  top_go        <- sig_go[order(sig_go$p_adjust_over_represented), ]
  top_go        <- head(top_go, 15)
  top_go$percDE <- 100 * top_go$numDEInCat / top_go$numInCat
  top_go$term   <- ifelse(is.na(top_go$term) | top_go$term == "",
                           top_go$category, top_go$term)

  p <- ggplot(top_go, aes(x = percDE, y = reorder(term, percDE))) +
    geom_point(aes(size = numDEInCat, color = p_adjust_over_represented)) +
    scale_color_gradient(low = "steelblue", high = "red") +
    labs(
      title    = "Top over-represented GO:BP categories",
      subtitle = "Wallenius method, BH correction",
      x        = "% DE genes in category",
      y        = "GO term",
      color    = "Adj. p-value",
      size     = "DE count"
    ) +
    theme_minimal(base_size = 11)

  ggsave(go_plot, plot = p, width = 9, height = 7)
} else {
  message("No significant GO categories found. Skipping bubble plot.")
  pdf(go_plot)
  plot.new()
  text(0.5, 0.5, "No significant GO:BP categories found", cex = 1.2)
  dev.off()
}
