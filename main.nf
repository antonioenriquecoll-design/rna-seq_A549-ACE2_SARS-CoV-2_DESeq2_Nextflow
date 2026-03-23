nextflow.enable.dsl=2

params.reads         = params.reads         ?: 'data/*.fastq.gz'
params.transcriptome = params.transcriptome ?: 'ref/Homo_sapiens.GRCh38.cdna.all.fa.gz'
params.annotation    = params.annotation    ?: 'ref/Homo_sapiens.GRCh38.109.gtf.gz'

/* -----------------------------
   CHANNELS
   ----------------------------- */
reads           = Channel.fromPath(params.reads)
transcriptomeCh = Channel.value(file(params.transcriptome))
annotationCh    = Channel.value(file(params.annotation))

/* Duplicar reads para QC y trimming */
reads.into { qc_reads; trim_reads }

/* -----------------------------
   PROCESSES
   ----------------------------- */

process FASTQC {
  tag { fastq.simpleName }
  input:  path fastq from qc_reads
  output: path "fastqc_${fastq.simpleName}"
  script:
  """
  fastqc ${fastq} -o fastqc_${fastq.simpleName}
  """
}

process TRIMMING {
  tag { fastq.simpleName }
  input:  path fastq from trim_reads
  output: path "trimmed_${fastq.simpleName}.fq.gz"
  script:
  """
  cutadapt -q 20 -m 20 -a AGATCGGAAGAGC \\
           -o trimmed_${fastq.simpleName}.fq.gz \\
           ${fastq}
  """
}

/* Construir índice Salmon (una vez) */
process BUILD_SALMON_INDEX {
  tag 'salmon_index'
  input:  path transcriptome from transcriptomeCh
  output: path "salmon_index"
  script:
  """
  salmon index -t ${transcriptome} -i salmon_index
  """
}

process SALMON {
  tag { trimmed.simpleName }
  input:
    path trimmed
    path index
    path ann
  output:
    path "quant_${trimmed.simpleName}"
  script:

  """
  salmon quant -i ${index} \\
      -l A \\
      -r ${trimmed} \\
      --validateMappings \\
      --seqBias --gcBias \\
      --fldMean 200 --fldSD 80 \\
      --geneMap ${ann} \\
      -o quant_${trimmed.simpleName}
  """
}

process EXTRACT_COUNTS {
  tag { quantdir.simpleName }
  input:  path quantdir
  output: path "counts_${quantdir.simpleName}.tsv"
  script:
  """
  awk 'BEGIN{FS=OFS="\\t"} NR==1{print "gene_id","NumReads"; next}{print \$1, \$5}' \\
      ${quantdir}/quant.genes.sf > counts_${quantdir.simpleName}.tsv
  """
}

process ROUND_COUNTS {
  tag { counts.simpleName }
  input:  path counts
  output: path "counts_int_${counts.simpleName}.tsv"
  script:
  """
  awk 'BEGIN {FS=OFS="\\t"} NR==1{print \$0; next}{\$2=int(\$2+0.5); print}' \\
      ${counts} > counts_int_${counts.simpleName}.tsv
  """
}

process DESEQ2 {
  publishDir "results", mode: 'copy'
  tag 'DESeq2'

  input:
    path infected_files
    path mock_files

  output:
    file "deseq2_results.tsv"
    file "plots.pdf"

  script:
  """
  Rscript run_deseq2.R \\
      --infected ${infected_files.join(" ")} \\
      --mock     ${mock_files.join(" ")} \\
      --out      deseq2_results.tsv \\
      --plots    plots.pdf
  """
}

/* -----------------------------
   WORKFLOW
   ----------------------------- */


workflow {

  /* QC (se ejecuta en paralelo, sin usar su salida) */
  FASTQC(qc_reads)

  /* trimming */
  trimmed = TRIMMING(trim_reads)

  /* índice + cuantificación */
  index = BUILD_SALMON_INDEX(transcriptomeCh)
  quant = SALMON(trimmed, index.out, annotationCh)

  /* counts */
  raw_counts   = EXTRACT_COUNTS(quant)
  clean_counts = ROUND_COUNTS(raw_counts)

  /* Selección por nombre de fichero: 3 y 3 */
  infected_list = clean_counts
      .filter{ it.name.contains("SRR11517720") || it.name.contains("SRR11517721") || it.name.contains("SRR11517722") }
      .collect()

  mock_list = clean_counts
      .filter{ it.name.contains("SRR11517724") || it.name.contains("SRR11517725") || it.name.contains("SRR11517726") }
      .collect()

  /* Lanzar DESeq2 con las LISTAS de archivos */
  DESEQ2(infected_list, mock_list)
}