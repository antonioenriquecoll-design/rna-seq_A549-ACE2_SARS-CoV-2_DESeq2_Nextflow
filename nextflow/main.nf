#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads         = 'data/*.fastq.gz'
params.transcriptome = 'ref/Homo_sapiens.GRCh38.cdna.all.fa.gz'
params.annotation    = 'ref/Homo_sapiens.GRCh38.109.gtf.gz'
params.deseq2_script = 'run_deseq2.R'
params.goseq_script  = 'run_goseq.R'

workflow {

    reads_ch         = Channel.fromPath(params.reads)
    transcriptome_ch = Channel.value(file(params.transcriptome))
    annotation_ch    = Channel.value(file(params.annotation))
    deseq2_script_ch = Channel.value(file(params.deseq2_script))
    goseq_script_ch  = Channel.value(file(params.goseq_script))

    // Control de calidad inicial
    FASTQC(reads_ch)

    // Trimming
    trimmed_ch = TRIMMING(reads_ch)

    // Índice de Salmon
    salmon_index_ch = BUILD_SALMON_INDEX(transcriptome_ch)

    // Cuantificación con Salmon
    quants_ch = SALMON(trimmed_ch, salmon_index_ch, annotation_ch)

    // Extracción y redondeo de conteos
    raw_counts_ch = EXTRACT_COUNTS(quants_ch)
    int_counts_ch = ROUND_COUNTS(raw_counts_ch)

    // Separar muestras por condición
    infected_ch = int_counts_ch
        .filter { it.name =~ /SRR11517720|SRR11517721|SRR11517722/ }
        .collect()

    mock_ch = int_counts_ch
        .filter { it.name =~ /SRR11517724|SRR11517725|SRR11517726/ }
        .collect()

    // DESeq2
    deseq_out_ch = DESEQ2(infected_ch, mock_ch, deseq2_script_ch)

    // GOseq
    GOSEQ(deseq_out_ch[0], goseq_script_ch, annotation_ch)
}

process FASTQC {
    tag { fastq.simpleName }

    input:
    path fastq

    output:
    path "fastqc_${fastq.simpleName}"

    script:
    """
    mkdir -p fastqc_${fastq.simpleName}
    fastqc ${fastq} -o fastqc_${fastq.simpleName}
    """
}

process TRIMMING {
    tag { fastq.simpleName }

    input:
    path fastq

    output:
    path "trimmed_${fastq.simpleName}.fq.gz"

    script:
    """
    cutadapt \
      -q 20 \
      -m 20 \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      -a AGATCGGAAGAGC \
      -o trimmed_${fastq.simpleName}.fq.gz \
      ${fastq}
    """
}

process BUILD_SALMON_INDEX {
    tag "salmon_index"

    input:
    path transcriptome

    output:
    path "salmon_index"

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
    salmon quant \
      -i ${index} \
      -l A \
      -r ${trimmed} \
      --validateMappings \
      --seqBias \
      --gcBias \
      --fldMean 200 \
      --fldSD 80 \
      --geneMap ${ann} \
      -o quant_${trimmed.simpleName}
    """
}

process EXTRACT_COUNTS {
    tag { quantdir.simpleName }

    input:
    path quantdir

    output:
    path "counts_${quantdir.simpleName}.tsv"

    script:
    """
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print "gene_id","NumReads"; next}{print \$1,\$5}' \
      ${quantdir}/quant.genes.sf > counts_${quantdir.simpleName}.tsv
    """
}

process ROUND_COUNTS {
    tag { counts.simpleName }

    input:
    path counts

    output:
    path "counts_int_${counts.simpleName}.tsv"

    script:
    """
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print \$0; next}{\$2=int(\$2+0.5); print}' \
      ${counts} > counts_int_${counts.simpleName}.tsv
    """
}

process DESEQ2 {
    publishDir "results/deseq2", mode: "copy"
    tag "DESeq2"

    input:
    path infected_files
    path mock_files
    path script

    output:
    path "deseq2_results.tsv"
    path "plots.pdf"

    script:
    """
    Rscript ${script} \
      --infected "${infected_files.join(' ')}" \
      --mock "${mock_files.join(' ')}" \
      --out deseq2_results.tsv \
      --plots plots.pdf
    """
}



process GOSEQ {
    publishDir "results/goseq", mode: "copy"
    tag "GOseq"

    input:
    path deseq_results
    path script
    path gtf

    output:
    path "goseq_results.tsv"
    path "goseq_de_genes.tsv"
    path "goseq_pwf.pdf"
    path "goseq_top_bp.pdf"

    script:
    """
    Rscript ${script} \
      --input ${deseq_results} \
      --results goseq_results.tsv \
      --genes goseq_de_genes.tsv \
      --pwfplot goseq_pwf.pdf \
      --goplot goseq_top_bp.pdf \
      --gtf ${gtf}
    """
}
