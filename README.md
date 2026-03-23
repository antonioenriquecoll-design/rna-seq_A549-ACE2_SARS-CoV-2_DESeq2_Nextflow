# RNA‑seq Differential Expression: A549‑ACE2 (SARS‑CoV‑2)
# Salmon + DESeq2 + Nextflow + Galaxy Europe (https://usegalaxy.eu/)

Repositorio del análisis de expresión diferencial (infected vs mock) con 3 réplicas por condición. Incluye un workflow reproducible en Nextflow, un script R para DESeq2, y los resultados/figuras generados en Galaxy Europe.

------------------------------------------------------------
1. DATOS UTILIZADOS
------------------------------------------------------------

Accesiones SRA analizadas:

Infected:
- SRR11517720
- SRR11517721
- SRR11517722

Mock:
- SRR11517724
- SRR11517725
- SRR11517726

Los archivos FASTQ no se incluyen en este repositorio por razones de tamaño. Deben colocarse en:

data/*.fastq.gz

------------------------------------------------------------
2. ARCHIVOS DE REFERENCIA
------------------------------------------------------------

Descargar a la carpeta ref/:

Homo_sapiens.GRCh38.cdna.all.fa.gz
Homo_sapiens.GRCh38.109.gtf.gz

Comandos de descarga:

wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

------------------------------------------------------------
3. ANÁLISIS REALIZADO EN GALAXY EUROPE (https://usegalaxy.eu/)
------------------------------------------------------------

Todo el análisis también fue ejecutado en Galaxy Europe (cutadapt, Salmon, unión de matrices, DESeq2, volcano plot, PCA, MA-plot, VST, heatmaps, etc.).

La carpeta results/ contiene todos los archivos generados en Galaxy Europe.
La carpeta plots/ contiene todas las figuras utilizadas en el informe.

Historial público de Galaxy (para evaluación):
https://usegalaxy.eu/u/antonioenrique/h/rna-seq-sars-cov-2-a549-ace2

------------------------------------------------------------
4. ENTORNO (CONDA)
------------------------------------------------------------

Antes de ejecutar el workflow Nextflow, debe crearse y activarse el entorno Conda.

IMPORTANTE:
Estos comandos se ejecutan en la TERMINAL. NO van dentro de main.nf, ni dentro de run_deseq2.R, ni dentro de ningún archivo.

Comandos:

mamba env create -f envs/environment.yml
conda activate rna_seq_nf

Esto instalará Nextflow, Salmon, Cutadapt, FastQC, R, DESeq2 y las dependencias necesarias.

------------------------------------------------------------
5. EJECUCIÓN DEL WORKFLOW NEXTFLOW
------------------------------------------------------------

El archivo principal del pipeline es:

main.nf

Para ejecutarlo, utilice este comando DESDE LA TERMINAL:

nextflow run main.nf \
  --reads "data/*.fastq.gz" \
  --transcriptome "ref/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
  --annotation "ref/Homo_sapiens.GRCh38.109.gtf.gz" \
  -with-report -with-trace -with-timeline -resume

IMPORTANTE:
Este comando NO va dentro de main.nf. Siempre debe ejecutarse desde el terminal del sistema.

Los resultados generados por Nextflow se guardarán en:

results/

------------------------------------------------------------
6. ESTRUCTURA DEL REPOSITORIO
------------------------------------------------------------

El repositorio debe tener esta estructura:

.
├── README.md                 (este archivo)
├── main.nf                   (workflow Nextflow)
├── run_deseq2.R              (script R para DESeq2)
├── nextflow.config           (configuración adicional)
├── envs/
│   └── environment.yml       (entorno Conda)
├── data/
│   └── README.txt            (explica dónde colocar los FASTQ)
├── ref/
│   ├── README.txt            (explica cómo descargar las referencias)
│   ├── Homo_sapiens.GRCh38.cdna.all.fa.gz
│   └── Homo_sapiens.GRCh38.109.gtf.gz
├── results/                  (resultados obtenidos en Galaxy)
└── plots/                    (todas las figuras finales)

------------------------------------------------------------
7. AUTORÍA
------------------------------------------------------------

Antonio Enrique Coll Meseguer
Belén Salar Benito
Máster en Análisis de Datos Ómicos (UMU)
Curso 2025–2026