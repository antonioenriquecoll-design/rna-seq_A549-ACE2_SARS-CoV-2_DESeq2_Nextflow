En esta carpeta deben colocarse los archivos FASTQ utilizados como entrada del workflow.

Los archivos NO se incluyen en este repositorio (son muy grandes).

El workflow Nextflow detectará automáticamente cualquier archivo que siga el patrón:

    data/*.fastq.gz

Ejemplo de descarga desde SRA (opcional):

    fasterq-dump SRR11517720
    gzip SRR11517720.fastq

Repetir para todos los accesiones utilizados en el análisis:

Infected:
  - SRR11517720
  - SRR11517721
  - SRR11517722

Mock:
  - SRR11517724
  - SRR11517725
  - SRR11517726