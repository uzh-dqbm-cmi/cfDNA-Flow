# cfDNA preprocessing pipeline

## Running the pipeline
`snakemake --configfile /cluster/dataset/medinfmk/pipeline-cfDNA/results/0memhg19decoy/0memhg19decoyTrue40-020210115.yaml -j16 preprocess`


## Output
`/BAM`

**sample.sorted.bam**: BAM files sorted by coordinates (temporary file)

**sample.markdup.bam**: BAM files sorted by coordinates, duplicates marked

**sample.sortedByCoord.bam**: BAM files sorted by coordinate, duplicates removed, final alignment file output, ready for downstream analysis
