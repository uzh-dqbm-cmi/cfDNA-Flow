# cfDNA-Flow

## Running the pipeline
`snakemake --configfile path_to_config/config.yaml -j16 do_all`


## Output
`/BAM`

**sample.sorted.bam**: BAM files sorted by coordinates (temporary file)

**sample.markdup.bam**: BAM files sorted by coordinates, duplicates marked

**sample.sortedByCoord.bam**: BAM files sorted by coordinate, duplicates removed, final alignment file output, ready for downstream analysis
