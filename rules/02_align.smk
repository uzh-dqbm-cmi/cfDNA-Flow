include: "common.smk"
include: "01_trim.smk"

global PARAMDIR
global REFGENOME
global WORKDIR
global THREADS

rule bwaIndex:
    input: REFGENOME
    output: REFGENOME + ".bwt"
    shell: "bwa index {input}"

rule bwaMem:
    input:
        ref = REFGENOME,
        frwd = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair1.fastq",
        rvr = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair2.fastq",
        index = REFGENOME + ".bwt"
    output:
        bam = "{WORKDIR}/BAM/{sample}.mem_sorted.bam",
        bai = "{WORKDIR}/BAM/{sample}.mem_sorted.bam.bai"
    benchmark:
        "{WORKDIR}/benchmarks/"+PARAMDIR+"/{sample}.bwaMem.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        flag = config['SAM_flag'],
        mapq = config['MAPQ'],
        cigar = config['cigar_flag'],
        tag = config['tag'],
        threads = config['threads']
    threads: THREADS
    shell:"""
        bwa mem -t {params.threads} -M {input.ref} {input.frwd} {input.rvr}  | \
            samtools sort --reference {input.ref} -l 9 -@ {params.threads} -O bam -o {output.bam} -T {output.bam}.tmp
        samtools index -@ {params.threads} {output.bam}
        """

rule bwaMem2:
    input:
        ref = REFGENOME,
        frwd = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair1.fastq",
        rvr = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair2.fastq",
        index = REFGENOME + ".bwt"
    output:
        bam = "{WORKDIR}/BAM/{sample}.mem2_sorted.bam",
        bai = "{WORKDIR}/BAM/{sample}.mem2_sorted.bam.bai"
    benchmark:
        "{WORKDIR}/benchmarks/"+PARAMDIR+"/{sample}.bwaMem2.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        flag = config['SAM_flag'],
        mapq = config['MAPQ'],
        cigar = config['cigar_flag'],
        tag = config['tag'],
        threads = config['threads']
    threads: THREADS
    shell:"""
        bwa-mem2 mem -t {params.threads} -M {input.ref} {input.frwd} {input.rvr}  | \
            samtools sort --reference {input.ref} -l 9 -@ {params.threads} -O bam -o {output.bam} -T {output.bam}.tmp
        samtools index -@ {params.threads} {output.bam}
        """

rule bwaAln:
    input:
        ref = REFGENOME,
        frwd = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair1.fastq",
        rvr = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair2.fastq"
    output:
        frwd_sai = "{WORKDIR}/BAM/{sample}-pair1.sai",
        rvr_sai = "{WORKDIR}/BAM/{sample}-pair2.sai"
    benchmark:
        "{WORKDIR}/benchmarks/"+PARAMDIR+"/{sample}.bwaAln.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        threads = config['threads']
    shell:
        "bwa aln -t {params.threads} {input.ref} {input.frwd} > {output.frwd_sai}; "
        "bwa aln -t {params.threads} {input.ref} {input.rvr} > {output.rvr_sai}"

rule bwaSampe:
    input:
        ref = REFGENOME,
        frwd = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair1.fastq",
        rvr = "{WORKDIR}/trimmed_reads/{sample}-trimmed-pair2.fastq",
        frwd_sai = "{WORKDIR}/BAM/{sample}-pair1.sai",
        rvr_sai = "{WORKDIR}/BAM/{sample}-pair2.sai"
    output:
        bam = "{WORKDIR}/BAM/{sample}.aln_sorted.bam",
        bai = "{WORKDIR}/BAM/{sample}.aln_sorted.bam.bai"
    benchmark:
        "{WORKDIR}/benchmarks/"+PARAMDIR+"/{sample}.bwaSampe.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        flag = config['SAM_flag'],
        mapq = config['MAPQ'],
        cigar = config['cigar_flag'],
        tag = config['tag']
    shell:
        "bwa sampe {input.ref} {input.frwd_sai} {input.rvr_sai} {input.frwd} {input.rvr} | "
        "samtools sort --reference {input.ref} "
        "-l 9 -@ 4 -O bam -o {output.bam} -T {output.bam}.tmp; "
        "samtools index -@ 4 {output.bam}"
