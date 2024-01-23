include: "common.smk"
include: "02_align.smk"

global REFGENOME
global PARAMDIR
global WORKDIR
global samples
global aligner

def filterCigar(wildcards):
    if config['cigar_flag'] == 'D|I':
        return "| awk '($0 ~ /^@/) || ($6 !~ /D|I/)' | "
    else:
        return '|'

def filterTag(wildcards):
    if config['tag'] == 'XA':
        return "-x XA "
    else:
        return ''

def filterSAM(wildcards):
    if config['SAM_flag']:
        return f"-F {config['SAM_flag']}"
    else:
        return ''

rule filterBams:
    input:
        bam = aligner,
        ref = REFGENOME
    output:
        filtered = temp(WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sorted.bam"),
        bai = temp(WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sorted.bam.bai")
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.filterBams.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        flag = filterSAM,
        mapq = config['MAPQ'],
        cigar = filterCigar,
        tag = filterTag
    shell:
        "samtools view -@ 4 -h {params.flag} -q {params.mapq} {params.tag} {input.bam} "
        "{params.cigar} "
        "samtools sort --reference {input.ref} -l 9 -@ 4 -O bam -o {output.filtered} -T {output.filtered}.tmp; "
        "samtools index {output.filtered}"

rule markDuplicates:
    input:
        bam = WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sorted.bam",
        bai = WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sorted.bam.bai"
    output:
        bam = temp(WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.markdup.bam"),
        bai = temp(WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.markdup.bai"),
        metrics = WORKDIR+"/QC/"+PARAMDIR+"/{sample}.markdup.metrics"
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.markDuplicates.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        tmpdir = WORKDIR+"/BAM/"+PARAMDIR,
        ram = config['MarkDuplicates']['ram']
    shell:
        "gatk --java-options '-Xmx{params.ram}G' MarkDuplicates "
        "--TMP_DIR {params.tmpdir} "
        "--CREATE_INDEX true "
        "-O {output.bam} -I {input.bam} -M {output.metrics}"

rule picardAlignmentMetrics:
    input: 
        bam = WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.markdup.bam",
        bai = WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.markdup.bai",
        ref = REFGENOME
    output:
        #WORKDIR+"/QC/Picard/"+PARAMDIR+"/{sample}.AlignmentMetrics.txt"
        WORKDIR+"/QC/" + PARAMDIR + "/{sample}.AlignmentMetrics.txt"
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.picardAlignmentMetrics.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "gatk CollectAlignmentSummaryMetrics "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output}"

rule picardCollectWgsMetrics:
    input: 
        bam = WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.markdup.bam",
        bai = WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.markdup.bai",
        ref = REFGENOME
    output:
        #WORKDIR+"/QC/Picard/"+PARAMDIR+"/{sample}.collect_wgs_metrics.txt"
        WORKDIR + "/QC/" + PARAMDIR + "/{sample}.collect_wgs_metrics.txt"
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.picardCollectWgsMetrics.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "gatk CollectWgsMetrics "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output}"


rule multiqc:
    input:
        # expand(WORKDIR+"/QC/Picard/"+PARAMDIR+"/{sample}.AlignmentMetrics.txt", sample=samples['sample_name'].values),
        # expand(WORKDIR+"/QC/Picard/"+PARAMDIR+"/{sample}.collect_wgs_metrics.txt", sample=samples['sample_name'].values)
        expand(WORKDIR + "/QC/" + PARAMDIR + "/{sample}.AlignmentMetrics.txt",sample=samples['sample_name'].values),
        expand(WORKDIR + "/QC/" + PARAMDIR + "/{sample}.collect_wgs_metrics.txt",sample=samples['sample_name'].values),
        expand(WORKDIR + "/QC/" + PARAMDIR + "/{sample}-trimmed-pair1_fastqc.html", sample=samples['sample_name'].values),
        expand(WORKDIR + "/QC/" + PARAMDIR + "/{sample}-trimmed-pair2_fastqc.html", sample=samples['sample_name'].values)
    output:
        #"{WORKDIR}/QC/Picard/{PARAMDIR}",
        # WORKDIR+"/QC/Picard/"+PARAMDIR+"/multiqc_report.html",
        directory("{WORKDIR}/QC/{PARAMDIR}/multiqc_data"),
        #WORKDIR + "/QC/" + PARAMDIR + "/multiqc_data/multiqc_report.html"
    shell:
        #"multiqc {WORKDIR}/QC/Picard/{PARAMDIR} -o {output}"
        "multiqc {WORKDIR}/QC/{PARAMDIR} -o {output}"

