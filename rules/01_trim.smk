include: "common.smk"

global WORKDIR
global PARAMDIR
global THREADS
global samples

def getFastPair(wildcards):
    loc = samples.loc[wildcards.sample]
    # print(f'in getFastPair: {config["sampleDir"]}, {wildcards.sample}')
    # print(f'loc: {loc}')
    forward = f'{config["sampleDir"]}{loc["fwd"]}'
    reverse = f'{config["sampleDir"]}{loc["rvr"]}'
    return forward, reverse

rule trimmed:
    input:  WORKDIR+"/QC/trimmed_reads/multiqc_report.html"

rule trim:
    input:
        pair = getFastPair
    params:
        q = config["q"],
        Q = config["Q"],
        l = config["l"],
        prefix = WORKDIR+"/trimmed_reads/{sample}"
    log:
        WORKDIR+"/logs/trimmed_reads/{sample}.log"
    output:
        temp(WORKDIR+"/trimmed_reads/{sample}-trimmed-pair1.fastq"),
        temp(WORKDIR+"/trimmed_reads/{sample}-trimmed-pair2.fastq")
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.trim.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    shell:"""
        mkdir -p {WORKDIR}/trimmed_reads
        skewer -Q {params.Q} \
            -q {params.q} \
            -l {params.l} \
            {input} \
            -o {params.prefix} &> {log}
        """

rule fastqc:
    input:
        WORKDIR+"/trimmed_reads/{sample}-trimmed-pair1.fastq",
        WORKDIR+"/trimmed_reads/{sample}-trimmed-pair2.fastq"
    params:
        #outdir = WORKDIR+"/QC/trimmed_reads"
        threads= THREADS,
        outdir= WORKDIR+"/QC/" + PARAMDIR
    output:
        # WORKDIR+"/QC/trimmed_reads/{sample}-trimmed-pair1_fastqc.html",
        # WORKDIR+"/QC/trimmed_reads/{sample}-trimmed-pair2_fastqc.html"
        WORKDIR + "/QC/" + PARAMDIR + "/{sample}-trimmed-pair1_fastqc.html",
        WORKDIR + "/QC/" + PARAMDIR + "/{sample}-trimmed-pair2_fastqc.html"
    wildcard_constraints:
        sample="[^/]+"
    shell: 
        "fastqc -t {params.threads} -o {params.outdir} {input}"

# rule multiqcOnFastqc:
#     input:
#         #expand(WORKDIR+"/QC/trimmed_reads/{sample}-trimmed-pair1_fastqc.html", sample=samples['sample_name'].values),
#         #expand(WORKDIR+"/QC/trimmed_reads/{sample}-trimmed-pair2_fastqc.html", sample=samples['sample_name'].values)
#         expand(WORKDIR + "/QC/{sample}-trimmed-pair1_fastqc.html", sample=samples['sample_name'].values),
#         expand(WORKDIR + "/QC/{sample}-trimmed-pair2_fastqc.html", sample=samples['sample_name'].values)
#     output:
#         #WORKDIR+"/QC/trimmed_reads/"+PARAMDIR+"/multiqc_report.html"
#         WORKDIR + "/QC/" + PARAMDIR + "/multiqc_report.html"
#     shell:
#         #"multiqc {WORKDIR}/QC/trimmed_reads/{PARAMDIR} -o {output}"
#         "multiqc {WORKDIR}/QC/{PARAMDIR} -o {output}"
