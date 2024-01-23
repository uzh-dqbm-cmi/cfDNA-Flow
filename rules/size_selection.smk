include: "common.smk"

global WORKDIR
global PARAMDIR
global THREADS
global SIZE_SELECTION_MIN
global SIZE_SELECTION_MAX
global SCRIPTS
global TEMP

rule bed_size_selection:
    input:
        #expand(WORKDIR + "/BED/" + PARAMDIR + "/{sample}.bed", sample = samples_df['sample_name'].values)
        WORKDIR +"/BED/"+PARAMDIR+"/{sample}.bed"
    output:
        #expand(WORKDIR + "/BED/" + PARAMDIR + "/szsel_"+SIZE_SELECTION_MIN+"_"+SIZE_SELECTION_MAX+"/{sample}.bed", sample = samples_df['sample_name'].values)
        WORKDIR + "/BED/" + PARAMDIR + "/szsel_{min}_{max}/{sample}.bed"
    wildcard_constraints:
        min="\d+",
        max="\d+",
        sample="[^/]+"
    params:
        min = "{min}",   # SIZE_SELECTION_MIN,
        max = "{max}",   # SIZE_SELECTION_MAX,
        n_cores = THREADS,
        n_partitions = THREADS,
        temp = TEMP
    shell: """
        start_time=$(date +%s)
        {SCRIPTS}/bedpe_size_selection.sh -i {input} --min {params.min} --max {params.max} -o {output} -p {params.n_cores} --temp {params.temp}
        end_time=$(date +%s)
        elapsed=$(( end_time - start_time ))
        eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H:%M:%S,%3N ')"
        """
#    shell: """
#        python {SCRIPTS}/bedpe_size_selection.py -i {input} -min {params.min} -max {params.max} -o {output} -n {params.n_cores} -p {params.n_partitions}
#        """


rule bam_size_selection:
    input:
        WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bam"
    output:
        bam=WORKDIR + "/BAM/" + PARAMDIR + "/szsel_{min}_{max}/{sample}.sortByCoord.bam",
        bai=WORKDIR + "/BAM/" + PARAMDIR + "/szsel_{min}_{max}/{sample}.sortByCoord.bai"
    wildcard_constraints:
        min="\d+",
        max="\d+",
        sample="[^/]+"
    params:
        min = "{min}",   # SIZE_SELECTION_MIN,
        max = "{max}",   # SIZE_SELECTION_MAX,
        temp = TEMP
    threads: THREADS
    shell: """
        samtools view -h {input} | awk 'substr($0,1,1)=="@" || ($9>= {params.min} && $9<=220) || ($9<=-{params.min} && $9>=-{params.max})' | samtools view -b > {output.bam}
        samtools index -@ {threads} {output.bam} {output.bai}
        """
