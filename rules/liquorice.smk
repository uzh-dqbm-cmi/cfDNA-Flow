from rules.config import ROI_LIST

include: "common.smk"

global WORKDIR
global PARAMDIR
global THREADS
global SCRIPTS
global SIZE_SELECTION_SUB_FOLDER
global REFGENOME
global REFGENOME_MAPPABILITY
global BLACKLIST
global ROI_LIST
global CONTROL_SAMPLES
global TEMP

rule LIQUORICE:
    input:
        bai = WORKDIR + "/BAM/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.sortByCoord.bam.bai",
        bam=WORKDIR + "/BAM/" + PARAMDIR + "/" + SIZE_SELECTION_SUB_FOLDER + "{sample}.sortByCoord.bam"
    output:
        expand(WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER+"/{sample}/*/fitted_gaussians.pdf",
            sample = samples_df['sample_name'].values)
    workdir:
        WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER
    wildcard_constraints:
        min="\d+",
        max="\d+",
        sample="[^/]+"
    params:
        refgenome_fasta = REFGENOME,
        mappability_bigwig = REFGENOME_MAPPABILITY,
        bedpathlist = ROI_LIST,
        blacklist = BLACKLIST,
        n_cpus = THREADS,
        out_dir = WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER
    shell: """
        LIQUORICE --bamfile "{input.bam}" \
            --refgenome_fasta "{refgenome_fasta}" \
            --mappability_bigwig "{mappability_bigwig}" \
            --bedpathlist {bedpathlist} \
            --blacklist {blacklist} \
            --n_cpus {n_cpus}
        """


rule LIQUORICE_summary:
    input:
        expand(WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER+"/{sample}/*/fitted_gaussians.pdf",
            sample = samples_df['sample_name'].values)
    output:
        WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER +"/summary_across_samples_and_ROIs.csv"
    workdir:
        WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER
    wildcard_constraints:
        min="\d+",
        max="\d+",
        sample="[^/]+"
    params:
        temp = TEMP,
        controls = CONTROL_SAMPLES,
        out_dir= WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER
    threads: THREADS
    shell: """
        LIQUORICE_summary --tmpdir {temp} --control_name_list {controls}
        """


