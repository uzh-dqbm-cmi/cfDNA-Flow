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
global TEMP
global samples
global control_samples
global roi_names


controls_sortByCoord = " ".join([f"{item}.sortByCoord" for item in control_samples[0].tolist()])

rule LIQUORICE:
    input:
        bam=WORKDIR + "/BAM/" + PARAMDIR + "/" + SIZE_SELECTION_SUB_FOLDER + "{sample}.sortByCoord.bam"
    output:
        WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER+"{sample}.sortByCoord/{roi}/fitted_gaussians.pdf"
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
        mkdir -p {params.out_dir} 
        cd {params.out_dir} && LIQUORICE --bamfile "{input.bam}" \
            --refgenome_fasta "{params.refgenome_fasta}" \
            --mappability_bigwig "{params.mappability_bigwig}" \
            --bedpathlist {params.bedpathlist} \
            --blacklist {params.blacklist} \
            --n_cpus {params.n_cpus}
        """


rule LIQUORICE_summary:
    input:
        expand(WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER+"{sample}.sortByCoord/{roi}/fitted_gaussians.pdf",
            sample = samples['sample_name'].values, roi=roi_names)
    output:
        WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER +"summary_across_samples_and_ROIS.csv"
    wildcard_constraints:
        min="\d+",
        max="\d+",
        sample="[^/]+"
    params:
        controls = controls_sortByCoord,
        out_dir= WORKDIR + "/feature/" + PARAMDIR + "/liquorice/" + SIZE_SELECTION_SUB_FOLDER
    threads: THREADS
    shell: """
        cd {params.out_dir} && LIQUORICE_summary --control_name_list {params.controls} --dirname {params.out_dir}
        """


