include: "common.smk"

global WORKDIR
global PARAMDIR
global THREADS
global SIZE_SELECTION_SUB_FOLDER
global SCRIPTS
global CONTROL_SAMPLES
global SAMPLES
global BIN_SIZE
global GENOME
global CNA_DO_CONTROL_REFS
global RC_NORMALIZATION
global hash_it
global control_RefName
global control_samples
global treated_samples_RefName
global treated_samples_names

rule cal_blacklist: # 01
    input:
        expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"{control_sample}_point.csv", control_sample = control_samples[0].values)
        # expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"{control_sample}.bed", control_sample = control_samples[0].values),
    output:
        #expand(WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.bed", sample = samples_df['sample_name'].values)
        WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"control_blacklist_"+BIN_SIZE+"Kbp_norm_TRUE.RDS"
    params:
        bed_sample_dir = WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER,
        control_samples = CONTROL_SAMPLES,   # file name: every line is sample name
        bin_size = BIN_SIZE,
        RCNormalization = RC_NORMALIZATION
    shell: """
        Rscript --vanilla {SCRIPTS}/tMAD/01_cal_blacklist.R {params.bed_sample_dir} {params.control_samples} {params.bin_size} {params.RCNormalization}
        """

rule cal_RefSample: # 02
    input:
        expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"{control_sample}_point.csv", control_sample = control_samples[0].values)
        # expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"{control_sample}.bed", control_sample = control_samples[0].values)
    output:
        WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+control_RefName+"_norm"+RC_NORMALIZATION+"_"+BIN_SIZE+"kbp_"+GENOME+"_count.csv"
        # there are more but skipped here: all combinations excluding one by one the control samples
    params:
        bed_sample_dir = WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER,
        control_samples = CONTROL_SAMPLES,   # file name: every line is sample name
        bin_size = BIN_SIZE,
        RCNormalization = RC_NORMALIZATION, # "FALSE",  # strange here is FALSE in the script  not as previous step 01_cal_blacklist
        genome = GENOME,
        do_control_refs = CNA_DO_CONTROL_REFS
    shell: """
        Rscript --vanilla {SCRIPTS}/tMAD/02_cal_RefSample.R {params.bed_sample_dir} {params.control_samples} {params.bin_size} {params.RCNormalization} {params.genome} {params.do_control_refs}
        """


rule cal_t_MAD_forall:  # 05
    input:
        expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"{control_sample}_point.csv", control_sample = control_samples[0].values)
    output:
        WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/tMAD_results.tsv"
    params:
        bed_sample_dir = WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER,
        control_samples = CONTROL_SAMPLES,   # file name: every line is sample name
        bin_size = BIN_SIZE,
        RCNormalization = RC_NORMALIZATION,  # strange here is FALSE in the script  not as previous step 01_cal_blacklist
        genome = GENOME,
        do_control_refs = CNA_DO_CONTROL_REFS,
        samples = SAMPLES,
        output_dir = WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/"
    shell: """
        Rscript --vanilla {SCRIPTS}/tMAD/05_cal_t_MAD_forall.R {params.bed_sample_dir} {params.control_samples} {params.bin_size} {params.RCNormalization} {params.genome} {params.do_control_refs} {params.samples} {params.output_dir}
        """

rule visualising_t_MAD_forall:  # 06_visualizing_forall
    input:
        WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/CNAData_"+BIN_SIZE+"_control_none_samples_"+treated_samples_RefName+".RDS",
        WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/CNAData_"+BIN_SIZE+"_control_"+control_RefName+"_samples_"+treated_samples_RefName+".RDS"
    output:
        expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/CNV_sample_{control_sample}.pdf", control_sample = control_samples[0].values),
        expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/CNV_sample_{treated_sample}.pdf", treated_sample = treated_samples_names)
        # and more output files for comparing with a reference
    params:
        bed_sample_dir = WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER,
        control_samples = CONTROL_SAMPLES,
        bin_size = BIN_SIZE,
        do_control_refs = CNA_DO_CONTROL_REFS,
        samples = SAMPLES,
        output_dir = WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/"
    shell: """
        Rscript --vanilla {SCRIPTS}/tMAD/06_visualizing_forall.R {params.bed_sample_dir} {params.control_samples} {params.bin_size} {params.do_control_refs} {params.samples} {params.output_dir}
        """

