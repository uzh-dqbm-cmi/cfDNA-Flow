include: "rules/common.smk"
include: "rules/01_trim.smk"
include: "rules/02_align.smk"
include: "rules/03_clean_bams.smk"
include: "rules/04_bam_to_bed.smk"
include: "rules/05_feature_extract.smk"
include: "rules/ichor_cna.smk"

global WORKDIR
global GENOMEBIT
global PARAMDIR
global SAMPLEDIR
global ANNO_NAME
global ANNO_FILE
global shendurelab_cfDNA
global ANNO_FILE
global SIZE_SELECTION_MIN
global SIZE_SELECTION_MAX
global SIZE_SELECTION_SUB_FOLDER
global BIN_SIZE
global CHROMS
global SCRIPTS
global THREADS
global TEMP
global SEED
# still hard coded in some places
global GENOME
global RC_NORMALIZATION
global hash_it
global MIXED_MODEL_FILE_SUFF
global HGMM_FILE_SUFF
global FFT_BINSIZE
global FFT_FREQ_FROM
global FFT_FREQ_TO
global CHROM_SIZES
global WPS_WINDOW
global WPS_MIN
global WPS_MAX
global WPS_BIGWIG_NORMALIZED
global WPS_BIN_SIZE
global WPS_BIN_NORMALIZATION
global WPS_OUT_SUFFIX

global CONTROL_SAMPLES
global WORKDIR
global hash_it

global samples
global control_samples
global control_RefName
global treated_samples
global treated_samples_names
global treated_samples_RefName

def print_config():
    print('Running cfDNA_pipeline with configuration parameters:')
    print(f'WORKDIR={WORKDIR}')
    print(f'PARAMDIR={PARAMDIR}')
    print(f'GENOMEBIT={GENOMEBIT}')
    print(f'PARAMDIR={PARAMDIR}')
    print(f'SAMPLEDIR={SAMPLEDIR}')
    print(f'ANNO_NAME={ANNO_NAME}')
    print(f'ANNO_FILE={ANNO_FILE}')
    print(f'shendurelab_cfDNA={shendurelab_cfDNA}')
    print(f'ANNO_FILE={ANNO_FILE}')
    print(f'SIZE_SELECTION_MIN={SIZE_SELECTION_MIN}')
    print(f'SIZE_SELECTION_MAX={SIZE_SELECTION_MAX}')
    print(f'SIZE_SELECTION_SUB_FOLDER={SIZE_SELECTION_SUB_FOLDER}')
    print(f'BIN_SIZE={BIN_SIZE}')
    print(f'SCRIPTS={SCRIPTS}')
    print(f'THREADS={THREADS}')
    print(f'TEMP={TEMP}')
    print(f'SEED={SEED}')
    print('chromosomes: ' + ','.join(CHROMS))
    print(f'CONTROL_SAMPLES:{CONTROL_SAMPLES}')
    print(f'control_RefName:{control_RefName}')
    print(f'treated_samples_RefName:{treated_samples_RefName}')
    print('Constants:')
    print(f'GENOME={GENOME}')
    print(f'CHROM_SIZES={CHROM_SIZES}')
    print(f'RC_NORMALIZATION={RC_NORMALIZATION}')
    print(f'hash_it={hash_it}')
    print(f'MIXED_MODEL_FILE_SUFF={MIXED_MODEL_FILE_SUFF}')
    print(f'HGMM_FILE_SUFF={HGMM_FILE_SUFF}')
    print(f'WPS_WINDOW={WPS_WINDOW}')
    print(f'WPS_MIN={WPS_MIN}')
    print(f'WPS_MAX={WPS_MAX}')
    print(f'WPS_BIGWIG_NORMALIZED={WPS_BIGWIG_NORMALIZED}')
    print(f'WPS_BIN_SIZE={WPS_BIN_SIZE}')
    print(f'WPS_BIN_NORMALIZATION={WPS_BIN_NORMALIZATION}')
    print(f'WPS_OUT_SUFFIX={WPS_OUT_SUFFIX}')
    print(f'FFT_BINSIZE={FFT_BINSIZE}')
    print(f'FFT_FREQ_FROM={FFT_FREQ_FROM}')
    print(f'FFT_FREQ_TO={FFT_FREQ_TO}')
    # print(f"correct_bias={config['correct_bias']}")

print_config()

rule do_align:
    input:
        expand(WORKDIR + "/BAM/{sample}.mem_sorted.bam",sample=samples["sample_name"]),
        expand(WORKDIR + "/BAM/{sample}.mem_sorted.bam.bai",sample=samples["sample_name"]),
        expand(WORKDIR + "/QC/" + PARAMDIR + "/{sample}-trimmed-pair1_fastqc.html",sample=samples['sample_name']),
        expand(WORKDIR + "/QC/" + PARAMDIR + "/{sample}-trimmed-pair2_fastqc.html",sample=samples['sample_name'])

rule do_bam:
    input:
        expand( WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bam",sample=samples["sample_name"])

rule do_bai:
    input:
        expand( WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bai",sample=samples["sample_name"])

rule do_qc:
    input:
        WORKDIR+ "/QC/" + PARAMDIR + "/multiqc_data"

rule do_preprocess:
    input:
        # WORKDIR+ "/QC/" + PARAMDIR + "/multiqc_data",
        expand(WORKDIR + "/BED/" + PARAMDIR + "/{sample}.bed",sample=samples["sample_name"])

# rule do_qc:
#     input:
#         WORKDIR + "/QC/trimmed_reads/" + PARAMDIR + "/multiqc_report.html",
#         WORKDIR + "/QC/Picard/" + PARAMDIR + "/multiqc_report.html"

# tMAD

rule do_cal_blacklist:  # tMAD 1
    input:
        WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"control_blacklist_"+BIN_SIZE+"Kbp_norm_TRUE.RDS"

rule do_cal_RefSample:  # tMAD 2
    input:
        WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+ control_RefName +"_norm"+RC_NORMALIZATION+"_"+BIN_SIZE+"kbp_"+GENOME+"_count.csv"


rule do_cal_t_MAD_forall:  # tMAD 3
    input:
        WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/tMAD_results.tsv"

rule do_visualising_t_MAD_forall: # tMAD 4 : last step for the visualisation
    input:
        expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/CNV_sample_{control_sample}.pdf", control_sample = control_samples[0].values),
        expand(WORKDIR +"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"tMAD/CNV_sample_{treated_sample}.pdf", treated_sample = treated_samples_names)

# ichorCNA

rule do_createPoN:
    input:
        WORKDIR + "/feature/" + PARAMDIR + "/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"pon/pon_median.rds"

rule do_ichorCNA:
    input:
        expand(WORKDIR + "/feature/" + PARAMDIR + "/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}/{sample}.cna.seg", sample=samples["sample_name"])

rule do_ichorCNA_results:
    input:
        WORKDIR + "/feature/" + PARAMDIR + "/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"ichorCNA_results.tsv"

# fragmentomics

rule do_length_binning:
    input:
        expand(WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_binned.tsv", sample=samples["sample_name"])

rule do_length_hist:
    input:
        expand(WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_hist.csv", sample=samples["sample_name"])

rule do_counts_sl_ratio:
    input:
        expand(WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_counts.tsv", sample=samples["sample_name"]),
        expand(WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_SLratio.tsv", sample=samples["sample_name"])
