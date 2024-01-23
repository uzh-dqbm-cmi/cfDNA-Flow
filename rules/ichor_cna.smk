include: "common.smk"
from typing import Any

global WORKDIR
global PARAMDIR
global ANNO_FILE
global ANNO_NAME
global THREADS
global shendurelab_cfDNA
global SIZE_SELECTION_SUB_FOLDER
global MIXED_MODEL_FILE_SUFF
global SCRIPTS
global CONTROL_SAMPLES
global samples

def get_option_config(option_prefix, config_section, config_name, default:Any=""):
    try:
        value = config[config_section][config_name]
    except:
        value = default
    if value == "" or value is None:
        return ""
    else:
        if type(value) is str:
            return f'{option_prefix} "{value}"'
        else:
            return f'{option_prefix} {value}'

# readCounter option for index creating (-b) doesn't work, hence:
# create .bam.bai from .bai files (samtools index output) so that it matches required input format
rule readCounterIndex:
    input:
        WORKDIR + "/BAM/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.sortByCoord.bai",
    output:
        WORKDIR + "/BAM/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.sortByCoord.bam.bai",
    shell:
     "cp {input} {output}"

rule readCounter:
    input:
        bais = WORKDIR + "/BAM/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.sortByCoord.bam.bai",
        bams = WORKDIR+"/BAM/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.sortByCoord.bam"
    output:
        wigs = "{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.wig"
    params:
        #chrs = ichorChrs,
        chrs=config["chromosomes"],
        binSize=config["ichorCNA"]["ichor_bin_size"],
        qual=config['ichorCNA']["ichor_qual"]
    shell:
        "readCounter {input.bams} "
        "-c {params.chrs} "
        "-w {params.binSize} "
        "-q {params.qual} > {output.wigs}"

rule createPoNWigs:
    input:
       controltxt = CONTROL_SAMPLES
    params:
        wigdir = "{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER,
        outputdir = "{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER
    output:
        "{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"pon_wigs.txt"
    shell:
        "python {SCRIPTS}/helper/pon_wigs.py "
        "--wig_dir {params.wigdir} "
        "--control_samples {input.controltxt} "
        "--output_destination {params.outputdir}"

rule createPoN:
    input:
        ponwigs = "{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"pon_wigs.txt",
        wigs = expand(WORKDIR + "/feature/" + PARAMDIR + "/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.wig",sample=samples["sample_name"]),
    output:
        pon_rds = "{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"pon/pon_median.rds",
        pon_txt = "{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"pon/pon_median.txt",
    params:
        gcwig = config["ichorCNA"]["gcwig"],
        mapwig = config["ichorCNA"]["mapwig"],
        centromere = config["ichorCNA"]["centromere"],
        basename = "{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"pon/pon"
    shell:
        "Rscript {SCRIPTS}/createPanelOfNormals.R "
        "--filelist {input.ponwigs} "
        "--gcWig {params.gcwig} "
        "--mapWig {params.mapwig} "
        "--centromere {params.centromere} "
        "--outfile {params.basename}"

def ichorPoN(wildcards):
    if config["ichorCNA"]["ichor_pon"]:
        return WORKDIR+"/feature/"+PARAMDIR+"/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"pon/pon_median.rds"
    else:
        return []

rule ichorCNA:
    input:
        wigs = rules.readCounter.output,
        pon = ichorPoN
    output:
        corrDepth="{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}/{sample}.correctedDepth.txt",
        param="{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}/{sample}.params.txt",
        cna="{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}/{sample}.cna.seg",
        segTxt="{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}/{sample}.seg.txt",
        seg="{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}/{sample}.seg",
        rdata="{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}/{sample}.RData"
    params:
        outDir="{WORKDIR}/feature/{PARAMDIR}/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}",
        id="{sample}",
        gcwig=config["ichorCNA"]["gcwig"],
        mapwig=config["ichorCNA"]["mapwig"],
        seqinfo=config["ichorCNA"]["seqinfo"],
        centromere=config["ichorCNA"]["centromere"],
        genomeBuild=config["ichorCNA"]["ichor_genome"],
        pon=ichorPoN,
        opt_normal=get_option_config("--normal","ichorCNA","normal"),
        opt_ploidy=get_option_config("--ploidy", "ichorCNA","ploidy"),
        opt_chrs=get_option_config("--chrs","ichorCNA","chrs"),
        opt_chrTrain=get_option_config("--chrTrain","ichorCNA","chrTrain"),
        opt_maxCN=get_option_config("--maxCN","ichorCNA","maxCN"),
        opt_scStates=get_option_config("--scStates","ichorCNA","scStates"),
        opt_estimateScPrevalence=get_option_config("--estimateScPrevalence","ichorCNA","estimateScPrevalence")
    shell:"""
        Rscript {SCRIPTS}/runIchorCNA.R --id '{params.id}' --WIG {input.wigs}  --gcWig {params.gcwig} \
            --mapWig {params.mapwig} --seqInfo {params.seqinfo} --centromere {params.centromere} \
            --normalPanel {input.pon} {params.opt_normal} {params.opt_ploidy} {params.opt_chrs} {params.opt_chrTrain} \
            {params.opt_maxCN} {params.opt_scStates} {params.opt_estimateScPrevalence} \
            --genomeBuild {params.genomeBuild} --outDir {params.outDir}
    """


rule ichorCNA_results:
    input:
        expand(WORKDIR + "/feature/" + PARAMDIR + "/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"{sample}/{sample}.params.txt",sample=samples["sample_name"])
    output:
        WORKDIR + "/feature/" + PARAMDIR + "/ichorCNA/"+SIZE_SELECTION_SUB_FOLDER+"ichorCNA_results.tsv"
    run:
        import do_ichorCNA_results
        do_ichorCNA_results.main(input_file=input, output=output)
