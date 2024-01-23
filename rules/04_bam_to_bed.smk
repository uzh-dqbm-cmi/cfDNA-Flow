include: "common.smk"
include: "03_clean_bams.smk"

global WORKDIR
global GENOMEBIT
global PARAMDIR
global SAMPLEDIR
global CHROM_SIZES
global REFGENOME
global samples
global TEMP
global THREADS
global SCRIPTS

rule removeDuplicates:
    input:
        bam = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.markdup.bam",
        bai = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.markdup.bai",
        qc = WORKDIR+ "/QC/" + PARAMDIR + "/multiqc_data"   # require it here because if done later the input files could be deleted
    output:
        bam = temp(WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bam"),
        metrics = temp(WORKDIR + "/QC/Picard/"+PARAMDIR+"/{sample}.remdup.metrics")
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.removeDuplicates.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        tmpdir = WORKDIR + "/BAM/"+PARAMDIR,
        ram = config['MarkDuplicates']['ram']
    shell:
        "gatk --java-options '-Xmx{params.ram}G' MarkDuplicates "
        "--TMP_DIR {params.tmpdir} "
        "--CREATE_INDEX true "
        "-O {output.bam} -I {input.bam} -M {output.metrics} "
        "--REMOVE_DUPLICATES true"

rule baiRemdup:
    input:
        WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bam"
    output:
        temp(WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bai")
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "samtools index {input} {output}"

rule ComputeGCBias:
    input:
        bam = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bam",
        bai = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bai",
        genome = GENOMEBIT
    output:
        temp(WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.GC.freq.txt")
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.ComputeGCBias.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        size = config['effective_genome_size'],
        blist = config['blacklist'],
        temp = TEMP
    shell:
        """
        TMPDIR="{params.temp}"
        export TMPDIR        
        computeGCBias -b {input.bam} --effectiveGenomeSize {params.size} -g {input.genome}  --GCbiasFrequenciesFile {output} --numberOfProcessors max/2 --blackListFileName {params.blist}
        """

rule CorrectGCBias:
    input:
        bam = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bam",
        bai = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bai",
        genome = GENOMEBIT,
        freq = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.GC.freq.txt"
    output:
        temp(WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.corrected.bam")
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.CorrectGCBias.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        size = config['effective_genome_size'],
        temp = TEMP
    shell: """
        TMPDIR="{params.temp}"
        export TMPDIR            
        correctGCBias -b {input.bam} --effectiveGenomeSize {params.size} -g {input.genome} --GCbiasFrequenciesFile {input.freq} -o {output}
    """

global BAM_SPLITS
global SPLIT_SUFF

BAM_SPLITS = []
SPLIT_SUFF = ""
try:
    if config['GCBias_correction_split']:
        frg_min = 1
        frg_max = 1000000000 # max
        for frg_len in config['GCBias_correction_split'].split(','):
            BAM_SPLITS.append(f'{frg_min}:{frg_len}')
            frg_min = int(frg_len)+1
        BAM_SPLITS.append(f'{frg_min}:{frg_max}')
        print('BAM_SPLITS by fragment sizes ' + ' '.join(BAM_SPLITS))
        SPLIT_SUFF=".split" + "-".join(config['GCBias_correction_split'].split(','))
    else:
        SPLIT_SUFF = "noSplit" # to prevent activating rule bam_split_GCBiasCorrect_merge
        print('No BAM_SPLITS ')
except:
    print('No BAM_SPLITS ')

rule bam_split_GCBiasCorrect_merge:
    input:
        bam = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bam",
        bai = WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.remdup.bai"
    output:
        temp(WORKDIR + "/BAM/" + PARAMDIR + "/{sample}" + SPLIT_SUFF + ".corrected.bam")
        # temp(WORKDIR + "/BAM/"+PARAMDIR+"/{sample}.corrected.bam")
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.split.correctGCBias.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        bam_prefix=WORKDIR + "/BAM/" + PARAMDIR + "/{sample}.split",
        genome = GENOMEBIT,
        genome_size = config['effective_genome_size'],
        blist= config['blacklist'],
        temp = TEMP
    threads: THREADS
    shell:"""
        echo "bam_split_GCBiasCorrect_merge: {input.bam}\n\t BAM_SPLITS: {BAM_SPLITS}" >> bam_split_GCBiasCorrect_merge.log
        TMPDIR="{params.temp}"
        export TMPDIR       

        BAMFILES=()
        for frgs in {BAM_SPLITS}
        do 
            BAMSPLIT={params.bam_prefix}${{frgs//:/-}}.bam
            BAMCORRECT={params.bam_prefix}${{frgs//:/-}}.corrected.bam
            echo "{SCRIPTS}/bam_split.sh -i {input.bam} -r $frgs -o $BAMSPLIT" >> bam_split_GCBiasCorrect_merge.log
            {SCRIPTS}/bam_split.sh -i {input.bam} -r $frgs -o $BAMSPLIT
            
            GCFREQ=${{BAMSPLIT}}.GC.freq.txt

            echo "computeGCBias -GCbiasFrequenciesFile $GCFREQ --blackListFileName {params.blist}" >> bam_split_GCBiasCorrect_merge.log
            # echo "computeGCBias -b "$BAMSPLIT" --effectiveGenomeSize {params.genome_size} -g {params.genome}  --GCbiasFrequenciesFile $GCFREQ --numberOfProcessors max/2 --blackListFileName {params.blist}" >> bam_split_GCBiasCorrect_merge.log
            computeGCBias -b $BAMSPLIT --effectiveGenomeSize {params.genome_size} -g {params.genome}  --GCbiasFrequenciesFile $GCFREQ --numberOfProcessors max/2 --blackListFileName {params.blist}
            
            # echo "correctGCBias -b $BAMSPLIT --effectiveGenomeSize {params.genome_size} -g {params.genome} --GCbiasFrequenciesFile $GCFREQ -o $BAMCORRECT" >> bam_split_GCBiasCorrect_merge.log
            correctGCBias -b $BAMSPLIT --effectiveGenomeSize {params.genome_size} -g {params.genome} --GCbiasFrequenciesFile $GCFREQ -o $BAMCORRECT

            BAMFILES+=( "$BAMCORRECT" )
            rm "$BAMSPLIT"
            rm "$GCFREQ"
        done
        
        echo "samtools merge -@ {THREADS} -o {output} $BAMFILES[@]" >> bam_split_GCBiasCorrect_merge.log
        samtools merge -@ {THREADS} -o {output} ${{BAMFILES[@]}}
        rm ${{BAMFILES[@]}}
        """

def GCbias(wildcards):
    # print(f"in GCbias: {WORKDIR}, {PARAMDIR}, {wildcards.sample}, config['correct_bias']")
    if config['correct_bias']:
        if config['GCBias_correction_split']:
            return expand(WORKDIR + "/BAM/" + PARAMDIR + "/{s}"+SPLIT_SUFF+".corrected.bam",s=wildcards.sample)
        else:
            return expand(WORKDIR+"/BAM/"+PARAMDIR+"/{s}.corrected.bam", s=wildcards.sample)
    else:
        return expand(WORKDIR+"/BAM/"+PARAMDIR+"/{s}.remdup.bam", s=wildcards.sample)

def getFinaleDB(wildcards):
    try:
        return SAMPLEDIR+"/"+samples.loc[wildcards.sample]['finale_db']
    except:
        return

try:
    REPLACE_INCOMPLETE_NUCLEOTIDES = config['finale_db']['replace_incomplete_nucleotides']
except:
    REPLACE_INCOMPLETE_NUCLEOTIDES = False

rule finale_db_to_bam:
    input:
        getFinaleDB
    output:
        temp(WORKDIR + "/BAM/" + PARAMDIR + "/{sample}.finaleDB.bam")
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.finale_db_to_bam.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    params:
        genome=REFGENOME,
        chrom_sizes=CHROM_SIZES,
        n_opt="-N" if REPLACE_INCOMPLETE_NUCLEOTIDES else "",
        threads=THREADS,
        temp_opt=f"--temp {TEMP}" if TEMP else ""
    shell:"""
    {SCRIPTS}/fragmentstein.sh -i {input} -o {output} -t {params.threads} -g {params.genome} -c {params.chrom_sizes} {params.n_opt} {params.temp_opt}
    """

def getBamFile(wildcards):
    # print(f"in GCbias: input_data_format={config['input_data_format']}")
    if config['input_data_format'] == 'bam':
        return f"{SAMPLEDIR}/{samples.loc[wildcards.sample]['bam']}"
    if config['input_data_format'] == 'finale_db':
        return f"{WORKDIR}/BAM/{PARAMDIR}/{wildcards.sample}.finaleDB.bam"
        # return expand(WORKDIR + "/BAM/" + PARAMDIR + "/{s}.finaleDB.bam",s=wildcards.sample)
    else:
        return GCbias(wildcards)

# to rename legacy named file ca use command:
# by command: "rename final sortByName *" (LeoMed) or "rename 's/final/sortByName/' *" (ScienceCloud), depending on Linux
rule sortByName:
    input:
        bam = getBamFile
    output:
        WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByName.bam"
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.sortByName.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    shell: # can do the same check as in sortByCoord_check but for "SO:queryname"
        "samtools sort -@ 4 -n {input.bam} -o {output}"

# TODO in sortByCoord: find a better way to decide when to use getBamFile() or rules.sortByName.output
rule sortByCoord:  # assuming when input_data_format is BAM it is sorted already by coordinate
    input: getBamFile if config['input_data_format'] == 'bam' else WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByName.bam"
    output: WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bam"
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.sortByCoord.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    shell:"""
        samtools sort -@ 4 {input} -o {output}
    """

# rule sortByCoord_check:
#     input:
#         bam = getBamFile
#     output:
#         WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bam"
#     shell:
#         """
#         SO=$(samtools view -H {input.bam} | grep "SO:" | sed 's/[\t ][\t ]*/ /g'| cut -d' ' -f3)
#         if ["$S0" == "SO:coordinate"]; then
#             mv {input.bam} {output}
#         else
#             samtools sort -@ 4 {input.bam} -o {output}
#         fi
#         """

rule bamToBed:
    input:
        bam = WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByName.bam"
    output:
        WORKDIR+"/BED/"+PARAMDIR+"/{sample}.bed"
    benchmark:
        WORKDIR+"/benchmarks/"+PARAMDIR+"/{sample}.bamToBed.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "bedtools bamtobed -bedpe -i {input.bam} > {output}"