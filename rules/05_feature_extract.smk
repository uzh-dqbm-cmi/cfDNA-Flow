include: "common.smk"
include: "04_bam_to_bed.smk"
include: "size_selection.smk"
include: "feature/anno.snake"
include: "feature/alignments.snake"
include: "copy_number_aberration.smk"

global WORKDIR
global PARAMDIR
global ANNO_FILE
global ANNO_NAME
global CHROMS
global CHROM_SIZES
global THREADS
global shendurelab_cfDNA
global SIZE_SELECTION_SUB_FOLDER
global HGMM_FILE_SUFF
global MIXED_MODEL_FILE_SUFF
global SCRIPTS
global samples

def add_chr(wildcards):
    if '19' in config['refGenome'] or '37' in config['refGenome']:
        path = config["outputDir"]+"/BED/"+config["paramDir"]+'/'
        return "for file in "+ path +"*.bed; do awk 'BEGIN{OFS=\"\\t\"}$1=\"chr\"$1' \"$file\" | awk 'BEGIN{OFS=\"\\t\"}$4=\"chr\"$4' > "+path+"tmp; cat "+path+"tmp > \"$file\"; rm "+path+"tmp; done;"
    else:
        return ''

rule contig_summary:
    input:
        bams=expand(WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bam", sample = samples['sample_name'].values),
        bais=expand(WORKDIR+"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bai", sample = samples['sample_name'].values)
    output:
        WORKDIR + "/feature/" + PARAMDIR + "/contig_summary.tsv",
    params:
        output = WORKDIR + "/feature/" + PARAMDIR + "/contig_summary.tsv",
        virus_contig = "NC_001526.4"
    shell:"""
        echo -e "file\\tall\\tX\\tY\\tmito\\tvirus" > "{params.output}"
        for bam in {input.bams}; do
            echo -n "$bam" >> "{params.output}"
            samtools idxstats "$bam" | awk -F '\\t' '$1 != "{params.virus_contig}" {{all_sum += $3}}; $1 == "{params.virus_contig}" {{virus = $3}}; $1 == "chrM" {{mito = $3}}; $1 == "chrX" {{X = $3}}; $1 == "chrY" {{Y = $3}} END {{print "\\t"all_sum"\\t"X"\\t"Y"\\t"mito"\\t"virus}}' >> "{params.output}"
        done
        """

rule length_binning:
    input:
        WORKDIR +"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bam"
    output:
        WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_binned.tsv",
    wildcard_constraints:
        sample="[^/]+"
    params:
        chrom_sizes=CHROM_SIZES,
        contigs=' '.join(CHROMS),
        binsize=100000,
    shell:"""
        contigs=("{params.contigs}")
        binsize={params.binsize}
        ## : > {output}
        while IFS=$'\\t' read -r chrom length
        do
            if [[ "$contigs" =~ (" "|^)$chrom(" "|$) && $((length/binsize)) > 0 ]]; then
                i=1
                while (( i <= $((length/binsize)) ))
                do
                    echo -e -n $chrom"\\t"$(( ((i-1)*binsize)+1 ))"\\t"$(( i*binsize))"\\t" >> "{output}"
                    (( i=$i+1 ))
                    samtools view "{input}" $chrom\\:$(( ((i-1)*binsize)+1 ))\\-$(( i*binsize)) | cut -f 9 | sort -n | uniq -c | awk '$2 > 99 && $2 <151{{short += $1}} $2 > 150 && $2 <221{{long += $1}} $2>0 && $2 <1000{{nfrags += $1}} $2>0 && $2 <1000{{printf $2":"$1", "}} END {{print "\\t"short"\\t"long"\\t"nfrags}}' >> "{output}"
                done
            fi
        done < {params.chrom_sizes}
        """

rule length_hist:
    input:
        WORKDIR +"/BAM/"+PARAMDIR+"/{sample}.sortByCoord.bam"
    output:
        WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_hist.csv"
    wildcard_constraints:
        sample="[^/]+"
    shell:"""
        samtools view {input} | grep -E "chr[1-9]" | cut -f 9 | sort | uniq -c | awk '$2>0{{print $2,$1}}' | sort -n > {output}        
        """

rule counts_sl_ratio:
    input:
        WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_binned.tsv",
    output:
        counts=WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_counts.tsv",
        sl_ratio=WORKDIR + "/feature/" + PARAMDIR + "/length/{sample}_SLratio.tsv"
    wildcard_constraints:
        sample="[^/]+"
    shell:"""
        cut -f 1,2,3,7 "{input}" >  "{output.counts}"
        awk -F"\t" '$6{{print $1, $2, $3, $5/$6; next}} {{print $1, $2, $3, 0}}' "{input}" > "{output.sl_ratio}"
      """

rule bedProcess:
    input:
        expand(WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}.bed", sample = samples['sample_name'].values)
    output:
        expand(WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}_len.csv",sample=samples['sample_name'].values),
        expand(WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}_point.csv",sample=samples['sample_name'].values),
        WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"mergeddf.csv"
    params:
        bed_dir = WORKDIR+"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER,
        beds = ' '.join(samples['sample_name'].values),
        addChr = add_chr
    shell:
        "{params.addChr} "
        "python {SCRIPTS}/fasterbedprocessing.py "
        "-d {params.bed_dir} -i {params.beds}"


rule globalLength:
    input:
        expand(WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}_len.csv",sample=samples['sample_name'].values)
    output: # TODO check if the output folder works so
        WORKDIR + "/feature/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"global_length.tsv"
    shell:
        "python {SCRIPTS}/global_length.py -i {input} -o {output}"


rule len_uniq_count:
    input:
        WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}_len.csv"
    output:
        WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"{sample}_lenuniqcount.csv"
    wildcard_constraints:
        sample="[^/]+"
    shell:"""
        sort -g {input} | uniq -c | awk 'BEGIN {{OFS="\t"}}; {{print $2,$1}}' > {output}
        """


rule uniq_counts:
    input:
        expand(WORKDIR+"/BED/"+PARAMDIR+"/"+SIZE_SELECTION_SUB_FOLDER+"{sample}_lenuniqcount.csv", sample=samples['sample_name'].values),
    output:
        WORKDIR + "/BED/" + PARAMDIR + "/"+SIZE_SELECTION_SUB_FOLDER+"uniqcounts.tsv"
    wildcard_constraints:
        sample="[^/]+"
    params:
        ext = "_lenuniqcount.csv",
        lb = 76,
        ub = 1000
    run:
        import pandas as pd
        import os

        master = pd.DataFrame(list(range(params.lb, params.ub)), columns= ["length"])
        for file in input:
            sample = os.path.basename(file).replace(params.ext, "")
            print(f'merge {sample}: {file}')
            csv = pd.read_table(file, sep="\t", names=["length", sample])
            master = pd.merge(master, csv, how="left", on=["length"])
        print(f'save in {output}')
        with pd.option_context('display.max_rows', 10,
                'display.max_columns', None,
                'display.precision', 3):
            print(master)
        master.to_csv(f'{output}', sep="\t", index=False)
