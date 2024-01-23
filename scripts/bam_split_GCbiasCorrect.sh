#!/bin/bash
usage()
{
    echo -e "\nDo size selection for a given fragment range\n"
    echo -e "Usage: $0 [options] -i input_file.bed --min 90 --max 150 -o output_file.bed\n"
    echo -e "required arguments:\n\t-i FILEPATH\t\t .bed file (BEDPE format)\n\t-o FILEPATH\t\t output .bed file\n"
    echo -e "\t-i, --input\t\tinput BAM file\n"
    echo -e "\t-o, --output\t\toutput corrected BAM file\n"
    echo -e "\t--bam_prefix\t\t bam_prefix\n"
    echo -e "\t--bam_splits\t\t BAM SPLITS\n"
    echo -e "\t--genome\t\t genome 2bit file\n"
    echo -e "\t--genome_size\t\t genome size file\n"
    echo -e "\t--blist\t\t black list\n"
	  echo -e "\t--temp\t\t\ttemp folder (default: the folder of the input bed file)\n"
    echo -e "\t-h, --help\t\tshow this help message and exit\n"
}

while [[ "$1" != "" ]]
do
    PARAM=$(echo $1)
    VALUE=$(echo $2)
    case ${PARAM} in
        -h|--help)
            usage
            exit
            ;;
        -i|--input)
            INPUT="${VALUE}"
            ;;
        -o|--output)
            OUTPUT="${VALUE}"
            ;;
        -l|--lines)
          START_LINE=$(echo ${VALUE} | cut -d: -f1)
          END_LINE=$(echo ${VALUE} | cut -d: -f2)
          ;;
        --min)
          MIN=${VALUE}
          ;;
        --max)
          MAX=${VALUE}
          ;;
        -c|--chunksize)
          CHUNKSIZE=${VALUE}
          ;;
        -p)
          CORES=${VALUE}
          ;;
        --temp)
          TEMP_FOLDER=${VALUE}
          ;;
            *)
                echo -e "\nERROR: unknown parameter \"${PARAM}\"\n"
                usage
                exit 1
                ;;
        esac
        shift 2
done

bam_split()
{
  local min="$1"
  local max="$2"
  if [[ "$min" = "" ]]
  then
    min=1
  fi

  if [[ "$max" = "" ]]
  then
    max=1000000000
  fi

  if [[ "${OUTPUT}" = "" ]]
  then
    OUTPUT="${BAMFILE}.f${min}-${max}.bam"
  fi

  # echo "samtools view ${BAMFILE} -h -@ 8 | awk 'function abs(x){return ((x < 0.0) ? -x : x)} (abs($9)>='${MIN}' && abs($9)<='${MAX}') || /^@/' | samtools view -S -b -o ${OUTPUT} -@ 8 -" >> bam_split_GCBiasCorrect_merge.log
  samtools view ${BAMFILE} -h -@ 8 | awk 'function abs(x){return ((x < 0.0) ? -x : x)} (abs($9)>='${min}' && abs($9)<='${max}') || /^@/' | samtools view -S -b -o ${OUTPUT} -@ 8 -
  # | samtools view -S -b -o ${OUTPUT} -
  # > ${OUTPUT}.tsv
  # || '/^(\#)/'

  samtools index ${OUTPUT}
}

if [[ -z "${BED_FILE}" ]]
then
	echo -e "\nERROR: missing parameter(s)"
	usage
	exit 1
fi

if [[ -z "${TEMP_FOLDER}" ]]
then
  TEMP_FILE="${BED_FILE}_TEMP"
else
  [ ! -d "$TEMP_FOLDER" ] && mkdir -p "${TEMP_FOLDER}"
  filename="$(basename -- ${BED_FILE})"
  TEMP_FILE="${TEMP_FOLDER}/${filename}_TEMP"
fi

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
