#!/bin/bash

usage()
{
  echo -e "\nSplit BAM file for a given fragment range\n"
  echo -e "Usage: $0 [options] -i input_bam_file -r region\n"
	echo -e "required arguments:\n"
	echo -e "\t-i FILEPATH\t\tsorted and indexed .bam file\n"
	echo -e "\t-r, --region\t\tmin:max fragment length (e.g. 1:150)\n\n"
  echo -e "optional arguments:\n"
  echo -e "\t-o, --output\t\toutput bed file  (default: '-' standard output)\n"
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
        -i)
            BAMFILE="${VALUE}"
            ;;
        -o|--output)
            OUTPUT="${VALUE}"
            ;;
        -r|--region)
          MIN=$(echo ${VALUE} | cut -d: -f1)
          MAX=$(echo ${VALUE} | cut -d: -f2)
          ;;
        *)
            echo -e "\nERROR: unknown parameter \"${PARAM}\"\n"
            usage
            exit 1
            ;;
    esac
    shift 2
done

if [[ "${BAMFILE}" = "" ]]
then
	echo -e "\nERROR: missing parameter(s)"
	usage
	exit 1
fi

if [[ "${MIN}" = "" ]]
then
  MIN=1
fi

if [[ "${MAX}" = "" ]]
then
  MAX=1000000000
fi

if [[ "${OUTPUT}" = "" ]]
then
  OUTPUT="${BAMFILE}.f${MIN}-${MAX}.bam"
fi

# echo "samtools view ${BAMFILE} -h -@ 8 | awk 'function abs(x){return ((x < 0.0) ? -x : x)} (abs($9)>='${MIN}' && abs($9)<='${MAX}') || /^@/' | samtools view -S -b -o ${OUTPUT} -@ 8 -" >> bam_split_GCBiasCorrect_merge.log
samtools view ${BAMFILE} -h -@ 8 | awk 'function abs(x){return ((x < 0.0) ? -x : x)} (abs($9)>='${MIN}' && abs($9)<='${MAX}') || /^@/' | samtools view -S -b -o ${OUTPUT} -@ 8 -
# | samtools view -S -b -o ${OUTPUT} -
# > ${OUTPUT}.tsv
# || '/^(\#)/'

samtools index ${OUTPUT}
