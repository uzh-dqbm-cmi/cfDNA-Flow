#!/bin/bash

usage()
{
  echo -e "\nCorrect GC bias for a given BAM file\n"
  echo -e "Usage: $0 [options] -i input_bam_file ...\n"
	echo -e "required arguments:\n"
	echo -e "\t-i FILEPATH\t\tsorted and indexed .bam file\n"
  echo -e "optional arguments:\n"
  echo -e "\t-o, --output\t\toutput bed file  (default: '-' standard output)\n"
	echo -e "\t-h, --help\t\tshow this help message and exit\n"
}

TMPDIR=TMP

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
        -g|--genome)
            GENOME="${VALUE}"
            ;;
        -s|--genome_size)
            GSIZE="${VALUE}"
            ;;
        -b|--blist)
            BLIST="${VALUE}"
            ;;
        -t|--temp)
            TMPDIR="${VALUE}"
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

if [[ "${OUTPUT}" = "" ]]
then
  OUTPUT="${BAMFILE}.correctGCbias.bam"
fi

if [[ "${GSIZE}" = "" ]]
then
  GSIZE=2805636331
fi

GCFREQ="${BAMFILE}.GC.freq.txt"

export TMPDIR

echo "computeGCBias -b ${BAMFILE} --effectiveGenomeSize ${GSIZE} -g ${GENOME}  --GCbiasFrequenciesFile ${GCFREQ} --numberOfProcessors max/2 --blackListFileName ${BLIST}"
computeGCBias -b ${BAMFILE} --effectiveGenomeSize ${GSIZE} -g ${GENOME}  --GCbiasFrequenciesFile ${GCFREQ} --numberOfProcessors max/2 --blackListFileName ${BLIST}

echo "correctGCBias -b ${BAMFILE} --effectiveGenomeSize ${GSIZE} -g ${GENOME} --GCbiasFrequenciesFile ${GCFREQ} -o ${OUTPUT}"
correctGCBias -b ${BAMFILE} --effectiveGenomeSize ${GSIZE} -g ${GENOME} --GCbiasFrequenciesFile ${GCFREQ} -o ${OUTPUT}
