#!/bin/bash
usage()
{
    echo -e "\nDo size selection for a given fragment range\n"
    echo -e "Usage: $0 [options] -i input_file.bed --min 90 --max 150 -o output_file.bed\n"
    echo -e "required arguments:\n\t-i FILEPATH\t\t .bed file (BEDPE format)\n\t-o FILEPATH\t\t output .bed file\n"
    echo -e "\t-o, --output\t\toutput bed file  (default: '-' standard output)\n"
    echo -e "\t-l, --lines\t\tlines region to be processed (e.g. 1000000:2000000)\n"
    echo -e "\t--min\t\t\tminimum length of fragments to be included (default: 90)\n"
    echo -e "\t--max\t\t\tmaximum length of fragments to be included (default: 150)\n"
    echo -e "\t-c, --chunksize\t\t\chunk size to split the input BED file (default: 1000000)\n"
    echo -e "\t-p\t\t\tNumber of parallel threads to be used. Job will be split into chunks by the given chunksize. (default: 1)\n"
	  echo -e "\t--temp\t\t\ttemp folder (default: the folder of the input bed file)\n"
    echo -e "\t-h, --help\t\tshow this help message and exit\n"
}

# defaults
MIN=90
MAX=150
CORES=1

START_LINE=1
END_LINE=-1
CHUNKSIZE=1000000
OUTPUT="-"
TEMP_FOLDER=""
PATH_SEP="/"
TEMP_FILE=""

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
            BED_FILE="${VALUE}"
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

#scriptdir=$(dirname "$0")
this_script=$0

if [[ "${END_LINE}" -lt 0 ]]
then
  END_LINE=$(( $(wc -l < "${BED_FILE}") ))
fi

# debug print
#echo -e "-i ${BED_FILE} --lines $START_LINE:$END_LINE --chunksize $CHUNKSIZE -o $OUTPUT" --temp $TEMP_FOLDER >> "${BED_FILE}"_TEMP.log

if [[ $(( END_LINE - START_LINE )) -gt $((CHUNKSIZE)) ]]
then
  # chunk it
	chunks=$(( (END_LINE-START_LINE+1) / CHUNKSIZE))
	START_ARR=($(seq 0 ${chunks} | awk -v OFS=" " '
			{if ('${START_LINE}' + (($1 + 1)  * '${CHUNKSIZE}') < '${END_LINE}') {
				print '${START_LINE}' + ($1  * '${CHUNKSIZE}') ":" '${START_LINE}' + (($1 + 1) * '${CHUNKSIZE}' - 1)}
			else {
				print '${START_LINE}' + ($1  * '${CHUNKSIZE}')":" '${END_LINE}'}
			}'))
	# debug print
	#echo -e "${this_script} --min ${MIN} --max ${MAX} -i ${BED_FILE}  (${END_LINE}-${START_LINE}+1)/${CHUNKSIZE} chunks:${chunks} => -l ::: " "${START_ARR[@]}" >> "${BED_FILE}"_TEMP.log
	if [[ -e "${TEMP_FILE}" ]]; then
	  rm "${TEMP_FILE}"
  fi
  parallel -k -j ${CORES} --eta "${this_script}" --min ${MIN} --max ${MAX} -i "${BED_FILE}" -l {} ::: "${START_ARR[@]}" >> "${TEMP_FILE}"
	if [[ "$OUTPUT" == "-" ]]; then
  	cat "${TEMP_FILE}"; rm "${TEMP_FILE}"
  else
    mv "${TEMP_FILE}" "${OUTPUT}"
  fi
  exit
fi

if [[ "$OUTPUT" == "-" ]]; then
  sed -n "${START_LINE},${END_LINE}p;$((END_LINE+1))q" "${BED_FILE}" | awk '($6-$2)>='${MIN}' && ($6-$2)<='${MAX}''
else
  sed -n "${START_LINE},${END_LINE}p;$((END_LINE+1))q" "${BED_FILE}" | awk '($6-$2)>='${MIN}' && ($6-$2)<='${MAX}'' > "${OUTPUT}"
fi
