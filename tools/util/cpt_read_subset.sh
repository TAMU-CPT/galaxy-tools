#!/bin/bash

FILE=$1
GENOME_LENGTH=$2
DESIRED_COVERAGE=$3

# Assume approximately same # of bases per read
NUMBER_OF_BASES_IN_FASTQ=$(cat ${FILE} | awk '(NR%4==2){print $0}' | wc | awk '{print $3-$1}')
NUMBER_OF_READS_IN_FASTQ=$(echo $(wc -l ${FILE} | sed 's/ .*//g')/4 | bc)
CURRENT_COVERAGE=$(echo "${NUMBER_OF_BASES_IN_FASTQ}/${GENOME_LENGTH}" | bc)

if [ "${CURRENT_COVERAGE}" -lt "${DESIRED_COVERAGE}" ];
then
    >&2 echo "Warning: desired coverage below actual coverage (${CURRENT_COVERAGE} < ${DESIRED_COVERAGE})"
    cat ${FILE}
    exit 0
else
    AVG_READ_LENGTH=$(echo ${NUMBER_OF_BASES_IN_FASTQ}/${NUMBER_OF_READS_IN_FASTQ} | bc)
    # Must be -l as usually <1
    REDUCTION_FACTOR=$(echo ${DESIRED_COVERAGE}/${CURRENT_COVERAGE} | bc -l)
    PRETTY_PERCENT=$(echo ${REDUCTION_FACTOR} | awk '{printf("%0.02f\n", 100*$1)}')
    NUMBER_OF_READS=$(echo ${REDUCTION_FACTOR}*${NUMBER_OF_READS_IN_FASTQ}/1 | bc)
    NUMBER_OF_LINES=$(echo ${NUMBER_OF_READS}*4 | bc)
    >&2 echo "$PRETTY_PERCENT% of FastQ file used ($NUMBER_OF_READS reads)"
    head -n ${NUMBER_OF_LINES} $FILE
    exit 0
fi


