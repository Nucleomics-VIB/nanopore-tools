#!/bin/bash

# script: ONT_merge_sequencing_summaries.sh
# merge all sequencing_summary_*.txt files in path into a single text file
# to be used for QC
#
# Stephane Plaisance (VIB-NC) 2019/03/11; v1.0

path=${1}
prefix=${2:-"merged-"}

read -d '' usage <<- EOF
Usage: ONT_merge_sequencing_summaries.sh <path> <prefix>
# path to the sequencing_summary_*.txt files
# prefix for the merged file (or 'merged-' as default)
EOF

# test minimal argument
if [ -z "${path}" ]; then
   echo "# no input path provided!"
   echo -e "${usage}" >&2
   exit 1
fi

( head -1 ${path}/sequencing_summary_1.txt; \
  find ${path} -name "sequencing_summary_*.txt" -exec sed -e '1d' {} \; ) \
  > ${prefix}sequencing_summary.txt &&
  ( mkdir -p RawData; bzip2 -z -c ${prefix}sequencing_summary.txt > Rawdata/data.bz2 )
