#!/bin/bash

# script: ONT_merge_fastq.sh
# merge all *.fastq files (pass and fail!) in path into a single compressed file
#
# Stephane Plaisance (VIB-NC) 2019/03/20; v1.0

path=${1}
prefix=${2:-"$(basename ${path})"}

read -d '' usage <<- EOF
Usage: ONT_merge_fastq.sh <path> <prefix>
# path containing the fastq files (will merge pass and fail files!!)
# prefix for the merged file (or 'path' as default)
EOF

# test minimal argument
if [ -z "${path}" ]; then
   echo "# no input path provided!"
   echo -e "${usage}" >&2
   exit 1
fi

# merge all fastq files in path and compress
cat /dev/null > ${prefix}_reads.fq.gz
find ${path} -name "*.fastq" -type f -exec sh -c "cat {} | \
  bgzip -c >> ${prefix}_reads.fq.gz" \;
