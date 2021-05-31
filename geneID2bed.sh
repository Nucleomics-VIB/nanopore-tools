#!/bin/bash
# script name: geneID2bed.sh
# extract gene regions from a EnsEMBL reference GTF file and a list of EnsEMBL gene IDs
# if -e added, extract exons regions instead
# write a BED file with results
#
# Stephane Plaisance (VIB-NC) 2021/05/31; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2021_05_31"

usage='# Usage: geneID2bed.sh
# -r <reference fasta>
# -g <gtf file>
# -t <targets (text ensgid list, one per line)>
# -f <features to extract (gene|exons)>
# script version '${version}'
# [-h for this help]'

while getopts "g:t:f:h" opt; do
  case $opt in
    g) optg=${OPTARG} ;;
    t) optt=${OPTARG} ;;
    f) optf=${OPTARG};;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check executables present
declare -a arr=( "gawk")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

if [ -z "${optg}" ]
then
  echo "# provide a GTF file matching the reference file"
  echo "${usage}" >&2
  exit 0
else
  GTF=${optg}
fi

if [ -z "${optt}" ]
then
  echo "# provide a list of targets (ensGIDs, one per line)"
  echo "${usage}" >&2
  exit 0
else
  TARGETS=${optt}
fi

if [ -z "${optf}" ]
then
  echo "# provide a feature type to extract. Choices are: 'gene' or 'exon'"
  echo "${usage}" >&2
  exit 0
else
  if [[ "${optf}" =~ ^(gene|exon)$ ]]; then
    FEATURES=${optf}
  else
    echo "only 'gene' or 'exon' are valid '-f' options"
    echo "${usage}" >&2
    exit 0
  fi
fi

# temp output folder
if [ -d "BED_tmp" ]; then 
rm -rf BED_tmp
fi

mkdir -p BED_tmp

while read ensg; do
grep ${ensg} ${GTF} | \
awk -v ensg="${ensg}" -v features="${FEATURES}" 'BEGIN{FS="\t"; OFS="\t"}
  {if ($3 == features) {
    print $1,$4,$5,ensg,$5-$4,$7
    }
  }' > BED_tmp/$(basename $ensg)_${FEATURES}.bed
done < ${TARGETS}

# merge all
cat BED_tmp/ENSG*.bed | \
sort -k 1V,1 -k 2n,2 -k 3n,3 | \
uniq > $(basename ${TARGETS} .txt)_${FEATURES}-results.bed

# print total extracted width
TOT_WIDTH=$(gawk 'BEGIN{FS="\t"; OFS="\t";tot=0}{tot=tot+$3-$2}END{print tot}' \
  $(basename ${TARGETS} .txt)_${FEATURES}-results.bed)
echo "# total ${FEATURES} reference width in $(basename ${TARGETS} .txt)_${FEATURES}-results.bed is $TOT_WIDTH bps" && \
rm -rf BED_tmp
