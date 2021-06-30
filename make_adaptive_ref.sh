#!/bin/bash
# script name: make_adaptive_ref.sh
# prepare reference file for adaptive sampling on ONT
#
# Stephane Plaisance (VIB-NC) 2021/03/23; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2021_03_23"

usage='# Usage: make_adaptive_ref.sh
# -r <reference fasta>
# -b <bed regions>
# -s <additional flank length (default to 5000)>
# script version '${version}'
# [-h for this help]'

while getopts "r:b:s:h" opt; do
  case $opt in
    r) optr=${OPTARG} ;;
    b) optb=${OPTARG} ;;
    s) opts=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check executable present
declare -a arr=( "samtools" "bedtools" "gawk")
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# check INPUTS
if [ -z "${optr}" ]
then
  echo "# provide a reference fasta file to be parsed"
  echo "${usage}" >&2
  exit 0
else
  REF=${optr}
fi

if [ -z "${optb}" ]
then
  echo "# provide a BED file with selected region(s)"
  echo "${usage}" >&2
  exit 0
else
  BED=${optb}
fi

# additional side length
BASES_TO_EXPAND_PER_SIDE=${opts:-5000}

# OUTPUTS
CHROM_SIZES=${REF}.chrom.sizes
SLOPPED_BED=${BED%.*}_slop_${BASES_TO_EXPAND_PER_SIDE}.bed

# This is the final file which you will upload into MinKNOW:
SUBSETTED_FASTA=$(basename ${REF%.*})_${SLOPPED_BED%.*}.fasta

# index fasta if absent
if [ ! -f "${REF}.fai" ] ; then
samtools faidx ${REF}
fi

# create chr-size file for bedtools
cut -f1,2 ${REF}.fai > ${CHROM_SIZES}

# sort input BED by chr then start
sort -k 1V,1 -k 2n,2 ${BED} \
  > sorted_${BED}

# create expanded BED
bedtools slop \
  -l ${BASES_TO_EXPAND_PER_SIDE} \
  -r ${BASES_TO_EXPAND_PER_SIDE} \
  -i sorted_${BED} \
  -g ${CHROM_SIZES} \
  > ini_${SLOPPED_BED}

# merge region overlaps where present and collapse their descriptions as a csv list
bedtools merge -i ini_${SLOPPED_BED} \
  -c 4 \
  -o collapse \
  > ${SLOPPED_BED}

# print total reference width
TOT_WIDTH=$(gawk 'BEGIN{FS="\t"; OFS="\t";tot=0}{tot=tot+$3-$2}END{print tot}' \
  ${SLOPPED_BED})
echo "# total reference width in ${SLOPPED_BED} is $TOT_WIDTH bps"

# extract fasta sequences
bedtools getfasta -fi ${REF} \
  -bed ${SLOPPED_BED} \
  -fo ${SUBSETTED_FASTA} \
  -name

echo "# the file ${SUBSETTED_FASTA} can be used in Minknow for adaptive sequencing"

exit 0

# create bed track from gff data for genes
# gawk 'BEGIN{FS="\t"; OFS="\t"}{
#   if($3=="gene") {
#     split($9,ann,";"); ensid=gensub(/gene_id \"(.*)\"/,"\\1","g",ann[1]); 
#     hugo=gensub(/ gene_name \"(.*)\"/,"\\1","g",ann[3]); 
#     print $1,$4,$5,ensid"|"hugo
#     }
#   }' *.gtf > genes.bed

# ref: https://community.nanoporetech.com/info_sheets/adaptive-sampling/v/ads_s1016_v1_revb_12nov2020