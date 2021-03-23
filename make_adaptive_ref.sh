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
# -s <additional edge length (default to 5000)>
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

# check executables present
declare -a arr=( "samtools" "bedtools")
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
SLOPPED_BED=${BED%.*}_slop-${BASES_TO_EXPAND_PER_SIDE}.bed

# This is the final file which you will upload into MinKNOW:
SUBSETTED_FASTA=${REF%.*}-${SLOPPED_BED%.*}.fasta

# index fasta
samtools faidx ${REF}

# create chr-size file for bedtools
cut -f1,2 ${REF}.fai > ${CHROM_SIZES}

# sort input BED by chr then start
sort -k 1V,1 -k 2n,2 ${BED} \
  > sorted_${BED}

# create expanded BED
bedtools slop -l ${BASES_TO_EXPAND_PER_SIDE} \
  -r ${BASES_TO_EXPAND_PER_SIDE} \
  -i sorted_${BED} \
  -g ${CHROM_SIZES} \
  > ${SLOPPED_BED}_ini

# merge region overlaps where present
bedtools merge -i ${SLOPPED_BED}_ini \
  > ${SLOPPED_BED}

# extract fasta sequences
bedtools getfasta -fi ${REF} \
  -bed ${BED} \
  -fo ${SUBSETTED_FASTA} \
  -name

echo "# the file ${SUBSETTED_FASTA} can be used in Minknow for adaptive sequencing"

exit 0

# ref: https://community.nanoporetech.com/info_sheets/adaptive-sampling/v/ads_s1016_v1_revb_12nov2020
# Navigate to the folder containing your FASTA and .bed files. If they are in different folders, create a link to the location of the FASTA file (if you are
# already in the folder containing the .bed file):
# ln -s /long/path/to/my/reference/in/different/folder/myref.fasta myref.fasta
# REF=myref.fasta
# or the .bed file (if you are already in the folder containing the FASTA file):
# ln -s /long/path/to/my/bed/in/different/folder/mybed.bed mybed.bed
# BED=mybed.bed
# 4. Edit the settings for your reference and your .bed with target regions:
# BASES_TO_EXPAND_PER_SIDE=half of your expected N50. For example, set 5000, i.e. 5000 bp either side for a N50 = 10 kb run.
# REF=your_reference.fasta
# BED=your_bed_with_target_regions.bed
# (The subsetted FASTA file will be called your_reference-your_bed_with_target_regions.fasta)
# These are the intermediate files that will be created:
# CHROM_SIZES=${REF}.chrom.sizes
# SLOPPED_BED=${BED%.*}_slop-${BASES_TO_EXPAND_PER_SIDE}.bed
# These will be saved in the same location where the commands are being run, and will remain until they are deleted.
# This is the final file which you will upload into MinKNOW:
# SUBSETTED_FASTA=${REF%.}-${SLOPPED_BED%.}.fasta
# This is also saved in the same location where the commands are being run.
# 5. Index the reference and get chromosome sizes:
# samtools faidx ${REF}
# cut -f1,2 ${REF}.fai > ${CHROM_SIZES}
# 6. Expand the .bed and extract the FASTA from the expanded .bed:
# bedtools slop -l ${BASES_TO_EXPAND_PER_SIDE} -r ${BASES_TO_EXPAND_PER_SIDE} -i ${BED} -g ${CHROM_SIZES} > ${SLOPPED_BED}
# bedtools getfasta -fi ${REF} -bed ${BED} -fo ${SUBSETTED_FASTA} -name
# These commands take the .bed file that contains the regions of interest and add a number of bases on either side of the ROI. This is to ensure that
# during adaptive sampling, reads which are on the boundary of the ROIs have the chance of overlapping the ROI when they are fully read.
# 7. Copy the subsetted FASTA file to your MinION Mk1C.