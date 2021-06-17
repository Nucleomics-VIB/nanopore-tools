#!/bin/bash
# script name: run_chopchop_insideout.sh
# enrich a both flanks from an inserted sequence with ONT Cas9
# run chopchop prediction within an inserted sequence (reporter gene)
# get predictions on the minus strand to enrich the left flank (5' from the insertion)
# get predictions on the plus strand to enrich the right flank (3' from the insertion)
# REM: use chopchop with a reference genome index complemented with the query sequence
#
# Stephane Plaisance (VIB-NC) 2021/05/04; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.1, 2021_05_10"

usage='# Usage: run_chopchop_insideout.sh
# -t <target inserted sequence name (eg. eGFP_vect:0-720)>
# -g <reference genome used by ChopChop (eg. hg38_orfs)>
# script version '${version}'
# [-h for this help]'

while getopts "r:g:h" opt; do
  case $opt in
    r) optr=${OPTARG} ;;
    g) optg=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# defaults
limitPrintResults=400

# check executables present
declare -a arr=( "chopchop.py" "bedtools" "tee" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# check INPUTS
if [ -z "${optr}" ]
then
  echo "# provide a target sequence region present in the reference index (genome; eg eGFP_vect:0-720)"
  echo "${usage}" >&2
  exit 0
fi

if [ -z "${optg}" ]
then
  echo "# specify the reference index (genome) to be search against and containing the target sequence (eg hg38_orfs)"
  echo "${usage}" >&2
  exit 0
fi

# add all logs to a global log file
datetag=$(date +%s)
basedir=$PWD
exec &> >(tee -a "${basedir}/run_chopchop_insideout_${datetag}.log")

# provide the target sequence name from which to capture)
# this name should be present in the adapted reference genome given by optg
target=${optr}

# define reference genome or take default hg38 (should be installed on your machine)
refg=${optg}

# default primer3 settings
primer3args='PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60'
backbone="AGGCTAGTCCGT"
scoring="DOENCH_2016"

# default
#scoring="G_20"

# create index name and folders for results
pfx="${basedir}/${optr//\:/_}_results_${datetag}"
mkdir -p ${pfx}

# search with chopchop
echo "# searching for guides"

# removed
#		--displaySeqFlanks 300 \
#		--limitPrintResults ${limitPrintResults} \
#		--makePrimers \
#		--primerFlanks 290 \
#		--jsonVisualize \

# adding --padSize 0 to avoid extension outside of the sequence
# Extra bases searched outside the exon. Defaults to the
#   size of the guide RNA for CRISPR and TALEN + maximum
#   spacer for TALEN
# https://bitbucket.org/valenlab/chopchop/issues/7/errors-creating-custom-reference-genome

# from web-run:/query.json
# ["-J"
#  "-BED"
#  "-GenBank"
#  "-G" "hg38"
#  "-filterGCmin  "20"
#  "-filterGCmax  "80"
#  "-filterSelfCompMax  "0"
##  for Fasta input "-rm1perfOff"
#  "-consensusUnion"
#  "-t" "WHOLE"
#  "-n" "N"
#  "-R" "4"
#  "-T" "1"
#  "-g" "20"
#  "-scoringMethod  "DOENCH_2016"
#  "-f" "NN"
#  "-v  "3"
#  "-M  "NGG"
#  "-BB  "AGGCTAGTCCGT"
#  "--nonOverlapping"]

chopchop.py \
		--MODE 1 \
		--PAM NGG \
		--maxMismatches 3 \
		--guideSize 20 \
		--padSize 0 \
		-scoringMethod ${scoring} \
		--fivePrimeEnd NN \
		--backbone ${backbone} \
		--target WHOLE \
		--enzymeCo N \
		--minResSiteLen 4 \
		--primer3options ${primer3args} \
		--guidePadding 20 \
		--BED \
		--GenBank \
		--offtargetsTable \
		--filterGCmin 20 \
		--filterGCmax 80 \
		--filterSelfCompMax 0 \
		--consensusUnion \
		--outputDir ${pfx} \
		--genome ${refg} \
		--targets ${target} \
		> ${pfx}/results.txt \
		2> ${pfx}/python.err

# filter best hits on - strand
gawk 'BEGIN{FS="\t"; OFS="\t"}{if( (NR==1) || ($4~/-/ && $6==0 && $11>0.3)) print $0}' \
	${pfx}/results.txt | column -s $'\t' -t \
	> ${pfx}/results-filtered_minus.txt

# filter best hits on + strand
gawk 'BEGIN{FS="\t"; OFS="\t"}{if( (NR==1) || ($4~/+/ && $6==0 && $11>0.3)) print $0}' \
	${pfx}/results.txt | column -s $'\t' -t \
	> ${pfx}/results-filtered_plus.txt

# report results
echo "# => found $(($(wc -l ${pfx}/results.txt|cut -d" " -f1)-1)) candidates, of which :"
echo "# - $(($(wc -l ${pfx}/results-filtered_minus.txt|cut -d" " -f1)-1)) minus-strand HQ hits"
echo "# - $(($(wc -l ${pfx}/results-filtered_plus.txt|cut -d" " -f1)-1)) plus-strand HQ hits"

# create BED from filtered
gawk -v tile=${tile} 'BEGIN{FS="\t"; OFS="\t"}{if (NR>1 && $6==0 && $11>0.3) {split($3, a, ":"); print a[1], a[2]-1, a[2]+22, tile"_dwn:"$1, $11, $4, a[2]-1, a[2]+22, "0,0,255"}}' \
  ${pfx}/results.txt | sort -k1V,1 -k2n,2 -k 3n,3 \
  > ${pfx}/filtered_guides.bed
