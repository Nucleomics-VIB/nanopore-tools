#!/bin/bash
# script name: run_chopchop_single.sh
# enrich a region with ONT Cas9
# run chopchop prediction upstream and downstream of the region
# REM: chopchop should have the reference genome
# produce a BED file with filtered sgRNAs for IGV visualisation
# REM: run the tile version of this script if teh region is larger tnan 30kb
#
# Stephane Plaisance (VIB-NC) 2019/06/12; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2019_06_12"

usage='# Usage: run_chopchop_single.sh
# -r <region (eg. chr1:16920000-16940000)>
# -g <reference genome used by ChopChop (default to hg38)>
# -w <prediction window (default to 3000bps either side)>
# -m <median DNA fragment size (default 30000)>
# script version '${version}'
# [-h for this help]'

while getopts "r:g:w:m:h" opt; do
  case $opt in
    r) optr=${OPTARG} ;;
    g) optg=${OPTARG} ;;
    w) optw=${OPTARG} ;;
    m) optm=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check executables present
declare -a arr=( "chopchop.py" "bedtools" "tee" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# add all logs to a global log file
datetag=$(date +%s)
basedir="run_chopchop_single_${datetag}_data"
mkdir -p ${basedir}
exec &> >(tee -a "${basedir}/run_chopchop_single.log")

# provide the region to be captured (start to end shorter than 20kb)
region=${optr:-"chr1:16920000-16940000"}
chr=${region%%:*}
coord=${region#*:}
start=${coord%-*}
end=${coord#*-}

# define window and step sizes
window=${optw:-3000}

# check width
width=$((${end}-${start}))
median=${optm:-30000}
if [ "${width}" -gt "${median}" ]; then
echo "# The region is too large (${width} bps) for a single sgRNA pair, consider tiling !"
exit 0
fi

# define reference genome or take default hg38 (should be installed on your machine)
refg=${optg:-"hg38"}

# default primer3 settings
primer3args='PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60'
backbone="AGGCTAGTCCGT"
scoring="DOENCH_2014"

# tile limits
tile="${chr}_${start}_${end}"

# create index name and folders for results
pfx="${basedir}/${chr}_${start}_${end}_results"
mkdir -p ${pfx}

# clear ${pfx}/filtered_guides.bed
cat /dev/null > ${pfx}/filtered_guides.bed

# define upstream & downstream for this tile
upstream="${chr}:$((${start}-${window}))-${start}"
downstream="${chr}:${end}-$((${end}+${window}))"

# search 3kb upstream
echo "# predicting ${window}bps upstream of ${tile}"
chopchop.py -J \
		-P \
		-T 1 \
		-M NGG \
		--maxMismatches 3 \
		-g 20 \
		--scoringMethod ${scoring} \
		-f NN \
		--backbone ${backbone} \
		--replace5P GG \
		-G ${refg} \
		-t WHOLE \
		-n N \
		-R 4 \
		-3 ${primer3args} \
		-A 290 \
		-a 20 \
		--rm1perfOff \
		-o ${pfx}/upstream/ \
		--filterGCmin 40 \
		--filterGCmax 80 \
		--filterSelfCompMax 0 \
		--BED \
		--GenBank \
		-Target ${upstream} \
		> ${pfx}/upstream_results.txt 2> ${pfx}/upstream_python.err

# filter best hits on + strand   
gawk 'BEGIN{FS="\t"; OFS="\t"}{if( (NR==1) || ($4~/+/ && $6==0 && $11>0.3)) print $0}' \
	${pfx}/upstream_results.txt | column -s $'\t' -t \
	> ${pfx}/upstream_results-filtered.txt
echo "# => found $(($(wc -l ${pfx}/upstream_results.txt|cut -d" " -f1)-1)) candidates of which $(($(wc -l ${pfx}/upstream_results-filtered.txt|cut -d" " -f1)-1)) HQ hits"

# create BED from filtered
gawk -v tile=${tile} 'BEGIN{FS="\t"; OFS="\t"}{if (NR>1 && $4~/+/ && $6==0 && $11>0.3) {split($3, a, ":"); print a[1], a[2], a[2]+23, tile"_ups:"$1, $11, $4, a[2], a[2]+23, "255,0,0"}}' \
${pfx}/upstream_results.txt \
> ${pfx}/filtered_guides.bed

# search 3kb downstream
echo "# predicting ${window}bps downstream of ${tile}"
chopchop.py -J \
		-P \
		-T 1 \
		-M NGG \
		--maxMismatches 3 \
		-g 20 \
		--scoringMethod ${scoring} \
		-f NN \
		--backbone ${backbone} \
		--replace5P GG \
		-G ${refg} \
		-t WHOLE \
		-n N \
		-R 4 \
		-3 ${primer3args} \
		-A 290 \
		-a 20 \
		--rm1perfOff \
		-o ${pfx}/downstream/ \
		--filterGCmin 40 \
		--filterGCmax 80 \
		--filterSelfCompMax 0 \
		--BED \
		--GenBank \
		-Target ${downstream} \
		> ${pfx}/downstream_results.txt 2> ${pfx}/downstream_python.err

# filter best hits on + strand   
gawk 'BEGIN{FS="\t"; OFS="\t"}{if( (NR==1) || ($4~/-/ && $6==0 && $11>0.3)) print $0}' \
	${pfx}/downstream_results.txt | column -s $'\t' -t \
	> ${pfx}/downstream_results-filtered.txt
echo "# => found $(($(wc -l ${pfx}/downstream_results.txt|cut -d" " -f1)-1)) candidates of which $(($(wc -l ${pfx}/downstream_results-filtered.txt|cut -d" " -f1)-1)) HQ hits"

# create BED from filtered
gawk -v tile=${tile} 'BEGIN{FS="\t"; OFS="\t"}{if (NR>1 && $4~/-/ && $6==0 && $11>0.3) {split($3, a, ":"); print a[1], a[2], a[2]+23, tile"_dwn:"$1, $11, $4, a[2], a[2]+23, "0,0,255"}}' \
${pfx}/downstream_results.txt \
>> ${pfx}/filtered_guides.bed
	
# merge all filtered-bed to one with track lines for IGV
mergedbed="${basedir}/tmp.bed"
sortedbed="${basedir}/filtered_guides.bed"
echo "browser position ${region}" > ${mergedbed}
echo -e "track name=CHOPCHOP description=${region} visibility=\"pack\" itemRgb=\"On\"" >> ${mergedbed}
echo -e "${chr}\t${start}\t${end}\tROI\t0\t+\t${start}\t${end}\t0,255,0" >> ${mergedbed}
# add all filtered bed files
find ${basedir} -name filtered_guides.bed -exec cat {} >> ${mergedbed} \;
bedtools sort -header -i ${mergedbed} > ${sortedbed} && rm ${basedir}/tmp.bed
