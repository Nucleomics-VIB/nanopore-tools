#!/bin/bash
# script name: run_chopchop.sh
# run chopchop prediction upstream and downstream of a region
#
# Stephane Plaisance (VIB-NC) 2019/06/12; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2019_06_12"

usage='# Usage: run_chopchop.sh
# -r <region (eg. chr1:10000-20000)> 
# -w <width for prediction at both ends (default to 3kb)>
# -m <median DNA fragment size (default 30000)>
# script version '${version}'
# [-h for this help]'

while getopts "r:w:m:h" opt; do
  case $opt in
    r) optr=${OPTARG} ;;
    w) optw=${OPTARG} ;;
    m) optm=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# provide the region to be captured (start to end shorter than 20kb)
region=${optr:-"chr1:16740273-16972964"}
chr=${region%%:*}
coord=${region#*:}
start=${coord%-*}
end=${coord#*-}

# check width
width=$((${end}-${start}))
median=${optm:-30000}
if [ "${width}" -gt "${median}" ]; then
echo "# The region is too large (${width} bps) for a single sgRNA pair, consider tiling !"
exit 0
fi

# define upstream and downstream regions
window=${optw:-3000}
upstream="${chr}:$((${start}-${window}))-${start}"
downstream="${chr}:${end}-$((${end}+${window}))"

# default primer3 settings
primer3args='PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60'

mkdir -p {upstream,downstream}

# search 3kb upstream
chopchop.py -J \
        -P \
        -T 1 \
        -M NGG \
        --maxMismatches 3 \
        -g 20 \
        --scoringMethod DOENCH_2014 \
        -f NN \
        --backbone AGGCTAGTCCGT \
        --replace5P GG \
        -G hg38 \
        -t WHOLE \
        -n N \
        -R 4 \
        -3 ${primer3args} \
        -A 290 \
        -a 20 \
        --rm1perfOff \
        -o upstream/ \
        --filterGCmin 40 \
        --filterGCmax 80 \
        --filterSelfCompMax 0 \
        --BED \
        --GenBank \
        -Target ${upstream} \
        > upstream/upstream_results.txt 2> upstream/python.err

# filter best hits on + strand   
gawk 'BEGIN{FS="\t"; OFS="\t"}{if( (NR==1) || ($4~/+/ && $6==0 && $11>0.3)) print $0}' \
	upstream/upstream_results.txt | column -s $'\t' -t \
	> upstream/upstream_results-filtered.txt

# search 3kb downstream
chopchop.py -J \
        -P \
        -T 1 \
        -M NGG \
        --maxMismatches 3 \
        -g 20 \
        --scoringMethod DOENCH_2014 \
        -f NN \
        --backbone AGGCTAGTCCGT \
        --replace5P GG \
        -G hg38 \
        -t WHOLE \
        -n N \
        -R 4 \
        -3 ${primer3args} \
        -A 290 \
        -a 20 \
        --rm1perfOff \
        -o downstream/ \
        --filterGCmin 40 \
        --filterGCmax 80 \
        --filterSelfCompMax 0 \
        --BED \
        --GenBank \
        -Target ${downstream} \
        > downstream/downstream_results.txt 2> downstream/python.err

# filter best hits on + strand   
gawk 'BEGIN{FS="\t"; OFS="\t"}{if( (NR==1) || ($4~/-/ && $6==0 && $11>0.3)) print $0}' \
	downstream/downstream_results.txt | column -s $'\t' -t \
	> downstream/downstream_results-filtered.txt
