#!/bin/bash

# script: read_until_subset_fast5.sh
#
# Aim: Subset stop_receiving fast5_pass reads
# get stop_receiving read IDs from other_reports/read_until_<FCID>_<runID>.csv
# feed the obtained list to fast5_subset (from ont_fast5_api in conda env)
# save to new fast5 files

###############
# Requirements:
###############

# miniconda installed
# conda environment ont_github created
# ont_fast5_api installed (https://github.com/nanoporetech/ont_fast5_api)

######################

usage='# Usage: read_until_subset_fast5.sh
# -r <path to the read_until_<FCID>_<runID>.csv file>
# -f <path to the fast5_pass folder>
# -t <threads (default to 1)>
# script version '${version}'
# [-h for this help]'

# activate conda env or die
myenv="ont_github"
source /etc/profile.d/conda.sh || (echo "conda profile not found!; exit 1");

conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

while getopts "r:f:t:h" opt; do
  case $opt in
    r) optr=${OPTARG} ;;
    f) optf=${OPTARG} ;;
    t) optt=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check executables present
declare -a arr=( "fast5_subset" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# check INPUTS
if [ -z "${optr}" ]
then
  echo "# provide the path to the read_until summary file (eg other_reports/read_until_<FCID>_<runID>.csv)"
  echo "${usage}" >&2
  exit 0
fi

if [ -z "${optf}" ]
then
  echo "# provide the path to the fast5_pass folder"
  echo "${usage}" >&2
  exit 0
fi

if [ -z "${optt}" ]
then
  thr=1
else
  thr=${optt}
fi

# create folder to store subset
outfolder="fast5_stop_receiving"
mkdir -p ${outfolder}

# create list of stop_receiving from other_reports/read_until_<FCID>_<runID>.csv
awk 'BEGIN{FS=","; OFS="\t"}{if ($8=="stop_receiving") print $6}' "${optr}" \
  > ${outfolder}/stop_receiving_ids.txt

# extract fast5 (20k batches) using fast5_subset from the conda env
fast5_subset -i ${optf} \
  -s ${outfolder} \
  -l ${outfolder}/stop_receiving_ids.txt \
  -n 20000 \
  -t ${thr} \
  -r
