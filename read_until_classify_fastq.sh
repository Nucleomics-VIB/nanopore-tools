#!/bin/bash

# script: read_until_classify_fastq.sh
#
# Aim: Identify fastq_pass reads based on read_until classification
# print read ID, read size, and classification to a text file for R plotting
# print reads in fastq format to new files with classification in the name

usage='# Usage: read_until_classify_fastq.sh
# -r <read_until_file (eg. other_reports/read_until_<FCID>_<runID>.csv)>
# -f <merged fastq file (can be gzipped)>
# script version '${version}'
# [-h for this help]'

while getopts "r:f:h" opt; do
  case $opt in
    r) optr=${OPTARG} ;;
    f) optf=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
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
  echo "# provide the path to a merged fastq file (can be gzipped)"
  echo "${usage}" >&2
  exit 0
fi

# store results in to folder with prefix
outpath="$(dirname ${optf})_split"
outpfx="$(basename ${optf} .fq.gz)"
mkdir -p ${outpath}

bioawk -c fastx -v rufile="${optr}" -v outp="${outpath}/${outpfx}" '
BEGIN{
  FS=","; 
  OFS="\t";
  countfile=outp"_read_until_counts.txt";
  while((getline line <rufile)>0) {
    split(line, f, ","); decision[f[6]]=f[8]};
  print "read_id","read_length","classified" > countfile
  }
{
  split($name, info, " ");
  if (decision[info[1]] == "") decision[info[1]]="absent";
  fastqfile=outp"_"decision[info[1]]".fq";
  print "@"$name" "$comment"\n"$seq"\n+\n"$qual >> fastqfile;
  print info[1], length($seq), decision[info[1]] >> countfile
  }' ${optf}
