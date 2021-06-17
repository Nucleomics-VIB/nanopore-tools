#!/bin/bash

# script: read_until_classify_fastq.sh
#
# Aim: Identify fastq_pass reads based on read_until classification
# print read ID, read size, and classification to a text file for R plotting
# print reads in fastq format to new files with classification in the name

usage='# Usage: read_until_classify_fastq.sh
# -r <read_until_file (eg. read_until_FAQ08767_7379d056.csv)>
# -f <fastq folder>
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
  echo "# provide a read_until summary file (eg read_until_FAQ08767_7379d056.csv)"
  echo "${usage}" >&2
  exit 0
fi

if [ -z "${optf}" ]
then
  echo "# provide a path containing fastq files (eg fastq_pass)"
  echo "${usage}" >&2
  exit 0
fi

bioawk -c fastx -v rufile="${optr}" '
BEGIN{
  FS=","; 
  OFS="\t";
  while((getline line <rufile)>0) {
    split(line, f, ","); decision[f[6]]=f[8]};
 print "read_id","read_length","classified" > "read_until_counts.txt"
 }
{
split($name, info, " ");
if (decision[info[1]] == "") decision[info[1]]="absent";
fastqfile=decision[info[1]]".fq";
print "@"$name" "$comment"\n"$seq"\n+\n"$qual >>fastqfile;
print info[1], length($seq), decision[info[1]] >>"read_until_counts.txt"
}' ${optf}/*.fastq
