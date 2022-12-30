#!/bin/bash
# script name: run_guppy_sup.sh
# (re)basecall fast5 data
# optionally demultiplex on the go
# REM: runs on the gridion to take advantage of the GPU
#
# Stephane Plaisance (VIB-NC) 2021/07/12; v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

version="1.0, 2021_07_12"

usage='# Usage: run_guppy_sup.sh
# -i <fast5 input path>
# -s <save path (default to "guppy_sup_out")>
# -c <config file (default to "dna_r9.4.1_450bps_sup.cfg")>
# -b <barcode kit (eg: SQK-RBK004; default without demultiplexing)>
# -p <allow_inferior_barcodes (default off)>
# script version '${version}'
# -l get the list of all current config files
# [-h for this help]'

while getopts "i:s:c:b:plh" opt; do
  case $opt in
    i) opti=${OPTARG} ;;
    s) opts=${OPTARG} ;;
    c) optc=${OPTARG} ;;
    b) optb=${OPTARG} ;;
    p) optp=true ;;
    l) echo "# list of all current workflows"; \
       guppy_basecaller --print_workflows; \
       exit 0 ;;
    h) echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *) echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# check executables present
declare -a arr=( "guppy_basecaller" )
for prog in "${arr[@]}"; do
$( hash ${prog} 2>/dev/null ) || ( echo "# required ${prog} not found in PATH"; exit 1 )
done

# test if minimal arguments were provided
if [ -z "${opti}" ]; then
echo "# no fast5 path provided!"
echo "${usage}"
exit 1
fi

# default arguments when absent
outfolder=${opts:-"guppy_sup_out"}
config=${optc:-"dna_r9.4.1_450bps_sup.cfg"}

#####################
# optional arguments
#####################

# only demux if optb is provided
barcode_kit=""

if [ -n "${optb}" ]; then
  barcode_kit="--barcode_kits ${optb}"
fi

# be more permissive for barcodes
be_permissive=""

if [ -n "${optp}" ]; then
  be_permissive="--allow_inferior_barcodes"
fi

# build guppy command
cmd="guppy_basecaller \
  -i ${opti} \
  -r \
  -s ${outfolder} \
  -c ${config} \
  -x 'cuda:0' \
  --compress_fastq \
  ${barcode_kit} \
  ${be_permissive}"

echo "# command: ${cmd}"
echo
eval ${cmd}
