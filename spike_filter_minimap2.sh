#!/bin/bash

# script: spike_filter_minimap2.sh
# search for lambda phage reads in MinIon data
#  align all fastq reads to the lambda reference (provided fasta) with minimap2
#  extract reads aligned / not-aligned from SAM (samtools view -F/-f 4)
#  convert aligned / not-aligned SAM data back to fastq (samtools fastq)
#
# Stephane Plaisance (VIB-NC) 2018/11/01; v1.1
# added min mapping quality filter

# required (tested versions)
# samtools (1.x, htslib)
# minimap2 (2.11.x)

read -d '' usage <<- EOF
Usage: spike_filter_minimap2.sh 
#   -i <nanopore_reads.fastq (required)>
#   -r <spiked reference (lambda or any other spiked genome, required)>
#   -q <minimal mapping quality to consider a reads aligment (default to 0)>
#   -t <threads to be used for alignment (default to 8)>
#   -x <full path to minimap2 if not in PATH>
#   -s <keep only spiked reads instead (reverse-mode)>
#   -h <show this help>
EOF

while getopts "i:r:q:t:x:sh" opt; do
	case $opt in
		i)
		  infile=${OPTARG}
		  ;;
		r)
		  reference=${OPTARG}
		  ;;
		q)
		  mapqual=${OPTARG}
		  ;;
		t)
		  threads=${OPTARG}
		  ;;
		x)
		  exe=${OPTARG}
		  ;;
		s)
		  keepspiked=1
		  ;;
		h)
		  echo "${usage}"
		  exit 0
		  ;;
		\?)
		  echo "${usage}"
		  exit 1
		  ;;
		*)
		  echo "${usage}" >&2
		  exit 1
		  ;;
	esac
done

# functions
function testexecutable ()
{
if [[ ! -x "$1" ]]
then
	echo "! # ${1} is not executable or absent"
	echo "${usage}"
	exit 1
else
	return 0
fi
}

# choose minimap2 exe
minimap_2=$(which minimap2)
minimap_exe=${exe:-${minimap_2}}

testexecutable "${minimap_exe}"

# test if minimal arguments were provided
if [ -z "${infile}" ]; then
echo "# no read data provided!"
echo "${usage}"
exit 1
fi

# test if minimal arguments were provided
if [ -z "${reference}" ]; then
echo "# no reference provided!"
echo "${usage}"
exit 1
fi

# minimal mapping quality to consider read alignment
minq=${mapqual:-0}

# filter out or keep the spiked alignments
if [ -n "${keepspiked}" ]; then
samfilter="-F 4"
prefix="spiked"
else
samfilter="-f 4"
prefix="mm2_filtered"
fi

# threads
nthr=${threads:-8}

# define filehandlers
outpath=$(dirname $infile);
outbase=$(basename ${infile});
outfile=${outbase%.f*};

##################################################
# aligning all reads to the spiked genome sequence

# echo "# (re-)creating reference index then aligning reads"
cmd1="${minimap_exe} \
	-t ${nthr} \
	-ax map-ont \
	-a \
	${reference} \
	${infile}"

echo "# aligning: ${cmd1} using ${minimap_exe}"

############################################
# filtering and creating de-spiked read file

# fetching unmapped reads from alignments
# samtools view ${samfilter} ${outpath}${outfile}_alignments.sam > ${outpath}/unmapped.sam
# converting back to fastq
# samtools fastq ${outpath}/unmapped.sam | bgzip -c > ${outpath}/${prefix}_${outbase%.gz}.gz

cmd2="samtools view ${samfilter} -bSq ${minq} - \
	| samtools fastq - \
	| bgzip -c > ${outpath}/${prefix}_${outbase%.gz}.gz"
echo "# filtering: ${cmd2}"

##################################################################
# merging both commands to limit output to the de-spiked read file
mrgcmd="${cmd1} | ${cmd2}"
echo "# merged command: ${mrgcmd}"
eval ${mrgcmd}

# count reads
incnt=$(zgrep -c "^@" ${infile})
outcnt=$(zgrep -c "^@" ${outpath}/${prefix}_${outbase%.gz}.gz)
difcnt=$(echo "${incnt}-${outcnt}" |bc)

echo
echo "# reads in the input : ${incnt}"
echo "# filtered reads     : ${outcnt}"
echo "# omitted reads      : ${difcnt}"
exit 0

################################################################################################################################
# Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]
# Options:
#   Indexing:
#     -H           use homopolymer-compressed k-mer (preferrable for PacBio)
#     -k INT       k-mer size (no larger than 28) [15]
#     -w INT       minizer window size [10]
#     -I NUM       split index for every ~NUM input bases [4G]
#     -d FILE      dump index to FILE []
#   Mapping:
#     -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [0.0002]
#     -g NUM       stop chain enlongation if there are no minimizers in INT-bp [5000]
#     -G NUM       max intron length (effective with -xsplice; changing -r) [200k]
#     -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]
#     -r NUM       bandwidth used in chaining and DP-based alignment [500]
#     -n INT       minimal number of minimizers on a chain [3]
#     -m INT       minimal chaining score (matching bases minus log gap penalty) [40]
#     -X           skip self and dual mappings (for the all-vs-all mode)
#     -p FLOAT     min secondary-to-primary score ratio [0.8]
#     -N INT       retain at most INT secondary alignments [5]
#   Alignment:
#     -A INT       matching score [2]
#     -B INT       mismatch penalty [4]
#     -O INT[,INT] gap open penalty [4,24]
#     -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]
#     -z INT[,INT] Z-drop score and inversion Z-drop score [400,200]
#     -s INT       minimal peak DP alignment score [80]
#     -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]
#   Input/Output:
#     -a           output in the SAM format (PAF by default)
#     -Q           don't output base quality in SAM
#     -L           write CIGAR with >65535 ops at the CG tag
#     -R STR       SAM read group line in a format like '@RG\tID:foo\tSM:bar' []
#     -c           output CIGAR in PAF
#     --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]
#     --MD         output the MD tag
#     --eqx        write =/X CIGAR operators
#     -Y           use soft clipping for supplementary alignments
#     -t INT       number of threads [3]
#     -K NUM       minibatch size for mapping [500M]
#     --version    show version number
#   Preset:
#     -x STR       preset (always applied before other options; see minimap2.1 for details) []
#                  - map-pb/map-ont: PacBio/Nanopore vs reference mapping
#                  - ava-pb/ava-ont: PacBio/Nanopore read overlap
#                  - asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
#                  - splice: long-read spliced alignment
#                  - sr: genomic short-read mapping
#
# See `man ./minimap2.1' for detailed description of these and other advanced command-line options.
################################################################################################################################
