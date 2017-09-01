#!/bin/bash

# script: spike_clean.sh
# search for lambda phage reads in MinIon data
#  align all fastq reads to the lambda reference (provided fasta) with graphmap
#  extract reads not-aligned from SAM (samtools view -f 4)
#  convert not-aligned SAM data back to fastq (samtools fastq)
#
# Stephane Plaisance (VIB-NC) 2017/09/01; v1.0

# required (tested versions)
# samtools (1.x, htslib)
# graphmap (0.5.x)

read -d '' usage <<- EOF
Usage: spike_clean.sh 
#   -i <nanopore_reads.fastq (required)>
#   -r <spiked reference (lambda or any other spiked genome, required)>
#   -t <graphmap alignment threshold (default to 1e-0)>
#   -e <read error rate (default to 0.45)>
#   -C <the spiked genome is circular (default OFF = linear)>
#   -h <show this help>
EOF

while getopts "i:r:t:e:Ch" opt; do
  case $opt in
    i)
      infile=${OPTARG}
      ;;
    r)
      reference=${OPTARG}
      ;;
	t)
	  threshold=${OPTARG}
	  ;;
	e)
	  errorrate=${OPTARG}
	  ;;
	C)
	  circular=1
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

# threshold for Evalue to call a read unmapped (1e-100!)
thresh=${threshold:-"1e-0"}

# Approximate error rate of the input read sequences
errate=${errorrate:-"0.45"}

# the spiked genome is circular
if [ -n "${circular}" ]; then
circ="-C"
fi

# threads
nthr=8

# define filehandlers
outpath=$(dirname $infile);
outbase=$(basename ${infile});
outfile=${outbase%.f*};

# echo "# (re-)creating reference index then aligning reads"
cmd="graphmap align --threads ${nthr}\
	--auto-rebuild-index ${circ}\
	--error-rate ${errate} \
	--evalue ${thresh} \
	-r ${reference} \
	-d ${infile} \
	-o ${outpath}"/"${outfile}_alignments.sam"
echo "# ${cmd}"
eval ${cmd}

# fetching unmapped reads from alignments
# samtools view -f 4 ${outfile}_alignments.sam > unmapped.sam
# converting back to fastq
# samtools fastq unmapped.sam | bgzip -c > cleaned_${outbase%.gz}.gz

cmd="samtools view -f 4 ${outpath}"/"${outfile}_alignments.sam | samtools fastq - | bgzip -c > ${outpath}"/"cleaned_${outbase%.gz}.gz"
echo "# ${cmd}"
eval ${cmd}

exit 0

################################################################################################################################
# GraphMap - A very accurate and sensitive long-read, high error-rate sequence mapper
# GraphMap Version: v0.5.2
# Build date: Aug 31 2017 at 21:08:57
# 
# GraphMap (c) by Ivan Sovic, Mile Sikic and Niranjan Nagarajan
# GraphMap is licensed under The MIT License.
# 
# Affiliations: Ivan Sovic (1, 3), Mile Sikic (2), Niranjan Nagarajan (3)
#   (1) Ruder Boskovic Institute, Zagreb, Croatia
#   (2) University of Zagreb, Faculty of Electrical Engineering and Computing
#   (3) Genome Institute of Singapore, A*STAR, Singapore
# 
# 
# Usage:
# 	graphmap [options] -r <reference_file> -d <reads_file> -o <output_sam_path>
# 
# 
# Usage:
#   graphmap [options]
# 
# Options
#   Input/Output options:
#     -r, --ref                STR   Path to the reference sequence (fastq or fasta).
#     -i, --index              STR   Path to the index of the reference sequence. If not specified, index is generated in
#                                    the same folder as the reference file, with .gmidx extension. For non-parsimonious
#                                    mode, secondary index .gmidxsec is also generated.
#     -d, --reads              STR   Path to the reads file.
#     -o, --out                STR   Path to the output file that will be generated.
#         --gtf                STR   Path to a General Transfer Format file. If specified, a transcriptome will be built
#                                    from the reference sequence and used for mapping. Output SAM alignments will be in
#                                    genome space (not transcriptome).
#     -K, --in-fmt             STR   Format in which to input reads. Options are:
#                                     auto  - Determines the format automatically from file extension.
#                                     fastq - Loads FASTQ or FASTA files.
#                                     fasta - Loads FASTQ or FASTA files.
#                                     gfa   - Graphical Fragment Assembly format.
#                                     sam   - Sequence Alignment/Mapping format. [auto]
#     -L, --out-fmt            STR   Format in which to output results. Options are:
#                                     sam  - Standard SAM output (in normal and '-w overlap' modes).
#                                     m5   - BLASR M5 format. [sam]
#     -I, --index-only          -    Build only the index from the given reference and exit. If not specified, index will
#                                    automatically be built if it does not exist, or loaded from file otherwise. [false]
#         --rebuild-index       -    Always rebuild index even if it already exists in given path. [false]
#         --auto-rebuild-index  -    Rebuild index only if an existing index is of an older version or corrupt. [false]
#     -u, --ordered             -    SAM alignments will be output after the processing has finished, in the order of
#                                    input reads. [false]
#     -B, --batch-mb           INT   Reads will be loaded in batches of the size specified in megabytes. Value <= 0 loads
#                                    the entire file. [1024]
# 
#   General-purpose pre-set options:
#     -x, --preset             STR   Pre-set parameters to increase sensitivity for different sequencing technologies.
#                                    Valid options are:
#                                     illumina  - Equivalent to: '-a gotoh -w normal -M 5 -X 4 -G 8 -E 6'
#                                     overlap   - Equivalent to: '-a anchor -w normal --overlapper --evalue 1e0
#                                    --ambiguity 0.50 --secondary'
#                                     sensitive - Equivalent to: '--freq-percentile 1.0 --minimizer-window 1'
# 
#   Alignment options:
#     -a, --alg                STR   Specifies which algorithm should be used for alignment. Options are:
#                                     sg       - Myers' bit-vector approach. Semiglobal. Edit dist. alignment.
#                                     sggotoh       - Gotoh alignment with affine gaps. Semiglobal.
#                                     anchor      - anchored alignment with end-to-end extension.
#                                                   Uses Myers' global alignment to align between anchors.
#                                     anchorgotoh - anchored alignment with Gotoh.
#                                                   Uses Gotoh global alignment to align between anchors. [anchor]
#     -w, --approach           STR   Additional alignment approaches. Changes the way alignment algorithm is applied.
#                                    Options are:
#                                     normal         - Normal alignment of reads to the reference.
#                                     (Currently no other options are provided. This is a placeholder for future features,
#                                    such as cDNA mapping) [normal]
#         --overlapper          -    Perform overlapping instead of mapping. Skips self-hits if reads and reference files
#                                    contain same sequences, and outputs lenient secondary alignments. [false]
#         --no-self-hits        -    Similar to overlapper, but skips mapping of sequences with same headers. Same
#                                    sequences can be located on different paths, and their overlap still skipped. [false]
#     -M, --match              INT   Match score for the DP alignment. Ignored for Myers alignment. [5]
#     -X, --mismatch           INT   Mismatch penalty for the DP alignment. Ignored for Myers alignment. [4]
#     -G, --gapopen            INT   Gap open penalty for the DP alignment. Ignored for Myers alignment. [8]
#     -E, --gapext             INT   Gap extend penalty for the DP alignment. Ignored for Myers alignment. [6]
#     -z, --evalue             FLT   Threshold for E-value. If E-value > FLT, read will be called unmapped. If FLT < 0.0,
#                                    thredhold not applied. [1e0]
#     -c, --mapq               INT   Threshold for mapping quality. If mapq < INT, read will be called unmapped. [1]
#         --extcigar            -    Use the extended CIGAR format for output alignments. [false]
#         --no-end2end          -    Disables extending of the alignments to the ends of the read. Works only for
#                                    anchored modes. [false]
#         --max-error          FLT   If an alignment has error rate (X+I+D) larger than this, it won't be taken into
#                                    account. If >= 1.0, this filter is disabled. [1.0]
#         --max-indel-error    FLT   If an alignment has indel error rate (I+D) larger than this, it won't be taken into
#                                    account. If >= 1.0, this filter is disabled. [1.0]
# 
#   Algorithmic options:
#     -k                       INT   Graph construction kmer size. [6]
#     -l                       INT   Number of edges per vertex. [9]
#     -A, --minbases           INT   Minimum number of match bases in an anchor. [12]
#     -e, --error-rate         FLT   Approximate error rate of the input read sequences. [0.45]
#     -g, --max-regions        INT   If the final number of regions exceeds this amount, the read will be called
#                                    unmapped. If 0, value will be dynamically determined. If < 0, no limit is set. [0]
#     -q, --reg-reduce         INT   Attempt to heuristically reduce the number of regions if it exceeds this amount.
#                                    Value <= 0 disables reduction but only if param -g is not 0. If -g is 0, the value of
#                                    this parameter is set to 1/5 of maximum number of regions. [0]
#     -C, --circular            -    Reference sequence is a circular genome. [false]
#     -F, --ambiguity          FLT   All mapping positions within the given fraction of the top score will be counted for
#                                    ambiguity (mapping quality). Value of 0.0 counts only identical mappings. [0.02]
#     -Z, --secondary           -    If specified, all (secondary) alignments within (-F FLT) will be output to a file.
#                                    Otherwise, only one alignment will be output. [false]
#     -P, --double-index        -    If false, only one gapped spaced index will be used in region selection. If true,
#                                    two such indexes (with different shapes) will be used (2x memory-hungry but more
#                                    powerful for very high error rates). [false]
#         --min-bin-perc       FLT   Consider only bins with counts above FLT * max_bin, where max_bin is the count of
#                                    the top scoring bin. [0.75]
#         --bin-step           FLT   After a chunk of bins with values above FLT * max_bin is processed, check if there
#                                    is one extremely dominant region, and stop the search. [0.25]
#         --min-read-len       INT   If a read is shorter than this, it will be marked as unmapped. This value can be
#                                    lowered if the reads are known to be accurate. [80]
#         --minimizer-window   INT   Length of the window to select a minimizer from. If equal to 1, minimizers will be
#                                    turned off. [5]
#         --freq-percentile    FLT   Filer the (1.0 - value) percent of most frequent seeds in the lookup process. [0.99]
#         --fly-index           -    Index will be constructed on the fly, without storing it to disk. If it already
#                                    exists on disk, it will be loaded unless --rebuild-index is specified. [false]
# 
#   Other options:
#     -t, --threads            INT   Number of threads to use. If '-1', number of threads will be equal to min(24, num_cores/2). [-1]
#     -v, --verbose            INT   Verbose level. If equal to 0 nothing except strict output will be placed on stdout. [5]
#     -s, --start              INT   Ordinal number of the read from which to start processing data. [0]
#     -n, --numreads           INT   Number of reads to process per batch. Value of '-1' processes all reads. [-1]
#     -h, --help                -    View this help. [false]
# 
#   Debug options:
#     -y, --debug-read         INT   ID of the read to give the detailed verbose output. [-1]
#     -Y, --debug-qname        STR   QNAME of the read to give the detailed verbose output. Has precedence over -y. Use
#                                    quotes to specify.
#     -b, --verbose-sam        INT   Helpful debug comments can be placed in SAM output lines (at the end). Comments can
#                                    be turned off by setting this parameter to 0. Different values increase/decrease
#                                    verbosity level.
#                                    0 - verbose off
#                                    1 - server mode, command line will be omitted to obfuscate paths.
#                                    2 - umm this one was skipped by accident. The same as 0.
#                                    >=3 - detailed verbose is added for each alignment, including timing measurements and
#                                    other.
#                                    4 - qnames and rnames will not be trimmed to the first space.
#                                    5 - QVs will be omitted (if available). [0]
################################################################################################################################
