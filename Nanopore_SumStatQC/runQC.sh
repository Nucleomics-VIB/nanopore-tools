#!/bin/bash
# run ONT QC on basecalling summary files(s)
# runQC <data.bz2> <FlowCell-ID> <barcodes.bz2>
# edit additional variables in the created config.yaml

working=$(pwd)
script_base=$(dirname $0)
script="Nanopore_SumStatQC.Rmd"

mkdir -p ${working}/RawData

echo "# processing $1"
echo "# flow-cell.ID: $2"
infile="${working}/RawData/$(basename $1).bz2"
#bgzip -c $1 > ${infile}

if [[ ! -z "$3" ]]; then
echo "# processing barcodes from $3"
bcfile="${working}/RawData/$(basename $3).bz2"
#bgzip -c $3 > ${bcfile}
else
bcfile=""
fi

echo "expRef: \"My fantastic experiment\"
Organism: \"NA\"
inputFile: \"${infile}\"
barcodeFile: \"${bcfile}\"
basecaller: \"Guppy 2.0.10 - X5-GPU\"
flowcellId: \"$2\"
flowcellType: \"FLO-MIN106\"
seqKit: \"SQK-LSK108\"
tutorialText: FALSE
sampleHours: 48
sampleIntervalMinutes: 60" | tee config.yaml

cp config.yaml $(dirname $0)/

echo "# edit $(dirname $0)/config.yaml"
echo "# then run with: R --slave -e 'rmarkdown::render(\"${script_base}/${script}\", \"html_document\")' "
