#!/bin/bash
# run ONT QC on basecalling summary files(s)

working=$(pwd)
script_base=$(dirname $0)
script="Nanopore_SumStatQC.Rmd"

mkdir -p ${working}/RawData

echo "# processing $1"
echo "# flow-cell.ID: $2"
infile="${working}/RawData/$(basename $1).bz2"
bgzip -c $1 > ${infile}

if [[ ! -z "$3" ]]; then
echo "# processing barcodes from $3"
bcfile="${working}/RawData/$(basename $3).bz2"
bgzip -c $3 > ${bcfile}
else
bcfile=""
fi

echo "inputFile: \"${infile}\"
barcodeFile: \"${bcfile}\"
basecaller:  \"Guppy 2.0.10 GPU\"
flowcellId:  \"$2\"
tutorialText: FALSE" > config.yaml

cp config.yaml $(dirname $0)/

R --slave -e "rmarkdown::render(\"${script_base}/${script}\", \"html_document\")"
