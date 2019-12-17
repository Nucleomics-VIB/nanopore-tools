#!/usr/bin/env bash
# script rundata2arch.sh
# author:Stephane Plaisance (VIB-NC), 2019-12-17, v1.0

# create backup files from a gridion run

if [ $# -eq 0 ]; then
    echo "# Please provide the path to the run folder as argument"
    exit 1
fi

# get runfolder as argument 1
runf=${1}
pfx=$(basename ${runf})

# check for valid runfolder
fqp=$(find ${runf} -type d -name fastq_pass -printf '%d\n' | sort -rn | head -1)

if (( ${fqp} != 1 )); then
	echo "# please provide a valid run-folder as input"
	echo "# current path is:"
	ls -lah ${runf}
	exit
fi

# get enclosing folder name
projf=$(basename $(dirname $(find ../ -name "${runf}")))

# create run folder_archive for transfer
tmpf=${projf}_data
mkdir -p ${tmpf}

# copy and rename gridion run-report
echo "# copying the ONT run report"
cp ${runf}/report.pdf ${tmpf}/${pfx}_report.pdf

# archive gridion run stat files
echo "# archiving run stat files to ${pfx}_run_info.tgz"
find ${runf} -maxdepth 1 -type f -exec tar -czvf ${tmpf}/${pfx}_run_info.tgz {} \;

# concat and archive fastq_pass
echo "# concatenating all reads to a single fastq.gz file"
cat $(find ${runf}/fastq_pass -name *.fastq | sort -k 1V,1) \
	| bgzip -c > ${tmpf}/${pfx}_pass.fq.gz

# copy Nanopore_SumStatQC files if present
if [ -d "Nanopore_SumStatQC" ]; then
	cp Nanopore_SumStatQC/config.yaml ${tmpf}/
	cp Nanopore_SumStatQC/Nanopore_SumStatQC.{pdf,html} ${tmpf}/
fi

echo "# data in ${tmpf} is ready for copy to the T: drive"

# check that the transfer drive is mounted
if grep -qs '/mnt/nuc-transfer ' /proc/mounts; then
    echo "# nuc-transfer (T:) is mounted. You can run the last command to copy the data to T:"
else
    echo "# nuc-transfer (T:) is NOT mounted, mount it first before issuing the next command."
fi

echo "# rsync -avr ${tmpf}/* /mnt/nuc-transfer/0003_Runs/GridION/${runf:0:4}-${runf:4:2}/${tmpf%%_data}"

echo "# after copy, please delete the local copy with:"

echo "# rm -rf ${tmpf}"
