#!/bin/bash

# create reference gene BED for adaptive sampling

# Sc64-1-1 from ensembl
reflink="http://ftp.ensembl.org/pub/release-104/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
gfflink="http://ftp.ensembl.org/pub/release-104/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.104.gff3.gz"
gtflink="http://ftp.ensembl.org/pub/release-104/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.104.gtf.gz"

for link in ${reflink} ${gfflink} ${gtflink}; do
wget ${link}
file=$(basename ${link})
echo "# unzipping ${file}"
gunzip ${file}
done

# create fasta index
samtools faidx $(basename ${reflink} .gz)

# create a bed track for genes from gff3 data
# create bed track from gff data for genes
gawk 'BEGIN{FS="\t"; OFS="\t"}{
  if($3=="gene") {
    split($9,ann,";"); ensid=gensub(/ID=gene:(.*)/,"\\1","g",ann[1]); 
    hugo=gensub(/Name=(.*)/,"\\1","g",ann[2]); 
    print $1,$4,$5,ensid"|"hugo
    }
  }' *.gff3 \
  | sort -k 1V,1 -k 2n,2 \
  > genes.bed

# count per chromosome and save
cut -f 1 genes.bed | sort | uniq -c | sort -k 2V,2 | awk '{print $2, $1}' > genecounts.txt

# count gene size
gawk 'BEGIN{FS="\t"; OFS="\t";tot=0}{tot=tot+$3-$2}END{print tot}' genes.bed >> genecounts.txt
