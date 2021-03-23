#!/bin/bash

# create reference gene BED for adaptive sampling

# hg38 from ensembl
reflink="http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
gfflink="http://ftp.ensembl.org/pub/release-103/gff3/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gff3.gz"
gtflink="http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz"

for link in ${reflink} ${gfflink} ${gtflink}; do
wget ${link}
file=$(basename ${link})
echo "# unzipping ${file}"
gunzip ${file}
done

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
