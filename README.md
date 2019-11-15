[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![ngs-tools](ngstools.png) - Nanopore-Tools
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl and cran for R will help you add them)*

Please refer to the accompanying **[wiki](https://github.com/Nucleomics-VIB/nanopore-tools/wiki)** for examples and workflows.

### **spike_filter.sh**

**[spike_filter.sh](spike_filter.sh)** aligns all reads from a MinION dataset to a reference genome [eg lambda:](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1?report=fasta) (NC_001416.1) with **[graphmap](https://github.com/isovic/graphmap)** to identify spiked sequences and creates a new fastq.gz file using **samtools** and either mapped or unmapped reads, that can be assembled using **canu** or used elsewhere.

```bash
Usage: spike_filter.sh 
#   -i <nanopore_reads.fastq (required)>
#   -r <spiked reference (lambda or any other spiked genome, required)>
#   -t <graphmap alignment threshold (default to 1e-0)>
#   -e <read error rate (default to 0.45)>
#   -C <the spiked genome is circular (default OFF = linear)>
#   -s <keep only spiked reads instead (reverse-mode)>
#   -h <show this help>
```
**[spike_filter_minimap2.sh](spike_filter_minimap2.sh)** aligns all reads from a MinION dataset to a reference genome [eg lambda:](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1?report=fasta) (NC_001416.1) with **[minimap2](https://github.com/lh3/minimap2)** to identify spiked sequences and creates a new fastq.gz file using **samtools** and either mapped or unmapped reads, that can be assembled using **canu** or used elsewhere.

```bash
Usage: spike_filter_minimap2.sh
#   -i <nanopore_reads.fastq (required)>
#   -r <spiked reference (lambda or any other spiked genome, required)>
#   -t <threads to be used for alignment (default to 8)>
#   -x <path to minimap2>
#   -s <keep only spiked reads instead (reverse-mode)>
#   -h <show this help>
```

**[run_chopchop_single.sh](run_chopchop_single.sh)** design sgRNA around a region of interest (ROI) for ONT Cas9 enrichment. The value of -m is indicative of the ONT molecule average size and only used to stop the code when the requested capture is smaller than the actual DNA molecules.

```bash
# Usage: run_chopchop_single.sh
# -r <region (eg. chr1:16920000-16940000)>
# -g <reference genome used by ChopChop (default to hg38)>
# -w <prediction window (default to 3000bps either side)>
# -m <median DNA fragment size (default 30000)>
# script version 1.0, 2019_06_12
# [-h for this help]
```

**[run_chopchop_tiled.sh](run_chopchop_tiled.sh)** design sgRNA around a large region of interest (ROI) for ONT Cas9 enrichment. Predict in tiles of width -w kb to enrich a region larger than -s kb. The value of -m is indicative of the ONT molecule average size and only used to stop the code when the requested capture is smaller than the actual DNA molecules.

Run this script twice with a shift of -s/2 kb in order to select guides for two separate reactions. The two Cas9 reactions performed on separate DNA aliquotes can be mixed and sequenced together to lead to overlapping read pilups as shown below.

```
setA:
-------- -------- -------- --------
setB:
    -------- -------- -------
```

```bash
# Usage: run_chopchop_tiled.sh
# -r <region (eg. chr1:16740273-16972964)>
# -g <reference genome used by ChopChop (default to hg38)>
# -s <tile width (default to 20000)>
# -w <prediction window (default to 1500 either side)>
# -m <median DNA fragment size (default 30000)>
# script version 1.0, 2019_06_12
# [-h for this help]
```

**[Nanopore_SumStatQC.Rmd](Nanopore_SumStatQC.Rmd)** is a slightly modified copy of the original **[ont_tutorial_basicqc](https://github.com/nanoporetech/ont_tutorial_basicqc)** R markdown. The file depends on the accessory config.yaml where the user edits the metadata. The report is then produced as PDF or HTML using either

* R --slave -e 'rmarkdown::render("Nanopore_SumStatQC.Rmd", "pdf_document")'
* R --slave -e 'rmarkdown::render("Nanopore_SumStatQC.Rmd", "html_document")'

A significant number of R dependencies are required (see top of the code)

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
