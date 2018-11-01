[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![ngs-tools](ngstools.png) - Nanopore-Tools
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl and cran for R will help you add them)*

Please refer to the accompanying **[wiki](https://github.com/Nucleomics-VIB/nanopore-tools/wiki)** for examples and workflows.

### **spike_filter.sh**

**[spike_filter.sh](spike_filter.sh)** aligns all reads from a MinION dataset to a reference genome [eg lambda:](https://www.ncbi.nlm.nih.gov/nuccore/J02459.1?report=fasta) with **[graphmap](https://github.com/isovic/graphmap)** to identify spiked sequences and creates a new fastq.gz file using **samtools** and either mapped or unmapped reads, that can be assembled using **canu** or used elsewhere.

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

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
