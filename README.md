[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![ngs-tools](ngstools.png) - NGS-Tools
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl and cran for R will help you add them)*

Please refer to the accompanying **[wiki](https://github.com/Nucleomics-VIB/ngs-tools/wiki)** for examples and workflows.

### **fastastats.pl**

The perl script **[fastastats.pl](fastastats.pl)** compute simple length statistics on a multi-fasta file.
```bash
## Usage: fastastats.pl <-i fasta-file (can be compressed)>
# Additional optional filtering parameters are:
# <-m minsize in kb (0)>
# <-x maxsize in kb (1E+6 | 1Gb)>
# <-h to display this help>
```

### **fastaCleanHeader.pl**

The perl script **[fastaCleanHeader.pl](fastaCleanHeader.pl)** cleans complex fasta record names that may interfere with some applications.

```bash
## Usage: fastaCleanHeader.pl <-i fasta_file (required)>
# <-o output file name (default to "cleaned_"<infile>; optional)>
# <-c keep only the leftmost word (display_id field; optional)>
# <-d delimiter (default to '|'; optional)>
# <-z to compress results with bgzip>
# <-h to display this help>
```

### **run_freebayes.sh**

The bash script **[run_freebayes.sh](run_freebayes.sh)** calls variant with freebayes from mappings and a reference genome.
```bash
# Usage: run_freebayes.sh -i <bam file> -r <fasta reference>
# script version 1.0, 2016_09_28
# [optional: -m <minmapq|20>]
# [optional: -q <minbaseq|20>]
# [optional: -F <minaltfrac|0.01>]
# [optional: -C <minaltcnt|10>]
```

### **makeENSref.sh**

The bash script **[makeENSref.sh](makeENSref.sh)** creates a series of reference files from ENSembl FTP downloads.
```bash
# Usage: makeENSref.sh
# -o <organism (default to <homo_sapiens>)> 
# -b <build number (default to <GRCh38>)> 
# -r <release number (default to 88)>
# script version 1.0, 2017_04_05
# [-h for this help]
```

### **gepard_plot.sh**

The bash script **[gepard_plot.sh](gepard_plot.sh)** creates a xy-plot from two related assemblies. The Java GUI tools does the same but this script is applicable to multiple inputs in batch. The original tool can be found at **https://github.com/univieCUBE/gepard**.
```bash
# Usage: gepard_plot.sh -x <reference assembly> -y <draft assembly> -p <path to gepard.jar and matrices>
# script version 1.0, 2017_04_21
# [optional: -o <result path/prefix>]
# [optional: -w <word size:10>]
# [optional: -W <window size:0>]
# [optional: -J <java extra parameters (eg -Xmx1G, put between double quotes if it contains spaces)>
# [optional: -h <this help text>]
```

### **mauve_reorder.sh**

The bash script **[mauve_reorder.sh](mauve_reorder.sh)** reorders a draft-assembly based on a reference assembly at CLI. The Java GUI tools does the same but this script is applicable to multiple inputs in batch. The original tool can be found at **http://darlinglab.org/mauve/download.html**
```bash
# Usage: mauve_reorder.sh -i <draft assembly> -r <reference assembly> -p <mauve path>
# script version 1.0, 2017_04_21
# [optional: -o <result folder>]
# [optional: -m <available RAM|1G>]
# [optional: -h <this help text>]
```

### **mappability.sh**

The bash script **[mappability.sh](mappability.sh)** creates a mappability track for a given single-read length and a reference genome. Such tracks used to be present in IGV from the web server and have been removed. A mappability track allows you to identify those regions of the genome which as so ubiquitous that reads have no chance to map top them and explain why you may not observe coverage in these regions *[(Derrien et al., 2012)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377)*. For more info about this data type, please look at our page **http://wiki.bits.vib.be/index.php/Create_a_mappability_track**. The script depends on a number of resources including the **[gem-libraries](http://algorithms.cnag.cat/wiki/The_GEM_library)** and several of Jim Kent's Utilities **[KentUtils](https://github.com/ENCODE-DCC/kentUtils)**.

```bash
# Usage: mappability.sh
# -i <reference assembly fasta>)> 
# -l <single read length for prediction (default to 100bps)>
# -p <prefix for output data (default to input file prefix)>
# -t <number of threads (default to 1)> 
# script version 1.0, 2017_05_19
# [-h for this help]
```

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
