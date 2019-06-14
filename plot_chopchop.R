#!/usr/bin/env Rscript

# script: plot_chopchop.R
# create plots from chopchop predicted sgRNA 
# (run_chopchop_single.sh or run_chopchop_tiled.sh)
# SP 2019-06-14 v1

# required libraries
suppressMessages(library("readr"))
suppressMessages(library("Sushi"))

# read user input
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("# please provide the path to filtered_results.bed\n", call.=FALSE)
} else if (length(args)==1) {
  infile <- args[1]
}

# load data in
filtered_guides <- suppressMessages(read_delim(infile,
                                               "\t", 
                                               escape_double = FALSE, 
                                               col_names = FALSE,
                                               trim_ws = TRUE,
                                               skip = 2))

colnames(filtered_guides) <- c("chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")

# convert strand to numeric
filtered_guides$strand[filtered_guides$strand=="+"] <- 1
filtered_guides$strand[filtered_guides$strand=="-"] <- -1
filtered_guides$strand <- as.numeric(filtered_guides$strand)

# set window
chrom <- unique(filtered_guides$chrom)
chromstart <- min(filtered_guides$start)
chromend <- max(filtered_guides$end)

# extract min data for Sushi
dat <- as.data.frame(filtered_guides[,1:6])

# plot and save to file
picture <- paste0(infile, "_plot.png", sep="")

png(picture,
    width = 1280, 
    height = 640, 
    units = "px", 
    pointsize = 12,
    bg = "white")

plotBed(beddata = dat, 
        chrom = chrom, 
        chromstart = chromstart,
        chromend = chromend,
        colorby = dat$strand,
        colorbycol = SushiColors(2),
        row = "auto",
        wiggle=0.001,
        splitstrand=TRUE)

labelgenome(chrom, chromstart, chromend, n=10, scale="Mb")

legend("topright", inset=0, legend=c("reverse","forward"),
       fill=SushiColors(2)(2),
       border=SushiColors(2)(2),
       text.font=2,
       cex=0.75)

null <- dev.off()
