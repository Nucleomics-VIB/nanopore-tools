#!/usr/bin/env Rscript

# script: Pinfish_QC-plots.R
# create plots from Pinfish cluster_memberships.tsv
# adapted from Isoseq3_QC-plots.R
# SP 2019-05-31 v1

# required libraries
suppressMessages(library("readr"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("scales"))
suppressMessages(library("grid"))
suppressMessages(library("gridExtra"))
suppressMessages(library("data.table"))

# read user input
args <- commandArgs(trailingOnly=TRUE)

## test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("# please provide the path to cluster_memberships.tsv\n", call.=FALSE)
} else if (length(args)==1) {
  infile <- args[1]
} else if (length(args)==2) {
  infile <- args[1]
  runpfx <- args[2]
} 

####################################################################################
# load cluster_memberships.tsv
####################################################################################

polished_cluster_report <- suppressMessages(read_tsv(infile))
colnames(polished_cluster_report) <- c("read_id", "cluster_id")

# set title
runID <- ifelse(exists("runpfx"), runpfx, "NA")
title <- paste0("Pinfish QC plots for run ", runID , sep="")

# plot Read support per transcript until 'lim' to avoid long tail
hist.data <- polished_cluster_report %>%
  group_by(cluster_id) %>% 
  tally() %>%
  arrange(desc(n)) %>%
  mutate(cum_sum = cumsum(n)) %>%
  mutate(percent = cum_sum/sum(n))

# add table for Read Nvalues (similar to N50 but whole range of N's)
# X% of the transcripts have Y or more supporting Reads
cum.data <- data.frame(percent=numeric(), Read_count=numeric())
for (lim in seq(0, 1, by=0.1)) {
  dat <- suppressWarnings(hist.data[min(which(hist.data$percent>lim)),])
  cum.data <- rbind(cum.data, c(100*(1-lim), dat$n))
}
colnames(cum.data) <- c("percent", "min_read_count")

lim <- 75
p1 <- suppressMessages(ggplot(hist.data[hist.data$n<=lim,], aes(n)) +
  geom_bar() +
  scale_x_continuous(breaks=c(1:lim)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45),
        axis.text.y=element_text(size=rel(1)),
        text = element_text(size=12)) +
  xlab("number of ONT reads in cluster") +
  ylab("Transcript count"))

# 10x subset 0 to 100% & count unique Transcripts
# remove transcripts with only one row & count unique Transcripts
plot.data <- data.frame(Read_sample=numeric(), Transcript_count=numeric(), min=factor())

# set limits
cluster.cutoffs <- c(5, 20)

for (iter in seq(1, 10, by=1)) {
for (i in seq(0, 1, by=0.1)) {
  # sample and count
  dat <- sample_frac(polished_cluster_report, i)
  plot.data <- rbind(plot.data, c( i, length(unique(dat$cluster_id)), 0) )
  # group and add counts
  dat <- dat %>%
    group_by(cluster_id) %>%
    summarise(n = n()) 
  # filter at FLNC=cluster.cutoffs[[1]] and count
  res <- dat %>%
    filter(n == cluster.cutoffs[[1]])
  plot.data <- rbind(plot.data, c( i, length(unique(res$cluster_id)), 1) )
  # filter at FLNC>cluster.cutoffs[[1]] and count
  res <- dat %>%
    filter(n >= cluster.cutoffs[[1]])
  plot.data <- rbind(plot.data, c( i, length(unique(res$cluster_id)), 2) )
  # filter at FLNC>cluster.cutoffs[[2]] and count
  res <- dat %>%
    filter(n >= cluster.cutoffs[[2]])
  plot.data <- rbind(plot.data, c( i, length(unique(res$cluster_id)), 3) )
  }
}

# merge the two dataframes and plot
colnames(plot.data) <- c("Read_sample", "Transcript_count", "min")

p2 <- suppressMessages(ggplot(plot.data, aes(x=Read_sample, y=Transcript_count, group=min, shape==factor(min), colour=factor(min))) +
  geom_smooth(method="loess", se=FALSE, fullrange=TRUE, size=1) +
  geom_point(size=3) +
  scale_shape_identity() +
  scale_x_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x=element_text(size=rel(1)),
        axis.text.y=element_text(size=rel(1)),
        text = element_text(size=12)) +
  xlab("Read sample from cluster_memberships.tsv") +
  ylab("Transcript count (10 random pulls)") +
  scale_color_manual(labels = c("one or more", 
                                paste("exactly ", cluster.cutoffs[[1]]), 
                                                     paste("more than ", cluster.cutoffs[[1]], sep=""), 
                                                     paste("more than ", cluster.cutoffs[[2]], sep="")), 
                       values = hue_pal()(4)) +
  guides(color=guide_legend("Read count \n/Transcript"))
)
p2
pdf("Pinfish_QC-plots.pdf", onefile = TRUE, width=12, height=8)
lay <- rbind(c(1,1,1,1,4,2,2,2),
             c(3,3,3,3,3,4,4,4))
table.legend <- "X% transcripts have Y or less supporting Reads's\n\n\n"
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.65)),
  colhead = list(fg_params=list(cex = 0.65)),
  rowhead = list(fg_params=list(cex = 0.65)))
grid.arrange(p1, 
             tableGrob(cum.data, 
                       rows=NULL, 
                       theme = mytheme), 
             p2, 
             ncol=2,
             top=textGrob(title, gp=gpar(fontsize=20, font=3)),
             layout_matrix = lay,
             vp=viewport(width=0.95, height=0.95))
# add legend
grid.text(table.legend, 
          x=unit(0.8, "npc"), 
          y=unit(0.87, "npc"),
          gp=gpar(fontsize=11, font=3))
          
#          family="Times New Roman"))
null <- dev.off()
