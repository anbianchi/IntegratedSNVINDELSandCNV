#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(DNAcopy)
library(cn.mops)
library(GenomeInfoDb)


getwd()
setwd("./final/")


BAMFiles <- list.files(pattern=".bam$")

segments <- read.table("../100bp_exon.bed", sep="\t", as.is=TRUE)


gr <- GRanges(segments[,1], IRanges(segments[,2], segments[,3]))

#seqlevels(gr) <- gsub("chr","",seqlevels(gr))
X <- getSegmentReadCountsFromBAM(BAMFiles,GR=gr)


resCNMOPS <- exomecn.mops(X, segAlgorithm="DNAcopy", norm =FALSE)

resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
CNVs <- as.data.frame(cnvs(resCNMOPS))
write.table(CNVs, "../results/results_cnmops.bed", append = TRUE, sep = "\t",row.names=FALSE, quote = FALSE)