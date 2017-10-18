#!/usr/bin/env Rscript
library(RSVSim)

library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)

print(args)

chrom <- args[1]
outdir <- args[2]
size_ins <- as.integer(args[3])
size_del <- as.integer(args[4])
size_dup <- as.integer(args[5])



print(chrom)
data(weightsMechanisms, package="RSVSim")
data(weightsRepeats, package="RSVSim")


sim = simulateSV(output=NA, chrs=c(chrom), ins=10, sizeIns=size_ins, percCopiedIns=1, dels=10, sizeDels=size_del, dups=10, sizeDups=size_dup, maxDups=4, verbose=TRUE, repeatBias=FALSE)



foo <- metadata(sim)
insert.csv <- paste(outdir,"insertions.csv", sep="/")
write.table(foo$insertions, file=insert.csv, quote=FALSE, sep="\t")
del.csv <- paste(outdir,"deletions.csv", sep="/")
write.table(foo$deletions, file=del.csv, quote=FALSE, sep="\t")
dups.csv <- paste(outdir,"tandemDuplications.csv", sep="/")
write.table(foo$tandemDuplications, file=dups.csv, quote=FALSE, sep="\t")


