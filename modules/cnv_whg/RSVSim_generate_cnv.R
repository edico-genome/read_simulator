#!/usr/bin/env Rscript

library(RSVSim)
library(GenomicRanges)
library(optparse)


##################################################
# load data
data(weightsMechanisms, package="RSVSim")
data(weightsRepeats, package="RSVSim")


##################################################
# cmd args
print("Parse options")

option_list = list(
  make_option(c("--outdir"), action="store", default="NULL", help=""),
  make_option(c("--size_ins"), action="store", help="", default="NULL"),
  make_option(c("--size_del"), action="store", help="", default="NULL"),
  make_option(c("--size_dup"), action="store", help="", default="NULL"),
  make_option(c("--nr_ins"), action="store", help="", default="NULL"),
  make_option(c("--nr_dels"), action="store", help="", default="NULL"),
  make_option(c("--nr_dups"), action="store", help="", default="NULL"),
  make_option(c("--target_chrs"), action="store",
	help="optional chr to narrow the simulation", default="NULL")
)


opt = parse_args(OptionParser(option_list=option_list))
required_args_mode1 = c("outdir", "size_ins", "size_del", 
		    "size_dup", "nr_ins", "nr_dels", "nr_dups", "target_chrs")

modes = c(required_args_mode1)

# check mode
print ("Validate options")
for (required_args in modes)
   {
    for (n in required_args_mode1)
        {
    	if ( opt[n] == "NULL" ) {
    	   msg = paste("parameter:", n, 
	       "must be provided. See script usage (--help)")
    	   stop(msg)
  	   }
    }
}

# output directory, should already exist
outdir = opt$outdir
target_chrs = opt$target_chrs

size_ins = as.numeric(opt$size_ins)
size_del = as.numeric(opt$size_del)
size_dup = as.numeric(opt$size_dup)
 
nr_ins = as.numeric(opt$nr_ins)
nr_dels = as.numeric(opt$nr_dels)
nr_dups = as.numeric(opt$nr_dups)


##################################################
# run cmds

# MODE 1
print("Run CMD")
sim = simulateSV(
        seed=246,
        output=outdir, chrs=c(target_chrs), ins=nr_ins, sizeIns=size_ins, 
    	dels=nr_dels, sizeDels=size_del, dups=nr_dups, sizeDups=size_dup,
        maxDups=10, verbose=TRUE, repeatBias=FALSE, percCopiedIns=1)


##################################################
# write output

foo <- metadata(sim)
insert.csv <- paste(outdir,"insertions.csv", sep="/")
write.table(foo$insertions, file=insert.csv, quote=FALSE, sep="\t")
del.csv <- paste(outdir,"deletions.csv", sep="/")
write.table(foo$deletions, file=del.csv, quote=FALSE, sep="\t")
dups.csv <- paste(outdir,"tandemDuplications.csv", sep="/")
write.table(foo$tandemDuplications, file=dups.csv, quote=FALSE, sep="\t")


