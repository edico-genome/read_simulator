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
  make_option(c("--cnv_db"), action="store", 
    	help="optionally specify positions where events occure (bed file format)",
	default="NULL"),
  make_option(c("--target_chrs"), action="store",
	help="optional chr to narrow the simulation", default="NULL")
)

opt = parse_args(OptionParser(option_list=option_list))
required_args_mode1 = c("outdir", "size_ins", "size_del", 
		    "size_dup", "nr_ins", "nr_dels", "nr_dups", "target_chrs")
required_args_mode2 = c("outdir", "cnv_db", "target_chrs")

modes = c(required_args_mode1, required_args_mode2)

# check mode
counter = 0
for (required_args in modes){
    counter = counter + 1
    for (n in required_args_mode1){
    	if ( opt[n] == "NULL" ) {
    	   msg = paste("parameter:", n, 
	       "must be provided. See script usage (--help)")
    	   stop(msg)
  	   }
     mode = counter
}

# output directory, should already exist
outdir = opt$outdir
target_chrs = opt$target_chrs


if ( mode == 1 ) {
   size_ins = as.numeric(opt$size_ins)
   size_del = as.numeric(opt$size_del)	
   size_dup = as.numeric(opt$size_dup)
	 
   nr_ins = as.numeric(opt$nr_ins)
   nr_del = as.numeric(opt$nr_del)
   nr_dup = as.numeric(opt$nr_dup)
} elseif ( mode == 2 ) {
   cnv_db = opt$cnv_db
} 


##################################################
# run cmds

# percCopiedIns=1 ( running with default = 0 )

# "chrB" and "startB"

# MODE 1
if ( mode == 1 ){
  sim = simulateSV(
    output=outdir, chrs=c(target_chrs), ins=nr_ins, sizeIns=size_ins, 
    dels=nr_dels, sizeDels=size_del, dups=nr_dups, sizeDups=size_dup,
    maxDups=10, verbose=TRUE, repeatBias=TRUE)

} else {

# MODE 2
  cnv_db = import(opt$cnv_db, format="bed") 
  cnv_db = cnv_db[seqnames(cnv_db) == opt$target_chrs]

  dups = cnv_db
  dups_sizes = width(ranges(dups))

  # percCopiedIns=1

  sim = simulateSV(	
      output=outdir, chrs=c(target_chrs), 
      ins=nr_ins,   sizeIns=size_ins,  regionsIns=ins, 
      dels=nr_dels, sizeDels=size_del, regionsDels=dels,
      dups=nr_dups, sizeDups=size_dup, regionsDups=dups, maxDups=10,
      verbose=TRUE, repeatBias=FALSE,  random=FALSE, bpSeqSize=50)
}


##################################################
# write output

foo <- metadata(sim)
insert.csv <- paste(outdir,"insertions.csv", sep="/")
write.table(foo$insertions, file=insert.csv, quote=FALSE, sep="\t")
del.csv <- paste(outdir,"deletions.csv", sep="/")
write.table(foo$deletions, file=del.csv, quote=FALSE, sep="\t")
dups.csv <- paste(outdir,"tandemDuplications.csv", sep="/")
write.table(foo$tandemDuplications, file=dups.csv, quote=FALSE, sep="\t")


