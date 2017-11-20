#!/usr/bin/env Rscript

# rm(list=ls())
# cat("\f")

library(IRanges)
library(RSVSim)
library(GenomicRanges)
library(BSgenome)
library(Biostrings)
library(optparse)


##################################################
# functions
deleteOverlaping = function(df){
  end_prev_region = df$end[1]
  lines_to_discard = vector()   
  
  for (idx in 2:length(df$chr)){
    if (df$start[idx] <= end_prev_region){
      lines_to_discard = c(lines_to_discard, idx)
    } else {
      end_prev_region = df$end[idx]
    }
  }
  df = df[,-idx]
  return(df)
}

get_sub_dataframe = function(df, names, sample_nr){
  df = df[df["variantsubtype"] == names,] 
  s = sample(seq(1,length(df$chr)), size=sample_nr)
  s = sort(s)
  return(df[s,])
}

getGRanges = function(df, target_chrs)
{
  ranges = IRanges(as.numeric(df$start), as.numeric(df$end), names=NULL)
  seqnames = Rle(target_chrs, c(length(df$chr)))
  gr <- GRanges(seqnames = seqnames, ranges = ranges)
  return(gr)
}

##################################################
# cmd args
option_list = list(
  make_option(c("--nrDeletions"), action="store", default="NULL", help=""),
  make_option(c("--nrDuplications"), action="store", help="", default="NULL"),
  make_option(c("--outdir"), action="store", help="should already exist", default="NULL"),
  make_option(c("--fasta"), action="store", default="NULL"),
  make_option(c("--target_chrs"), action="store", help="required", default="NULL"),
  make_option(c("--target_bed"), action="store", help="required", default="NULL"),
  make_option(c("--DGV"), action="store", default="NULL")
)

required_args = c("nrDeletions", "nrDuplications",
                  "outdir", "fasta", "target_chrs", "target_bed", "DGV")

opt = parse_args(OptionParser(option_list=option_list))

# check opts are defined
for (n in required_args){
  if ( opt[n] == "NULL" ) {
    msg = paste("parameter:", n, "must be provided. See script usage (--help)")
    stop(msg)
  }
}

# parse args
outdir = opt$outdir
fasta = opt$fasta
genome = readDNAStringSet(fasta)
DGV = opt$DGV
target_chrs = opt$target_chrs
target_bed = opt$target_bed


minEventLength = 3000
maxEventLength = 100000


##################################################
# load data
data(weightsMechanisms, package="RSVSim")
data(weightsRepeats, package="RSVSim")

# arg parser
if (! length(target_chrs) == 1){
  stop("Invalid target chromosome selection")
}

target_chrs_stripped = gsub("chr", "", target_chrs)
if (! target_chrs_stripped %in% as.character(seq(1,22))){
  stop("Invalid target chromosome selection")
}


##################################################
# parse the csv file
print("parse the csv file")
cc = rep("NULL", 20)
cc[c(1,2,3,4,5,6)] = NA
df = read.csv(file = DGV, sep = "", colClasses=cc, stringsAsFactors = FALSE) 

print("get chromosme variants")
df = df[df["chr"] == target_chrs_stripped,] 

df$start = as.numeric(df$start)
df$end = as.numeric(df$end)
df$width = df$end - df$start
df = df[df$width < maxEventLength,]
df = df[df$width > minEventLength,]

print("convert DGV to bed like format")
df_ranges = df[, c('chr', 'start', 'end')]
df_ranges = with(df_ranges, GRanges(chr, IRanges(start+1, end)))

print("get exome target bed")
exome_bed_df = read.table(target_bed)
colnames(exome_bed_df) = c('chr', 'start', 'end')
exons = with(exome_bed_df, GRanges(chr, IRanges(start+1, end)))

# get cnvs that overlap at least one target region
df = df[countOverlaps(df_ranges, exons) > 1,]

# now remove any cnv's that overlap each other
df = deleteOverlaping(df)

dels = get_sub_dataframe(df = df, names = c("deletion", "loss"), sample_nr = 10) 
ins = get_sub_dataframe(df = df, names = c("insertion", "gain", "tandem"), sample_nr = 10)

dels_gr = getGRanges(df=dels, target_chrs=target_chrs)
tandems_gr = getGRanges(ins, target_chrs=target_chrs)
# ins_gr = getGRanges(df=ins, target_chrs=target_chrs)

# name the variables
names(dels_gr) <- paste('deletion', 1:length(dels_gr), sep='')
names(tandems_gr) <- paste('tandemDuplication', 1:length(tandems_gr), sep='')

##################################################
# run simulation
sim = simulateSV(output=outdir, chrs=c(target_chrs), genome=genome, 
                 regionsDels=dels_gr,
                 regionsDups=tandems_gr,
                 # regionsIns=ins_gr,	
		 maxDups=10,
                 verbose=TRUE, repeatBias=FALSE,  random=FALSE)


##################################################
# write output

foo <- metadata(sim)
insert.csv <- paste(outdir,"insertions.csv", sep="/")
write.table(foo$insertions, file=insert.csv, quote=FALSE, sep="\t")
del.csv <- paste(outdir,"deletions.csv", sep="/")
write.table(foo$deletions, file=del.csv, quote=FALSE, sep="\t")
dups.csv <- paste(outdir,"tandemDuplications.csv", sep="/")
write.table(foo$tandemDuplications, file=dups.csv, quote=FALSE, sep="\t")

