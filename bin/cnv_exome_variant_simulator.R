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
  df = df[-idx,]
  return(df)
}

get_sub_dataframe = function(df, names, sample_nr){
		
  df = df[df["variantsubtype"] == names,] 
  msg0 = paste(names, collapse = ",")
  msg = paste("Number of samples found in GDB: ", msg0," ", length(df$chr))	
  print(msg)
  
  if ( length(df$chr) < sample_nr ){
     msg = paste0("Samples found in GDB: ",length(df$chr)," is less than required: ",sample_nr)
     print(msg)
     stop("abort")
  }
  
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
nrDeletions = as.numeric(as.character(opt$nrDeletions))
nrDuplications = as.numeric(as.character(opt$nrDuplications))

minEventLength = 100
maxEventLength = 100000


##################################################
# load data
data(weightsMechanisms, package="RSVSim")
data(weightsRepeats, package="RSVSim")

# arg parser
if (! length(target_chrs) == 1){
  stop("Invalid target chromosome selection")
}


##################################################
# parse the csv file
print("parse the csv file")
cc = rep("NULL", 20)
cc[c(1,2,3,4,5,6)] = NA
df = read.csv(file = DGV, sep = "", colClasses=cc, stringsAsFactors = FALSE) 

print("get varians in target chromosome")
target_chrs_stripped = gsub("chr", "", target_chrs)
if (! target_chrs_stripped %in% as.character(seq(1,22))){
  stop("Invalid target chromosome selection")
}
df = df[df["chr"] == target_chrs_stripped,] 
head(df)

print("convert ranges to numbers")
df$start = as.numeric(as.character(df$start))
df$end = as.numeric(as.character(df$end))
df$width = df$end - df$start

print("filter on min and max event lengths")
df = df[df$width < maxEventLength,]
df = df[df$width > minEventLength,]
head(df)

print("get exome target bed")
exome_bed_df = read.table(target_bed)

# be consistent between target bed and DGV
chr_is_in_string = grepl("chr", exome_bed_df[1,1], fixed=TRUE)
if (chr_is_in_string){
 df$chr <- paste("chr", df$chr, sep="")
}

print("convert DGV to bed like format")
df_ranges = df[, c('chr', 'start', 'end')]
df_ranges = with(df_ranges, GRanges(chr, IRanges(start+1, end)))

print("parse target bed as GRanges object")
colnames(exome_bed_df) = c('chr', 'start', 'end')
exons = with(exome_bed_df, GRanges(chr, IRanges(start+1, end)))
head(exons)

print("get cnvs that overlap at least one target region")
df = df[countOverlaps(df_ranges, exons) >= 1,]
head(df)

print("now remove any cnv's that overlap each other")
df = deleteOverlaping(df)
head(df)


dels = get_sub_dataframe(
     df = df, names = c("deletion", "loss"), sample_nr = nrDeletions) 
ins = get_sub_dataframe(
    df = df, names = c("insertion", "gain", "tandem"), sample_nr = nrDuplications)

dels_gr = getGRanges(df=dels, target_chrs=target_chrs)
tandems_gr = getGRanges(ins, target_chrs=target_chrs)
# ins_gr = getGRanges(df=ins, target_chrs=target_chrs)

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

