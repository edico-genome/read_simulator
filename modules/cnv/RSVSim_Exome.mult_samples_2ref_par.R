
# github.com/fritzsedlazeck/SURVIVOR
library(RSVSim)
if(packageVersion("RSVSim") < "1.17.0") {
    stop("Need RSVSim 1.17.0")
}
library(GenomicRanges)
library(rtracklayer)
library(doParallel)
library(foreach)

# Parallelize generating the references 
numCores <- detectCores()
cl <- makeCluster(max(numCores-1, 1)) # use numCores-1
registerDoParallel(cl)

data(weightsMechanisms, package="RSVSim")
data(weightsRepeats, package="RSVSim")

# function removeOverlappingRegions
# this function repeatedly removes the region that overlaps the most other regions
# until there are no overlaps
removeOverlappingRegions <- function(bed){
    my_bed = bed
    maxVal = 10
    c=countOverlaps(my_bed, my_bed)
    removed <- 0
    while(maxVal>1)
    {
        maxIndex = which.max(c)
        #cat("maxIndex = ", maxIndex, "\n")
        my_bed=my_bed[-maxIndex]
        c=countOverlaps(my_bed, my_bed)
        maxVal = max(c)
        removed <- removed+1
    }
    cat("Removed ", removed, " overlapping regions, ", length(my_bed), "remaining regions \n")
    return(my_bed)
}

# specify count for each event
nHomozygousDeletions <- 10
nHeterozygousDeletions <- 10
nTandemDuplications <- 20

# limit CNV event length
maxEventLength = 40000
minEventLength = 1000
#https://bioconductor.statistik.tu-dortmund.de/packages/3.2/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf 
# Start the clock!
ptm <- proc.time()
inputs = 1:3    # which samp* to generate
# generate samples 1..
#for(i in 1:nSamplesToGenerate){ # was serial, now parallelized
cat("Parallelizing genome generation. There is no standard output from individual jobs. Please wait \n")
foreach(i=inputs,
        .packages=c('GenomicRanges', 'rtracklayer', 'RSVSim')) %dopar% {
    cat("Working on sample ", i, "\n")
    target_bed <- read.table("/mnt/archive/gavinp/1000_genomes/20120518.consensus_add50bp.chrom.bed")
    colnames(target_bed) <- c('chr', 'start', 'end')
    exons <- with(target_bed, GRanges(chr, IRanges(start+1, end)))
    outdir <- paste("/staging/gavinp/struct_vars/CNV/simul/RSVSim/exome/samp", i, sep="")
    fa_file <- "/mnt/archive/gavinp/1000_genomes/hs37d5.mod.fa"
    # load the cnv database, then filter for target chromosome & length, retain regions that overlap 
    # with target exons then remove largest-sized overlapping region until there are no overlaps
    bed = import("/home/gavinp/Downloads/GRCh37_hg19_variants_2016-05-15.allCNVs_in_20120518_targets.bed", format="bed")
    bed = bed[seqnames(bed)==20]    # filter for target chromosomes
    bed = bed[width(ranges(bed))<maxEventLength]
    bed = bed[width(ranges(bed))>minEventLength]
    bed = bed[countOverlaps(bed, exons) > 1]
    initial_bed = bed
    bed = removeOverlappingRegions(bed)
    # reset overlap flag
    foundOverlaps = TRUE
    if (length(bed)<(nHomozygousDeletions+nHeterozygousDeletions+nTandemDuplications)){
        stop("Non-overlapping CNV events with length < maxEventLength are too few to support specified number of deletions/duplications ")
    }
    # ensure homozygous deletions do not overlap 
    while(foundOverlaps){
        cat("Generating homozygous Deletions\n")
        homo_indices = sample(length(bed), size=nHomozygousDeletions, replace=FALSE)
        homoDels = sort(bed[homo_indices,])
        bed = bed[-homo_indices,]
        foundOverlaps=length(findOverlaps(homoDels, ignoreSelf=TRUE))>0
        cat("Found overlaps", foundOverlaps,"\n")
    }
    foundOverlaps = TRUE
    while(foundOverlaps){
        cat("Generating heterozygous Deletions\n")
        hetero_indices = sample(length(bed), size=nHeterozygousDeletions, replace=FALSE)
        heteroDels = sort(bed[hetero_indices,])
        bed = bed[-hetero_indices,]
        foundOverlaps=(sum(countOverlaps(homoDels, heteroDels)>0) || length(findOverlaps(heteroDels, ignoreSelf=TRUE))>0)
        cat("Found overlaps", foundOverlaps,"\n")
    }
    heteroDels1_indices = sample(length(heteroDels), size=floor(nHeterozygousDeletions/2), replace=FALSE)
    heteroDels1 = heteroDels[heteroDels1_indices]
    heteroDels2 = heteroDels[-heteroDels1_indices]
    foundOverlaps = TRUE
    while(foundOverlaps){
        cat("Generating tandem duplications\n")
        dup_indices = sample(length(bed), size=nTandemDuplications, replace=FALSE)
        dups = sort(bed[dup_indices,])
        bed = bed[-dup_indices,]
        # do we need dups1 and dups2?
        foundOverlaps=(sum(countOverlaps(homoDels, dups)>0) || sum(countOverlaps(heteroDels, dups)>0) || length(findOverlaps(dups, ignoreSelf=TRUE))>0 ) 
        cat("Found overlaps", foundOverlaps,"\n")
    }

    dups_sizes = width(ranges(dups))
    # generate TWO reference genomes so that we can simulate heterozygous deletions
    for (p in 1:2){
        if (p==1){
                gRangeDels <- c(homoDels, heteroDels1)
                gRangeDels <- sort(gRangeDels)
        } else {
                gRangeDels <- c(homoDels, heteroDels2)
                gRangeDels <- sort(gRangeDels)
        }
        names(gRangeDels) <- paste('deletion', 1:length(gRangeDels), sep='')
        names(dups) <- paste('tandemDuplication', 1:length(gRangeDels), sep='')
        #deletions_bed_out <- file.path(outdir, "deletions1.bed")
        #duplications_bed_out <- file.path(outdir, "duplications1.bed")
        #export.bed(gRangeDels, deletions_bed_out)
        #export.bed(gRangeDups, duplications_bed_out)
        # one thing to note is that the tandem duplications are random multiplicities for each reference, so it is likely the resulting copy numbers will not match between the genomes
        sim = simulateSV(output=outdir, genome=fa_file, chrs=c("20"), random=FALSE, dups=length(dups), regionsDels=gRangeDels, maxDups=4, sizeDups=dups_sizes, regionsDups=dups, verbose=TRUE, repeatBias=FALSE, bpSeqSize=50)
        
        foo <- metadata(sim)
        insert.csv <- paste(outdir,"/insertions", p, ".csv", sep="")
        write.table(foo$insertions, file=insert.csv, quote=FALSE, sep="\t")
        del.csv <- paste(outdir,"/deletions", p, ".csv", sep="")
        write.table(foo$deletions, file=del.csv, quote=FALSE, sep="\t")
        dups.csv <- paste(outdir,"/tandemDuplications", p, ".csv", sep="")
        write.table(foo$tandemDuplications, file=dups.csv, quote=FALSE, sep="\t")
        orig_name <- paste(c(outdir,"/genome_rearranged.fasta"), collapse="")
        ref_name <- paste(c(outdir,"/genome_rearranged",p,".fasta"), collapse="")
        sprintf("Renaming genome_rearranged.fasta to ", ref_name)
        file.rename(orig_name, ref_name)
    }

    sprintf("Wrote output to %s ", outdir)
}

cat("Total elapsed time ", proc.time() - ptm, "\n")
