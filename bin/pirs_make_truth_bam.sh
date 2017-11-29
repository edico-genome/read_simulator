#!/bin/bash

# This script will convert a pirs read_info.gz file to a BAM-like truth alignment file.
# The BAM-like file will not have valid CIGAR or SEQ or QUAL fields.

# $1 is a valid SAM header txt file corresponding to the original FASTA used as input to
#    bcftools consensus..
# $2 is a read_info.gz file produced by pirs.
# $3 is the basename (full path) of the liftover file(s) produced by bcftools consensus.
# $4 is the output truth BAM file (unsorted).  To make the sorted truth BAM and its index,
# do the following:
#    sambamba sort --tmpdir=/staging/tmp -p -t 8 /path/to/unsorted.truth.bam

SCRIPT_DIR=`dirname $0`
XFORM_READ_INFO="${SCRIPT_DIR}/xform_read_info.pl"

cat $1 <($XFORM_READ_INFO <(zcat $2) $3)  | sambamba view -S -f bam /dev/stdin > $4
