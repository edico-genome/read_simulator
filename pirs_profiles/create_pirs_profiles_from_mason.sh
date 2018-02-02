#!/usr/bin/env bash

base=/opt/pirs-2.0.1
fasta=/mnt/vault/reference_genomes/Hsapiens/hg19/seq/hg19.fa
bam=/mnt/archive/TEST_RUNNER/mason_PE_2x250/33210/2/mason_PE_2x250.bam
outdir=`pwd`/mason_250/mason

echo "outdir: "$outdir

small_bam=$bam
# samtools view -h $bam | head -2000000 | samtools view -h  -b /dev/stdin > $small_bam

##################################################
# echo "Base-Calling Profile"
echo "$base/baseCalling_Matrix_calculator -b -r $fasta -i $bam -o $outdir"
$base/baseCalling_Matrix_calculator -b -r $fasta -i $small_bam -o $outdir 

if [[ ! $? ]]; then
    echo "failed"
    exit 1
fi

##################################################
echo "InDel Profile"
echo "BAM: $bam"
cmd="$base/indelstat_sam_bam $small_bam $outdir"
$cmd

if [[ ! $? ]]; then
    echo "failed"
    exit 1
fi


