#!/usr/bin/env bash

cnv_db=/home/gavinp/Downloads/GRCh37_hg19_variants_2016-05-15.allCNVs_in_20120518_targets.bed

target_bed=/mnt/archive/gavinp/1000_genomes/20120518.consensus_add50bp.chrom.bed

cmd="./RSVSim_generate_modified_genome.R"
cmd+=" --nHomozygousDeletions 10"
cmd+=" --nHeterozygousDeletions 10"
cmd+=" --nTandemDuplications 20"
cmd+=" --maxEventLength 40000"
cmd+=" --minEventLength 1000"
cmd+=" --outdir /staging/tmp"
cmd+=" --fa_file /mnt/archive/gavinp/1000_genomes/hs37d5.mod.fa"
cmd+=" --cnv_db $cnv_db" 
cmd+=" --target_chrs 20"
cmd+=" --target_bed $target_bed"

echo $cmd
$cmd