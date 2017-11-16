#!/usr/bin/env bash

cnv_db=/home/gavinp/Downloads/GRCh37_hg19_variants_2016-05-15.allCNVs_in_20120518_targets.bed
target_bed=/mnt/archive/gavinp/1000_genomes/20120518.consensus_add50bp.chrom.bed

# MODE 1
cmd="./RSVSim_generate_modified_genome.R"
cmd+=" --outdir /staging/tmp"
cmd+=" --size_ins"
cmd+=" --size_del"
cmd+=" --size_dup"
cmd+=" --nr_ins"
cmd+=" --nr_dels"
cmd+=" --nr_dups"
cmd+=" --target_chrs chr20"

echo $cmd
$cmd


# MODE 2
cmd="./RSVSim_generate_modified_genome.R"
cmd+=" --outdir /staging/tmp"
cmd+=" --target_chrs chr20"
cmd+=" --cnv_db"

echo $cmd
$cmd
