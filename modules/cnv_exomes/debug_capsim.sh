#!/usr/bin/env bash

bed="/mnt/archive/gavinp/1000_genomes/20120518.consensus_add50bp.chrom.bed"

cmd="/home/theoh/git/read_simulator/modules/cnv_exomes/RSVSim_generate_modified_genome.R " 
cmd+=" --nHomozygousDeletions 10"
cmd+=" --nHeterozygousDeletions 10" 
cmd+=" --nTandemDuplications 20" 
cmd+=" --maxEventLength 40000" 
cmd+=" --minEventLength 1000" 
cmd+=" --outdir /mnt/archive/sim_data/dSim/cnv_exome/CNVExomeModFastas"
cmd+=" --fa_file /mnt/vault/reference_genomes/Hsapiens/hs37d5/seq/hs37d5.fa"
# cmd+=" --fa_file /mnt/archive/gavinp/1000_genomes/hs37d5.mod.fa"
cmd+=" --cnv_db /home/gavinp/Downloads/GRCh37_hg19_variants_2016-05-15.allCNVs_in_20120518_targets.bed"
cmd+=" --target_chrs 20"
cmd+=" --target_bed $bed"

echo $cmd
$cmd

diff /mnt/archive/gavinp/1000_genomes/hs37d5.mod.fa /mnt/vault/reference_genomes/Hsapiens/hs37d5/seq/hs37d5.fa