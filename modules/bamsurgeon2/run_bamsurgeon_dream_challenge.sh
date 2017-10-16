#!/usr/bin/env bash

# Constants
NUM_PROCESSES=4
BASEDIR=/home/gavinp/Downloads/new_bamsurgeon/bamsurgeon/bams_p27/bin/
FASTA=/mnt/vault/reference_genomes/Hsapiens/hs37d5/seq/hs37d5.fa
FASTA=/home/gavinp/Downloads/new_bamsurgeon/bamsurgeon/hs37d5.fa

# activate python venv
cmd="source /home/theoh/git_master/read_simulator/venv/bin/activate"
# cmd="source /home/gavinp/Downloads/new_bamsurgeon/bamsurgeon/bams_p27/bin/activate.sh"
echo "Activate python venv"
echo $cmd
$cmd

# filter secondary alignments from BAM
# samtools view -F 3840 PrecisionFDA_HG001.bam -b > $NOSECONDARY
# samtools view -F 3840 filters out secondary/supplementary alignments
# samtools view -F 3840 PrecisionFDA_HG001.bam -b > $NOSECONDARY
# NOSECONDARY=/mnt/archive/gavinp/PrecisionFDA_filtered.bam


INBAM=/mnt/archive/edico_s/results/cancer_data/Dream_challenge_set4/hs37d5/Dream_challenge_set4_normal_dragen_rhdr.bam
NOSECONDARY=/staging/bamsurgeon/dream_challenge_bamsurgeon_no_secondary_aligns.bam
cmd="samtools view -F 3840 $INBAM -b > $NOSECONDARY"
echo $cmd

if [[ ! -d /staging/bamsurgeon/ ]]; then 
    mkdir /staging/bamsurgeon/
fi



if [[ ! -f $NOSECONDARY ]]; then 
    samtools view -F 3840 $INBAM -b > $NOSECONDARY
fi 

echo "INDEX BAM"
echo "samtools index $NOSECONDARY"
# samtools index $NOSECONDARY

NEW_BAM=/staging/bamsurgeon/dream_challenge_bamsurgeon_tumor.bam

# input variants
# test_dir=/home/gavinp/Downloads/new_bamsurgeon/bamsurgeon/test_data
# var_in=$test_dir/random_snvs.txt
var_in="/mnt/archive/sim_data/dSim/bamsurgeon/Dream_challenge_bamsurgeon/dream_challenge_bamsurgeon.bed"

echo RUN BAMSURGEON
#https://github.com/adamewing/bamsurgeon/issues/87
cmd="$BASEDIR/addsnv.py -f $NOSECONDARY -r $FASTA -o $NEW_BAM -v $var_in -p $NUM_PROCESSES"
echo $cmd
$cmd

echo "move to nas"
cmd="cp $NEW_BAM /mnt/archive/sim_data/dSim/bamsurgeon/Dream_challenge_bamsurgeon/dream_challenge_bamsurgeon.bam"
echo $cmd
$cmd