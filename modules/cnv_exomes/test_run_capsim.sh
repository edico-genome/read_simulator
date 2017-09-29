#!/usr/bin/env bash

fasta1=/mnt/archive/sim_data/dSim/cnv_exome/CNVExomeModFastas/genome_rearranged1.fasta
fasta2=/mnt/archive/sim_data/dSim/cnv_exome/CNVExomeModFastas/genome_rearranged2.fasta

cmd="./run_capsim.sh -n 800000 -l 300 -f 500 -o /staging/tmp -1 $fasta1 -2 $fasta2"
echo $cmd
$cmd
