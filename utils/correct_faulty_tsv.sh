#!/usr/bin/env bash

for dataset in cnv_whg_sim_mixed_100bp_100kbp_30X_chr1 \
    cnv_whg_sim_mixed_100bp_100kbp_20X_chr1 \
    cnv_whg_sim_mixed_100bp_100kbp_10X_chr1 \
    cnv_whg_sim_mixed_100bp_100kbp_5X_chr1; do 

    base="/mnt/archive/sim_data/dSim/$dataset/CNVgdbVCF/"
    echo $base
    vcf=$base/RSVSim_truth.vcf
    tsv_new=$base/RSVsim_truth_new.tsv
    tsv_new_sorted=$base/RSVsim_truth_sorted.tsv
    order_file=$base/order_file.txt

    # create truth tsv
    vcf_to_truth_table_script=/home/theoh/git/read_simulator/bin/rsvsim.vcf_to_cnv_truth_table.pl

    echo "$vcf_to_truth_table_script $vcf > $tsv_new"
    $vcf_to_truth_table_script $vcf > $tsv_new

    # sort 
    head -1 $tsv_new > $tsv_new_sorted 

    awk '{ print $1 "\t" $2 "\t" $2 "\t" $3 "\t" $4 "\t"$5 "\t" $6 }' $tsv_new | sortBed -i /dev/stdin -faidx $order_file | cut -f 1-2,4-7 >> $tsv_new_sorted

    echo $tsv_new_sorted

done