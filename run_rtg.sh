#!/bin/bash

if [[ -d /staging/tmp/rtg ]]; then
    rm -r /staging/tmp/rtg
fi
mkdir /staging/tmp/rtg

dragen_vcf=/mnt/suite_logs/28058/suite_def_tmp/mrjd_pirs/EId2/prefilter/5011-dragen_map_align_sort_dedup_C_bed_save_bam-DB_vlrd_chr1_normalNoise_varRate0.002.INV28058-EId2_vlrd.vcf

dragen_vcf=/mnt/suite_logs/28066/suite_def_tmp/mrjd_pirs/EId2/prefilter/5011-dragen_map_align_sort_dedup_C_bed_save_bam-DB_vlrd_chr1_normalNoise_varRate0.002.INV28066-EId2_vlrd.vcf

truth_vcf=/mnt/archive/sim_data/pirs_simulator/vlrd_chr1_normalNoise_varRate0.002_indelErrors1_errorRate1_varRate0.002/truth.vcf

cmd="/home/theoh/p4/sw/trunk/test/bin/Dragen_tools_package/dragen_tools" \
cmd+=" --tool-name VCF-COMPARE" \
cmd+=" --algorithm RTG" \
cmd+=" --tmp-work-dir /staging/tmp/rtg" \
cmd+=" --output-mode annotate" \
cmd+=" --vcf-score-field QUAL " \
cmd+=" --enable-gzip-output false -a true" \
cmd+=" --output-dir /staging/tmp/rtg" \
cmd+=" --dataset-name vlrd_chr1_normalNoise_varRate0.002 -iid 28058 -eid 2 " \
cmd+=" --target-bed /mnt/vault/theoh/mrjd/mrjd_test_chr1.bed " \
cmd+=" --dragen-vcf $dragen_vcf" \
cmd+=" --truth-vcf $truth_vcf" \
cmd+=" --reference-fasta /mnt/vault/reference_genomes/Hsapiens/hg19_onlyChr1/seq/hg19_chr1.fa"

echo $cmd
$cmd
