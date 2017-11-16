#!/bin/bash

ROOT=/home/theoh/p4/sw/trunk/test/third_party/gatk-protected-3.1.1/gatk-protected
GATK_EXE=$ROOT/target/GenomeAnalysisTK.jar
KEY=/mnt/bioinfotools/loni/tools/gatk_key/severine_edicogenome.com.key


vcf1=/mnt/vault/NA12878_NIST_gold_set/refdata/giab/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz 

vcf2=/mnt/vault/reference_genomes/Hsapiens/GRCh37/validation/dream-syn4/truth_small_variants.vcf.gz

ref=/mnt/vault/reference_genomes/Hsapiens/hs37d5/seq/hs37d5.fa

# cmd="java -jar GenomeAnalysisTK.jar 

cmd="java -Xms2g -Xmx5500m -jar $GATK_EXE "
cmd+=" -T CombineVariants -R $ref "
cmd+=" --variant:foo $vcf1 "
cmd+=" --variant:bar $vcf2 "
cmd+=" -o output.vcf "
cmd+=" -genotypeMergeOptions PRIORITIZE "
cmd+=" -priority foo,bar "

echo $cmd
$cmd

