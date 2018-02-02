#!/usr/bin/env bash

base=/opt/pirs-2.0.1

fasta=/mnt/vault/reference_genomes/Hsapiens/hs37d5/seq/hs37d5.fa
bam=/mnt/archive/TEST_RUNNER/SRR82646_hs37d5/28650/11//SRR82646_hs37d5.bam
vcf=/mnt/vault/NA12878_NIST_gold_set/refdata/giab/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
outdir=`pwd`/250

echo "outdir: "$outdir

# bam=/mnt/archive/sanket_experiments/GATK_References/00012_20160904_severinec/01_dragen/file_SRA056922_heads_deterministic/SRA056922_heads_deterministic_dragen.bam
# vcf=/mnt/vault/NA12878_NIST_gold_set/refdata/giab/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer_hg19.vcf.gz
# fasta=/mnt/vault/reference_genomes/Hsapiens/hg19/seq/hg19.fa

echo "create small vcf"
small_vcf=/staging/tmp/small.vcf
zcat $vcf > $small_vcf

echo "create small bam"
small_bam=$bam.small.bam
samtools view -h $bam | head -1000000 | samtools view -h  -b /dev/stdin > $small_bam
echo "SMALL SAM"
echo $small_bam



##################################################
# echo "Base-Calling Profile"
echo "$base/baseCalling_Matrix_calculator -b -r $fasta -i $small_bam -s $small_vcf -o $outdir"
$base/baseCalling_Matrix_calculator -b -r $fasta -i $small_bam -s $small_vcf -o $outdir 

if [[ ! $? ]]; then
    echo "failed"
    exit 1
fi

##################################################
echo "InDel Profile"
echo "BAM: $bam"
cmd="$base/indelstat_sam_bam $small_bam SRR250bp_indel_profile"
$cmd

if [[ ! $? ]]; then
    echo "failed"
    exit 1
fi



##################################################
# echo GC%-Depth Profile
# cmd="$base/soap.coverage -sam -cvg -onlyuniq -i ./soap/testA_2.soap ./soap/testA_2.single -refsingle ./ref/human.fa -o testA_soap.dresult -depthsingle testA_soap.depth > testA_soap.deplog 2>testA_soap.deperr"
# run "$cmd"
# ./gc_coverage_bias -r ./ref/human.fa -o testA_soap -w 100,150,200 testA_soap.depth



