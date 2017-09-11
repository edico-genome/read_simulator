
SIMDIR="/mnt/archive/jcr/projects/mapalign/alt_aware/sim/hla/chr6_GL000250v2_alt--chr6_GL000251v2_alt"
TOOLDIR="/home/cooper/p4/sw/trunk/test/bin/ma_qual_pipeline/alt_aware"

FASTA_to_VCF="$TOOLDIR/fasta_sam_to_vcf.pl"
XFORM_PIRS_ALIGNS="$TOOLDIR/xform_pirs_read_info.pl"

# create truth vcf for our ALT contig simulation
$FASTA_to_VCF /mnt/vault/reference_genomes/Hsapiens/hg38/seq/hg38.fa /opt/bwakit-0.7.12-0/resource-GRCh38/hs38DH.fa.alt chr6_GL000250v2_alt:1-4672374 chr6_GL000251v2_alt > chr6_GL000250v2_alt--chr6_GL000251v2_alt.vcf

# create sim reads
$SIMDIR/pirs.250--251.sh

# create truth sam
$XFORM_PIRS_ALIGNS <(zcat pirs_reads_100_400.read.info.gz) /opt/bwakit-0.7.12-0/resource-GRCh38/hs38DH.fa.alt > pirs_reads_100_400.truth.sam
