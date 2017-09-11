
CONTIG_1="chr6_GL000250v2_alt"
CONTIG_2="chr6_GL000251v2_alt"

HG38_FA="/mnt/vault/reference_genomes/Hsapiens/hg38/seq/hg38.fa"

GET_FASTA="/home/cooper/products/UCSC_tools/faSomeRecords"

FDIR="/mnt/archive/jcr/projects/mapalign/alt_aware/sim/hla/chr6_GL000250v2_alt--chr6_GL000251v2_alt"

FASTA_1="$FDIR/${CONTIG_1}.fa"
FASTA_2="$FDIR/${CONTIG_2}.fa"

$GET_FASTA $HG38_FA <(echo "$CONTIG_1") $FASTA_1
$GET_FASTA $HG38_FA <(echo "$CONTIG_2") $FASTA_2

OUT_PRE="pirs.GL000250_GL000251"

pirs simulate -l 100 -x 35 --insert-len-mean=400 --insert-len-sd=40 --diploid \
     --base-calling-profile=/opt/pirs-2.0.1/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz \
     --indel-error-profile=/opt/pirs-2.0.1/Profiles/InDel_Profiles/phixv2.InDel.matrix \
     --gc-bias-profile=/opt/pirs-2.0.1/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat \
     --phred-offset=33 --no-indel-errors --no-gc-bias -c "gzip" -t 24 -o $OUT_PRE $FASTA_1 $FASTA_2

