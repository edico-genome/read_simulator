#!/bin/bash

function usage {
   cat << EOF
Usage: cat rsvsim_table.tsv | unified_rsv_table.sh

           UNIFIED_RSV_TABLE.SH converts RSVSim structural variant tables to 
           a BED-like format with consistent fields, regardless of SV type.
EOF
   exit 1
}

#if [ $# -ne 2 ]; then usage; fi

#if [ "$1" != "-f" ]; then usage; fi

#  Header rows in stdin start with 'Name'.
#  First 4 fields of each non-header row of stdin must be 'NUMBER   NAME   CHROM   POS'.
#  We will sort records by CHROM and then POS, and use NAME field for vcf record ID. 

#TMPFILE=`mktemp /tmp/variants.XXXXXX`

# write header line
printf "#CHROM\tSTART\tEND\tID\tCHROM_2\tSTART_2\tEND_2\tLENGTH\tREPEATS\n"

# Convert the different RSVSim formats of insertions, deletions, and tandem dups to
# a consistent 9 column format. The REPEATS field is 1 for deletions and insertions,
# and is greater than 1 for tandem dups. REPEATS is not the same as ALLELE_COUNT!
awk '$1 != "Name"' | cut -f2- | awk -v OFS="\t" '{tmp1=$1;$1=$2;$2=$3;$3=$4;$4=tmp1}1' | awk -v OFS="\t" '{if ($4 ~ /^insertion_/) $9=1; if ($4 ~ /^tandemDuplication/) {$8=$5;$9=$6;$5=$6=$7="NA";} if ($4 ~ /^deletion/) {$8=$5;$9=1;$5=$6=$7="NA";} NF=9; print $0;}'

