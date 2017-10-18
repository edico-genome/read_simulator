#!/bin/bash

function usage {
   cat << EOF
Usage: cat *.csv | convert_RSVSim_to_VCF.sh -i /dev/stdin -f <fasta_file> -c <contig_order_file> -v <vcf_header_file>

           CONVERT_RSVSIM_TO_VCF.SH converts RSVSim structural variant tables to
           a VCF format and writes the result to STDOUT.

           bedtools is required to be in the user's PATH.
EOF
   exit 1
}

if [ $# -ne 8 ]; then usage; fi

TOOL_DIR=`dirname $0`

while [[ $# -gt 1 ]]; do

   key="$1"

   case $key in
      -i|--input-sv-table)
      SV_TABLE="$2"
      shift # past argument
      ;;
      -f|--fasta)
      FASTA="$2"
      shift # past argument
      ;;
      -c|--contig-order)
      ORDER_FILE="$2"
      shift # past argument
      ;;
      -v|--vcf-header)
      VCF_HDR="$2"
      shift # past argument
      ;;
      *)
           # unknown arg
           usage 
      ;;
  esac
  shift # past argument or value
done

cat $SV_TABLE | $TOOL_DIR/unified_rsv_table.sh | $TOOL_DIR/RSVSim_to_VCF.pl /dev/stdin $FASTA | sortBed -i /dev/stdin -faidx $ORDER_FILE | cat $VCF_HDR /dev/stdin

