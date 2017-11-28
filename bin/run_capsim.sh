#!/usr/bin/env bash

# exit the scipt if any statement returns a non-true value
set -e

#####################################################
# constants

BASE="/home/gavinp/.usr/local/bin/"
JAPSA_CAPSIM="$BASE/jsa.sim.capsim"

# list of processed probes files to select from
# "/home/gavinp/Downloads/sequencing_kits/S03723314/S03723314_Probes_chr20.txt.fa", \

#seed random generator
RANDOM=$$$(date +%s)

#####################################################
# Read the options this was called with
usage=$(cat <<EOF
usage: $(basename $0) 
  -n			          number of reads
  -l                              read length
  -f                              fragment length
  -o                              outdir
  -1                              modified fasta1
  -2                              modified fasta2
  -p                              probe_file
EOF
)

while getopts ':1:2:n:l:f:o:p:' opt ; do
    case "$opt" in
    1) mod_fasta1="${OPTARG}" ;;
    2) mod_fasta2="${OPTARG}" ;;	
    n) N_READS="${OPTARG}" ;;
    l) READ_LENGTH="${OPTARG}" ;;
    f) FRAG_LENGTH="${OPTARG}" ;;
    p) probe_file="${OPTARG}" ;;
    o) outdir="${OPTARG}"
       if [[ ! -d $outdir  ]]; then
	   echo "Outdir: $outdir does must exist"
	   exit 1
       fi ;;
    :) echo "$usage"
       exit 1 ;;    
    *) echo "$usage"
       exit 1 ;;
  esac
done


echo ""
echo "SETTINGS"
echo "=============================="
for key in mod_fasta1 mod_fasta2 N_READS READ_LENGTH FRAG_LENGTH outdir probe_file; do
    eval value=\${$key}
    if [[ -z "$value" ]]; then
	echo "Missing argument: $key"
	echo ""
	echo "$usage"
	exit 1
    else
	printf "%-20s:  %s\n" $key $value
    fi
done
echo "=============================="


#####################################################
function run {
    cmd=$1
    echo "Run cmd: $cmd"
    $cmd
    
    ret_val=$?
    if [[ $ret_val -ne 0 ]]; then
        echo "Cmd had a non-zero return value: $ret_val"
	exit 1
    fi
}

#####################################################
# MAIN

# avoid nextera - it is not working properly in the 
# simulator - seems to be a capsim issue, maybe related to length

echo "Processing $probe_file "
cmd="bowtie2-build $mod_fasta1 $outdir/bowtie_ref1"; run "$cmd"
cmd="bowtie2-build $mod_fasta2 $outdir/bowtie_ref2"; run "$cmd"

cmd="bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x $outdir/bowtie_ref1 -U $probe_file -S $outdir/probes.sam"; run "$cmd"
cmd="bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x $outdir/bowtie_ref2 -U $probe_file -S $outdir/probes.sam"; run "$cmd"

# cmd="samtools view -b $outdir/probes.sam | samtools sort -o $outdir/probes.bam"; run "$cmd"
cmd="samtools sort $outdir/probes.sam $outdir/probes_sorted"; run "$cmd"
cmd="samtools index $outdir/probes_sorted.bam"; run "$cmd"

capsim_options=" --probe $outdir/probes_sorted.bam --ID someid --fmedian $FRAG_LENGTH --miseq $outdir/output --illen $READ_LENGTH --num $N_READS "

# set java
ORIGINAL_PATH="$PATH"
echo "Set Java Version: 1.8.0"
export PATH="/usr/java/jdk1.8.0_72/jre/bin:$PATH"
echo "Export PATH : $PATH"


# run for 1st haplotype
cmd="$JAPSA_CAPSIM --reference $mod_fasta1 $capsim_options --logFile $outdir/capsim1.log "; run "$cmd"
# cmd="mv $outdir/output_1.fastq.gz $outdir/output_ref1_1.fastq.gz"; run "$cmd"
# cmd="mv $outdir/output_2.fastq.gz $outdir/output_ref1_2.fastq.gz"; run "$cmd"

# run for 2nd haplotype
# cmd="$JAPSA_CAPSIM --reference $mod_fasta2 $capsim_options --logFile $outdir/capsim2.log "; run "$cmd"
# cmd="mv $outdir/output_1.fastq.gz $outdir/output_ref2_1.fastq.gz"; run "$cmd"
# cmd="mv $outdir/output_2.fastq.gz $outdir/output_ref2_2.fastq.gz"; run "$cmd"

# reset java
echo "Reset Java to original version"
export PATH="$ORIGINAL_PATH"


# merge results
# fq1
# cmd="cat $outdir/output_ref1_1.fastq.gz $outdir/output_ref2_1.fastq.gz"
# echo "$cmd > $outdir/capsim_1.fastq.gz"
# $cmd > $outdir/capsim_1.fastq.gz

# fq2
# cmd="cat $outdir/output_ref1_2.fastq.gz $outdir/output_ref2_2.fastq.gz"
# echo "$cmd > $outdir/capsim_2.fastq.gz"
# $cmd > $outdir/capsim_2.fastq.gz

echo "done"






