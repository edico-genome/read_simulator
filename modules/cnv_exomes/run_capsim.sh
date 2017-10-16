#!/usr/bin/env bash

# exit the scipt if any statement returns a non-true value
set -e

#####################################################
# constants

JAPSA_CAPSIM="/home/gavinp/.usr/local/bin/jsa.sim.capsim"

# list of processed probes files to select from
# "/home/gavinp/Downloads/sequencing_kits/S03723314/S03723314_Probes_chr20.txt.fa", \
PROBE_FILES=(
    "/home/gavinp/Downloads/sequencing_kits/SeqCapEZ/exome_probes.txt", \
    "/home/gavinp/Downloads/sequencing_kits/SeqCapEZ/results3.txt", \
    "/home/gavinp/Downloads/sequencing_kits/S07604624/results.txt")

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
EOF
)

while getopts ':1:2:n:l:f:o:' opt ; do
    case "$opt" in
    1) mod_fasta1="${OPTARG}" ;;
    2) mod_fasta2="${OPTARG}" ;;	
    n) N_READS="${OPTARG}" ;;
    l) READ_LENGTH="${OPTARG}" ;;
    f) FRAG_LENGTH="${OPTARG}" ;;
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
for key in mod_fasta1 mod_fasta2 N_READS READ_LENGTH FRAG_LENGTH outdir; do
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
    echo $cmd
    $cmd
}

#####################################################
# MAIN


# avoid nextera - it is not working properly in the 
# simulator - seems to be a capsim issue, maybe related to length
selected_probe=${PROBE_FILES[$RANDOM % ${#PROBE_FILES[@]} ]}

echo "Processing $selected_probe "
cmd="bowtie2-build $mod_fasta1 $outdir/bowtie_ref1"; run "$cmd"
cmd="bowtie2-build $mod_fasta2 $outdir/bowtie_ref2"; run "$cmd"

cmd="bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x $outdir/bowtie_ref1 -U $selected_probe -S $outdir/probes.sam"; run "$cmd"

cmd="bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x $outdir/bowtie_ref2 -U $selected_probe -S $outdir/probes.sam"; run "$cmd"

#bowtie2 --local --fast-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x bowtie_ref -U $selected_probe -S probes.sam # samp9 -> 67575 lines

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
cmd="mv $outdir/output_1.fastq.gz $outdir/output_ref1_1.fastq.gz"; run "$cmd"
cmd="mv $outdir/output_2.fastq.gz $outdir/output_ref1_2.fastq.gz"; run "$cmd"

# run for 2nd haplotype
cmd="$JAPSA_CAPSIM --reference $mod_fasta2 $capsim_options --logFile $outdir/capsim2.log "; run "$cmd"
cmd="mv $outdir/output_1.fastq.gz $outdir/output_ref2_1.fastq.gz"; run "$cmd"
cmd="mv $outdir/output_2.fastq.gz $outdir/output_ref2_2.fastq.gz"; run "$cmd"

# reset java
echo "Reset Java to original version"
export PATH="$ORIGINAL_PATH"


# merge results
# fq1
cmd="cat $outdir/output_ref1_1.fastq.gz $outdir/output_ref2_1.fastq.gz"
echo "$cmd > $outdir/output_1.fastq.gz"
$cmd > $outdir/output_1.fastq.gz

# fq2
cmd="cat $outdir/output_ref1_2.fastq.gz $outdir/output_ref2_2.fastq.gz"
echo "$cmd > $outdir/output_2.fastq.gz"
$cmd > $outdir/output_2.fastq.gz

echo "done"






