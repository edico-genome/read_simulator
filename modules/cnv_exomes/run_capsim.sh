#!/usr/bin/env bash

# exit the scipt if any statement returns a non-true value
set -e

#####################################################
# constants

JAPSA_CAPSIM=/home/gavinp/.usr/local/bin/jsa.sim.capsim

# list of processed probes files to select from
PROBE_FILES=(\
  "/home/gavinp/Downloads/sequencing_kits/S03723314/S03723314_Probes_chr20.txt.fa", \
  "/home/gavinp/Downloads/sequencing_kits/SeqCapEZ/results3.txt", \
  "/home/gavinp/Downloads/sequencing_kits/S07604624/results.txt")

#seed random generator
RANDOM=$$$(date +%s)

#####################################################
# Read the options this was called with
usage=$(cat <<EOF
usage: $(basename $0) 
  -e                              probe entry (directory)
  -n			          number of reads
  -l                              read length
  -f                              fragment length
  -o                              outdir
  -1                              modified fasta1
  -2                              modified fasta2
EOF
)

while getopts ':n:l:f:e:' opt ; do
    case "$opt" in
    1) mod_fasta1="${OPTARG}" ;;
    2) mod_fasta2="${OPTARG}" ;;	
    n) N_READS="${OPTARG}" ;;
    l) READ_LENGTH="${OPTARG}" ;;
    f) FRAG_LENGTH="${OPTARG}" ;;
    e) entry="D_${OPTARG}"
       if [[ ! -d $entry  ]]; then
	   echo "Entry: $entry must be a directory"
	   # exit 1
       fi
       ;;
    o) outdir="${OPTARG}" ;;
    :) echo "$usage"
       exit 1 ;;    
    *) echo "$usage"
       exit 1 ;;
  esac
done


echo ""
echo "SETTINGS"
echo "=============================="
for key in entry N_READS READ_LENGTH FRAG_LENGTH outdir; do
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

echo "DONE"
exit 0

#####################################################
function run {
    cmd=$1
    echo $cmd
    $cmd
}

#####################################################
# MAIN

echo "Working on $entry"

# avoid nextera - it is not working properly in the 
# simulator - seems to be a capsim issue, maybe related to length
selected_probe=${PROBE_FILES[$RANDOM % ${#PROBE_FILES[@]} ]}

echo "Processing $entry $selected_probe "
SAMP_DIR=$PWD/$entry

pushd $entry
echo $entry $selected_probe
echo "pwd = $PWD"

cmd="bowtie2-build $mod_fasta1 $outdir/bowtie_ref1"; run "$cmd"
cmd="bowtie2-build $mod_fasta2 $outdir/bowtie_ref2"; run "$cmd"

# add -p NTHREADS for multiple threads
# samp9->164851 lines
# samp9->164851 lines
cmd="bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x bowtie_ref2 -U $selected_probe -S $outdir/probes.sam"; run "$cmd"

#bowtie2 --local --fast-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x bowtie_ref -U $selected_probe -S probes.sam # samp9 -> 67575 lines

cmd="/home/gavinp/bin/bin/samtools view -b $outdir/probes.sam | /home/gavinp/bin/bin/samtools sort -o $outdir/probes.bam"; run "$cmd"

cmd="samtools index $outdir/probes.bam"; run "$cmd"

capsim_options=" --probe probes.bam --ID someid --fmedian $FRAG_LENGTH --miseq output --illen $READ_LENGTH --num $N_READS"

# run for 1st haplotype
cmd="$JAPSA_CAPSIM --reference genome_rearranged1.fasta $capsim_options --logFile $outdir/capsim1.log"
cmd="mv $outdir/output_1.fastq.gz $outdir/output_ref1_1.fastq.gz"; run "$cmd"
cmd="mv $outdir/output_2.fastq.gz $outdir/output_ref1_2.fastq.gz"; run "$cmd"

# run for 2nd haplotype
cmd="$JAPSA_CAPSIM --reference genome_rearranged2.fasta $capsim_options --logFile $outdir/capsim2.log"
cmd="mv $outdir/output_1.fastq.gz $outdir/output_ref2_1.fastq.gz"; run "$cmd"
cmd="mv $outdir/output_2.fastq.gz $outdir/output_ref2_2.fastq.gz"; run "$cmd"

# merge results
# fq1
cmd="cat $outdir/output_ref1_1.fastq.gz $outdir/output_ref2_1.fastq.gz"
echo "$cmd > $outdir/output_1.fastq.gz"
$cmd > $outdir/output_1.fastq.gz
# fq2
cmd="cat $outdir/output_ref1_2.fastq.gz $outdir/output_ref2_2.fastq.gz"
echo "$cmd > $outdir/output_2.fastq.gz"
$cmd > $outdir/output_2.fastq.gz

echo "Finished $entry"
