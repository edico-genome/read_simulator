set -e  # exit the scipt if any statement returns a non-true value
OUTPUT_LOG=./sample_log_file.txt

#seed random generator
RANDOM=$$$(date +%s)
#clear log file
:> $OUTPUT_LOG

PARALLEL_EXE=/home/gavinp/bin/parallel
PARALLEL_JOBS=8

makecapture(){
    entry=${1}
    echo "Working on $entry"
    if [[ -f $entry  ]]; then
        echo "Skipping $entry since it is a file"
        return 1
    fi
    JAPSA_CAPSIM=/home/gavinp/.usr/local/bin/jsa.sim.capsim
    #N_READS=400000  # N = C.G/L = 50 x 1e6 / 100 for chr20
    N_READS=800000  # N = C.G/L = 50 x 1e6 / 100 for chr20
                # capsim says N_READS is 'Number of Fragments'
    READ_LENGTH=300 
    FRAG_LENGTH=500
    # list of processed probes files to select from
    PROBE_FILES=("/home/gavinp/Downloads/sequencing_kits/S03723314/S03723314_Probes_chr20.txt.fa" "/home/gavinp/Downloads/sequencing_kits/SeqCapEZ/results3.txt", "/home/gavinp/Downloads/sequencing_kits/S07604624/results.txt")
    # avoid nextera - it is not working properly in the 
    # simulator - seems to be a capsim issue, maybe related to length
    selected_probe=${PROBE_FILES[$RANDOM % ${#PROBE_FILES[@]} ]}
    OUTPUT_LOG=./sample_log_file.txt
    echo $entry $selected_probe >> $OUTPUT_LOG
    echo "Processing $entry $selected_probe "
    SAMP_DIR=$PWD/$entry
    REF=$SAMP_DIR/genome_rearranged.fasta
    PROBE_FA_FILE=$selected_probe
    pushd $entry
    echo $entry $selected_probe > $OUTPUT_LOG
    echo "pwd = $PWD"
    bowtie2-build genome_rearranged1.fasta bowtie_ref1
    bowtie2-build genome_rearranged2.fasta bowtie_ref2
    # add -p NTHREADS for multiple threads
    bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x bowtie_ref2 -U $selected_probe -S probes.sam # samp9->164851 lines
    bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x bowtie_ref2 -U $selected_probe -S probes.sam # samp9->164851 lines
    #bowtie2 --local --fast-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x bowtie_ref -U $selected_probe -S probes.sam # samp9 -> 67575 lines
    /home/gavinp/bin/bin/samtools view -b probes.sam | /home/gavinp/bin/bin/samtools sort -o probes.bam 
    samtools index probes.bam
    $JAPSA_CAPSIM --reference genome_rearranged1.fasta --probe probes.bam --ID someid --fmedian $FRAG_LENGTH --miseq output --illen $READ_LENGTH --num $N_READS --logFile capsim.log
    mv output_1.fastq.gz output_ref1_1.fastq.gz
    mv output_2.fastq.gz output_ref1_2.fastq.gz
    $JAPSA_CAPSIM --reference genome_rearranged2.fasta --probe probes.bam --ID someid --fmedian $FRAG_LENGTH --miseq output --illen $READ_LENGTH --num $N_READS --logFile capsim.log
    mv output_1.fastq.gz output_ref2_1.fastq.gz
    mv output_2.fastq.gz output_ref2_2.fastq.gz
    cat output_ref1_1.fastq.gz output_ref2_1.fastq.gz > output_1.fastq.gz
    cat output_ref1_2.fastq.gz output_ref2_2.fastq.gz > output_2.fastq.gz
    echo "Finished $entry"
}
export -f makecapture

shopt -s nullglob
list=(./*samp*[0-9])
echo $list
echo "Starting $PARALLEL_JOBS parallel jobs"
$PARALLEL_EXE --line-buffer -j $PARALLEL_JOBS makecapture ::: "${list[@]}"
