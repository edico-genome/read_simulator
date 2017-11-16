#!/usr/bin/env bash

##################################################
# REF
# FASTA=/mnt/vault/reference_genomes/Hsapiens/hs37d5/seq/hs37d5.fa
# FASTA=/home/gavinp/Downloads/new_bamsurgeon/bamsurgeon/hs37d5.fa


##################################################
# Constants
NUM_PROCESSES=4


##################################################
# CLI
usage=$(cat <<EOF
usage: $(basename $0) 
  -f			          fasta
  -w                              workdir
  -o                              outdir
  -b                              tumor bam name
  -n                              normal bam
  -e                              bamsurgeon bed
EOF
)

while getopts ':f:w:o:b:' opt ; do
    case "$opt" in
    f) FASTA="${OPTARG}" ;;
    w) WORKDIR="${OPTARG}" ;;	
    b) TUMOR_BAM_NAME="${OPTARG}" ;;	
    o) OUTDIR="${OPTARG}"
       if [[ ! -d $OUTDIR  ]]; then
	   echo "Outdir: $OUTDIR must exist"
	   exit 1
       fi ;;
    n) NORMAL_BAM="${OPTARG}" ;;
    e) BAMSUREON_BED="${OPTARG}" ;;
    :) echo "$usage"
       exit 1 ;;    
    *) echo "$usage"
       exit 1 ;;
  esac
done


INBAM=$NORMAL_BAM
NOSECONDARY=$WORKDIR/no_secondary_aligns_test.bam
SORTED_BAM=$WORKDIR/reads_sorted.bam
RG_UPDATED_BAM=$WORKDIR/new_rgids.bam
TUMOR_BAM=$WORKDIR/$TUMOR_BAM_NAME


echo ""
echo "SETTINGS"
echo "=============================="
for key in FASTA WORKDIR OUTDIR NORMAL_BAM INBAM NOSECONDARY SORTED_BAM RG_UPDATED_BAM TUMOR_BAM BAMSUREON_BED; do
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


##################################################
abort()
{
if [[ $? -ne 0 ]]; then 
    echo "An error occurred. Exiting..." >&2
    exit 1
else  
    echo "Complete"
    exit 0
fi 
}

trap 'abort' 0
set -e


##################################################
filter_secondary_alignments_from_BAM() 
{
    echo "Filter secondary alignments from BAM"
    echo "samtools view -F 3840 $INBAM -b > $NOSECONDARY"
    samtools view -F 3840 $INBAM -b > $NOSECONDARY
    echo "Created: $NOSECONDARY"
}

##################################################
sort_bam() 
{
    _bam=$1
    _bam_tmp=$WORKDIR/_tmp.bam

    echo "SORT BAM"
    cmd="picard SortSam INPUT=$_bam OUTPUT=$_bam_tmp SORT_ORDER=coordinate CREATE_INDEX=true "
    echo $cmd
    $cmd
    echo "mv $_bam_tmp $_bam"
    mv $_bam_tmp $_bam
}

##################################################
update_rgids()
{
    echo "UPDATE RGIDS"

    cmd="picard AddOrReplaceReadGroups INPUT=$NOSECONDARY" 
    cmd+=" OUTPUT=$RG_UPDATED_BAM " 
    cmd+=" RGID=tumor "
    cmd+=" RGLB=library1 "
    cmd+=" RGPL=illumina " 
    cmd+=" RGPU=H0164ALXX140820.2 "
    cmd+=" RGSM=tumor_sample " 
    cmd+=" SORT_ORDER=coordinate " 
    cmd+=" CREATE_INDEX=true "
    cmd+=" TMP_DIR=$WORKDIR "
    
    echo $cmd
    $cmd
    echo "CREATED $RG_UPDATED_BAM"

    echo tabix $RG_UPDATED_BAM
    tabix $RG_UPDATED_BAM
}


##################################################
run_bamsurgeon()
{
    # activate python venv
    cmd="source /home/theoh/git_master/read_simulator/venv/bin/activate"
    echo "Activate python venv"
    echo $cmd
    $cmd
    
    echo RUN BAMSURGEON
    echo index bam 
    samtools index $RG_UPDATED_BAM

    script="/home/theoh/git/bamsurgeon/bin/addsnv.py"
    cmd="python $script -f $RG_UPDATED_BAM -r $FASTA -o $TUMOR_BAM -v $BAMSUREON_BED --mindepth 15 -p $NUM_PROCESSES"
    echo "$cmd > $OUTDIR/bamsurgeon.log"
    $cmd > $OUTDIR/bamsurgeon.log
}

##################################################
move_to_nas()
{
    echo "move BAM to nas"
    echo "cp $TUMOR_BAM\* $OUTDIR"
    cp ${TUMOR_BAM}* $OUTDIR
}


##################################################
# main 

filter_secondary_alignments_from_BAM
sort_bam $NOSECONDARY
update_rgids
run_bamsurgeon
sort_bam $TUMOR_BAM 
move_to_nas
