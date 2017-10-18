#!/bin/bash

function 1-generateCNV()
{
    #mkdir -p $cnvdir
    mkdir -p $D_SCRATCHDIR/gencnv
    gencnvdir=$D_SCRATCHDIR/gencnv
    rocdir=$TEST_RUNNER_ROOT/bin/structural_vars/ROC
    cnvsimultoolsdir=$TEST_RUNNER_ROOT/bin/structural_vars/simulation
    addDNAdbdir=$TEST_RUNNER_ROOT/utility_scripts
    faSomeRecDir=/home/sameerg/UCSC_tools

    GLOBAL_dataset_name="RSVSim_"$chrom"_"$test_case_count

    mkdir -p /mnt/archive/TEST_RUNNER/$invocation_id/$GLOBAL_dataset_name
    #Generate truth CSVs
    Rscript $cnvsimultoolsdir/RSVSim_generate_cnv.R $chrom $gencnvdir $sizeins $sizedel $sizedup

    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

    #Combine truth CSVs into a truth vcf
    cat $gencnvdir/*.csv | $cnvsimultoolsdir/convert_RSVSim_to_VCF.sh -i /dev/stdin -f $D_VC_REF -c $cnvsimultoolsdir/hg19.contig_names.txt -v $cnvsimultoolsdir/hg19_header.vcf > $gencnvdir/cnvtruth.vcf

    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

    perl $rocdir/rsvsim.vcf_to_cnv_truth_table.pl $gencnvdir/cnvtruth.vcf > $gencnvdir/$GLOBAL_dataset_name.truth.cnv.tsv

    cd $gencnvdir
    bgzip -c $gencnvdir/cnvtruth.vcf > $gencnvdir/cnvtruth.vcf.gz
    tabix $gencnvdir/cnvtruth.vcf.gz

    echo "$chrom" >> $gencnvdir/$chrom.bed

    $faSomeRecDir/faSomeRecords $D_VC_REF $gencnvdir/$chrom.bed $gencnvdir/$chrom".fa"

    bcftools consensus -c $gencnvdir/liftover_1.txt -H 1 -f $gencnvdir/$chrom".fa" $gencnvdir/cnvtruth.vcf.gz 1> $gencnvdir/bcffasta_1.fa

    bcftools consensus -c $gencnvdir/liftover_2.txt -H 2 -f $gencnvdir/$chrom".fa" $gencnvdir/cnvtruth.vcf.gz 1> $gencnvdir/bcffasta_2.fa

    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
 
    
    pirs simulate -l 100 -x 35 --insert-len-mean=400 --insert-len-sd=40 --diploid \
        --base-calling-profile=/opt/pirs-2.0.1/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz \
        --indel-error-profile=/opt/pirs-2.0.1/Profiles/InDel_Profiles/phixv2.InDel.matrix \
        --gc-bias-profile=/opt/pirs-2.0.1/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat \
        --phred-offset=33 --no-indel-errors --no-gc-bias -c "gzip" -t 24 $gencnvdir/bcffasta_1.fa $gencnvdir/bcffasta_2.fa 

    rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
    
    rm -f $gencnvdir/*.fa

    
    cp $gencnvdir/* /mnt/archive/TEST_RUNNER/$invocation_id/$GLOBAL_dataset_name
    
    
    echo $GLOBAL_dataset_name >> /mnt/archive/TEST_RUNNER/$invocation_id/input.txt
    echo $dg_ref >> /mnt/archive/TEST_RUNNER/$invocation_id/input.txt
    cd /mnt/archive/TEST_RUNNER/$invocation_id/$GLOBAL_dataset_name
    ls $PWD/*fq.gz >> /mnt/archive/TEST_RUNNER/$invocation_id/input.txt
    ls -d $PWD/$GLOBAL_dataset_name.truth.cnv.tsv
    echo /mnt/archive/TEST_RUNNER/$invocation_id/$GLOBAL_dataset_name/$GLOBAL_dataset_name.truth.cnv.tsv >> /mnt/archive/TEST_RUNNER/$invocation_id/input.txt


}
