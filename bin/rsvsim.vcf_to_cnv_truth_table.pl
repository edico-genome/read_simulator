#!/usr/bin/perl

use strict; use warnings;

sub usage {

   print "Usage:\n  rsvsim.vcf_to_cnv_truth_table.pl <rsvsim.vcf> \n";
   print "        Extracts CNV's from RSVSim VCF and writes to truth table.\n";
   return;
}

if ( scalar(@ARGV) != 1 ) {
   usage();
   exit 1;
}

my ($vcfFN) = @ARGV;

printf STDOUT "#CHROM\tPOS\tTYPE\tLENGTH\tSEGMENT_VALUE\tCOPY_NUMBER_FLAG\n";
open my $vcfFH, "$vcfFN" or die $!;
while(<$vcfFH>)
{
      next if /^#/;
      chomp $_;
      my ($chrom,$pos,$id,$ref,$alt,$qual,$filt,$info,$format,$geno) = split(/\t/,$_);

      $pos += 1;  # POS in VCF is for the base before the start of a CNV

      my $cnv_length = length($alt) - length($ref);
      $cnv_length = -$cnv_length if ($cnv_length < 0);

      my $type = "CNV";
      $type = "DEL" if ( $id =~ /^deletion/ );
      $type = "TANDUP" if ( $id =~ /^tandemDuplication/ );
      $type = "INS" if ( $id =~ /^insertion/ );

      my $alt_count = "NA";
      $alt_count = 1 if ($geno eq "0/1" or $geno eq "1/0");
      $alt_count = 2 if ($geno eq "1/1");

      if ($alt_count eq "NA") {
          printf "Unexpected genotype found: $geno ...halting\n";
          exit 1;
      }

      my $allele_count = -1;
      $allele_count = 2-$alt_count if ($type eq "DEL"); 
      $allele_count = 2+$alt_count if ($type eq "INS");

      if ($type eq "TANDUP") {
          my @fields = split /;/, $info;
          for my $fld (@fields) {
              my ($tag,$val) = split /=/, $fld;
               $allele_count = 2+$val if ($tag eq "TDUP");

          }
          $cnv_length /= ($allele_count-2);   # we want the length of the ref seq that was copied $val times
          $pos -= $cnv_length;  # For a TANDUP vcf record, the POS was moved to end of first dup'ed segment.
      }
      if ($type eq "INS") {
          my @fields = split /;/, $info;
          for my $fld (@fields) {
             my ($tag,$val) = split /=/, $fld;
             if ($tag eq "SRCINS") {
                 # Parse string containing reference location of insertion seq source
                 my ($cntg, $istart, $iend) = split /[:-]/, $val;
                 $chrom = $cntg;
                 $pos = $istart;
             }
         }
      }

      
      if ($allele_count < 0 || $allele_count == 2) {
          printf "Unexpected variant record detected:\n   $_\n";
          exit 1;
      }
      my $cnv_flag = ".";
      $cnv_flag = "+" if ($allele_count > 2);
      $cnv_flag = "-" if ($allele_count < 2);
      if ($cnv_flag eq ".") {
          printf "Invalid CNV_FLAG detected ...halting\n";
          exit 1; 
      }
      printf STDOUT "$chrom\t$pos\t$type\t$cnv_length\t$allele_count\t$cnv_flag\n";

}
close $vcfFH;
