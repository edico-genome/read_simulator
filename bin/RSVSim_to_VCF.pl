#!/usr/bin/perl

# Input variant format will be the output from order_variants.sh. 
# Will also require the associated FASTA as input.

# Will use samtools faidx input.fasta chr5:238972-338754 to fetch
# reference sequence for insertion into REF or ALT field of VCF record.

use strict; use warnings;

#use File::Cat;

sub usage {

   print "Usage:\n  RSVSim_to_VCF.pl <unified_SV_table.csv> <ref.fasta>\n";
   return;
}

if (@ARGV != 2) {
   usage();
   exit(1);
}

my $infile = shift @ARGV;
my $fasta = shift @ARGV;

#$tmpfile=`mktemp /tmp/simvar.XXXXXX`;

my $prob_hom_del = 0.1;

my $QUAL = 999;
my $FILTER="PASS";
my $INFO = "DP=20";
my $FORMAT = "GT";
my $GENO = "0/1";
# theo hack
# my $GENO = "1/1";

open (INFILE, $infile) or die $!;
while (my $line = <INFILE>) {
   chomp($line);
   my @field = split /\t/, $line;
   next if ($field[0] =~ /^#/);

   my $chrom = $field[0];
   my $pos = $field[1]-1;   # for struct vars, VCF spec wants the POS of the base before the variant.
   my $var_id = $field[3];

   die "Found a START POS = 1!" if ($pos == 0);

   my $seq_range = "$chrom:$pos-$field[2]";
   my $seq = uc `samtools faidx $fasta $seq_range | tail -n +2 | tr -d '\n'`;
   if (length($seq)-1 != $field[7]) {
      die "Error detected related to variant length: $line";
   }
   my $pos_base = substr($seq, 0 , 1);

   my $refseq = "";  my $altseq = "";

   my $info = $INFO;

   # genotype for current variant
   my $gtype = "0/1";
   $gtype = "1/0" if (rand(1) > 0.5);
   # theo hack
   # $gtype = "1/1";
   
   if ($var_id =~ /^deletion/) {

       $refseq = $seq;
       $altseq = $pos_base;
       $gtype = "1/1" if (rand(1) <= $prob_hom_del);

   } elsif ($var_id =~ /^insertion/) {

       # Inserted SEQ is defined by 1st 3 input fields and position of insertion
       # is defined by fields 4-6 (0-based field numbering). A CNV detection tool
       # attempts to detect the source pos and length of the inserted SEQ (fields 0-2).
       # We include this truth info in the INFO field of insertion VCF records.
       my $src_chrom = $chrom;
       my $src_pos = $field[1];
       $info .= ";SRCINS=$src_chrom:$src_pos-$field[2]";

       # The position of the insertion is reported in the CHROM and POS vcf fields.
       $chrom = $field[4];
       $pos = $field[5]-1;
       die "Found a START POS = 1! Have not dealt with this corner case yet." if ($pos == 0);
       $pos_base = uc `samtools faidx $fasta $chrom:$pos-$pos | tail -n +2 | tr -d '\n'`;
       $refseq = $pos_base;
       $altseq = $pos_base . substr($seq, 1); 

   } elsif ($var_id =~ /^tandemDuplication/) {

       my $ndups = $field[8];
       my $dupseq = substr($seq, 1);  # $seq includes the prior base
       $pos += length($dupseq);       # set $pos to POS of the base at end of dupseq.
       $pos_base = substr($dupseq, -1);  # ref base at POS position
       $refseq = $pos_base;
       $altseq =  $pos_base . $dupseq x $ndups; # pos_base + dup segs
       $info .= ";TDUP=$ndups";
   } else {
     warn "Variant ID type is invalid: $var_id ...skipping.";
     next;
   }

   print "$chrom\t$pos\t$var_id\t$refseq\t$altseq\t$QUAL\t$FILTER\t$info\t$FORMAT\t$gtype\n"
}
close INFILE;

# Sort the VCF:  sortBed -faidx /home/cooper/p4/sw/trunk/test/bin/structural_vars/simulation/hg19.contig_names.txt -i tester.vcf > tester.sorted.vcf
# And add the VCF header:  cat vcf_hdr.txt tester.sorted.vcf > good.vcf

