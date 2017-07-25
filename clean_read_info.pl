#!/usr/bin/perl

# Usage:
#   clean_read_info.pl < pirs_simulation.read.info > clean.read.info

# This utility condenses a Pirs truth alignment "read.info" file in two ways:
#   1. The reference FASTA name "/staging/tmp/theoh/simReads/sra[12].fa" is replaced with just the number [12]
#   2. Substitution entries with unchanged nucleotide, usually changing case such as 'a'->'A', are deleted

while(my $line = <>) {
    ((print $line), next) if $line =~ /^\W/;
    chomp $line;
    my @fields = split(/\t/,$line);
  # Replace "/staging/tmp/theoh/simReads/sra2.fa" with "2"
  $fields[1] =~ s/^.*(\d)\D*$/$1/;
  # Delete no-change substitutions such as 'a'->'A'
    my $subs = uc $fields[7];
  $subs =~ s/\d+,(A->A|C->C|G->G|T->T);//g;
  $subs = "-" if $subs eq "";
    $fields[7] = $subs;
    print join("\t",@fields) . "\n";
}
