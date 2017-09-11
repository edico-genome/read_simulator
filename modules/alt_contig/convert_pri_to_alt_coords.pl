#!/usr/bin/perl -w

use strict;
use warnings;

use v5.10.1;

# This tool converts user-specified primary contig start and end coords to the
# corresponding coords of the specified overlapping ALT contig.

my ($samFN, $altContigName, $startPos, $endPos) = @ARGV;

sub usage {

    printf STDERR "Usage:\n\tconvert_pri_to_alt_coords.pl <alt.contig.sam> <alt.contig.name> <priContigStartPos> <priContigEndPos>\n";
    exit(1); 
}

if (!defined($startPos) or !defined($endPos) or $startPos < 1 or $endPos < $startPos) {
   printf STDERR "startPos must be greater than zero and endPos cannot be less than startPos.\n\n";

   usage;
}

sub cigar_array {

    my @result = ();
    my ($cigar) = @_;
    my @op_arr = $cigar =~ /[\d]+[MDISHN=X]/g;

    for my $op (@op_arr) {
        my @len = ($op =~ /^[\d]+/g);
        my @symbol = ($op =~ /[MDISHN=X]$/g);
        push (@result,["@symbol","@len"]);
    }
    return @result;

}

sub pos_from_cigar {

   my $altSeqStart;
   my $altSeqEnd;
   my $altSeqPos = 0;
   my $cigar = $_[0];
   my $priSeqPos = $_[1] - 1;
   my $priSeqStart = $_[2];
   my $priSeqEnd = $_[3];

   # print "$cigar\t$priSeqPos\t$priSeqStart\t$priSeqEnd\n";

   if ($priSeqStart <= $priSeqPos) {
      die("Start of specified range on primary seq must be >= ALT contig alignment POS.");
   }

   my @cigar_arr = cigar_array($cigar);

   foreach my $cig_op (@cigar_arr) {

      my ($op, $len) = @$cig_op;

      given ($op) {
          when ($_ eq "M") {
              $priSeqPos += $len;
              $altSeqPos += $len;
          }
          when ($_ eq "I") {
              $altSeqPos += $len;
          }
          when ($_ eq "D") {
              $priSeqPos += $len;
          }
          when ($_ eq "S" or $_ eq "H") {
              $altSeqPos += $len;
          }
      }

      if (not defined $altSeqStart and $priSeqPos >= $priSeqStart) {
          $altSeqStart = $altSeqPos;
          if ($op eq "M") {
              $altSeqStart -= ($priSeqPos-$priSeqStart);
          } else {
              die "Specified primary contig start position is not within a 'M' CIGAR op."
          }
      }
      if (not defined $altSeqEnd and $priSeqPos >= $priSeqEnd) {
          $altSeqEnd = $altSeqPos;
          if ($op eq "M") {
              $altSeqEnd -= ($priSeqPos-$priSeqEnd);
          } else {
              die "Specified primary contig stop position is not within a 'M' CIGAR op."
          }
          last;
      }
   }

   if ($altSeqStart < 1 or $altSeqEnd < 1) {
      printf STDERR "Could not determine altSeqStart or altSeqEnd position.\n";
      exit(1);
   }

   return ($altSeqStart,$altSeqEnd);
}

open my $samFH, "$samFN" or die $!;
while(<$samFH>)
{
   next if /^@/;
   my @t = split;
   next if ($t[0] ne $altContigName);

   my $priSeqName = $t[2];
   my $priSeqPos1 = $t[3];
   my $cigar = $t[5];

   my ($altSeqStartPos, $altSeqEndPos) = pos_from_cigar($cigar, $priSeqPos1, $startPos, $endPos);

   printf STDOUT "$priSeqName\t$startPos\t$endPos\t$altContigName\t$altSeqStartPos\t$altSeqEndPos\n";

   last;
}
close $samFH;

