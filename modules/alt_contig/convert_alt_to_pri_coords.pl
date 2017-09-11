#!/usr/bin/perl -w

use strict;
use warnings;

use v5.10.1;

# This tool converts user-specified ALT contig start and end coords to the
# corresponding coords of the overlapping primary contig. There may be
# unexpected results when startPos or endPos are not located with 'M' ops.

my ($samFN, $altContigName, $startPos, $endPos) = @ARGV;

sub usage {

    print STDERR "Usage:\n\tconvert_alt_to_pri_coords.pl alt.contig.sam alt.contig.name [altContigStartPos] [altContigEndPos]\n";
    exit(1); 
}

if (!defined($samFN) or !defined($altContigName)) {

   print STDERR "ALT SAM file and a contig name must be specified.\n";
   usage;
}

if (defined($startPos) and $startPos < 1) {

   print STDERR "startPos must be positive.\n\n";
   usage;
}

if (defined($endPos) and $endPos < $startPos) {

   print STDERR "endPos is less than startPos.\n";
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

   my $priSeqStart;
   my $priSeqEnd;
   my $cigar = $_[0];
   my $priSeqPos1 = $_[1];
   my $altSeqStart = $_[2];
   my $altSeqEnd = $_[3];

   # print "$cigar\t$priSeqPos1\t$altSeqStart\t$altSeqEnd\n";

   my @cigar_arr = cigar_array($cigar);

   # The following two variables are always sync'ed. Init them
   # to the zeroeth base of the ALT contig since we have not
   # yet started walking through the CIGAR.
   my $priSeqPos = $priSeqPos1 - 1;
   my $currAltPos = 0;

   foreach my $cig_op (@cigar_arr) {

      my ($op, $len) = @$cig_op;

      given ($op) {
          when ($_ eq "M") {
              $priSeqPos += $len;
              $currAltPos += $len;
          }
          when ($_ eq "I") {
              $currAltPos += $len;
          }
          when ($_ eq "D") {
              $priSeqPos += $len;
          }
          when ($_ eq "S") {
              $currAltPos += $len;
          }
      }

      if (not defined $priSeqStart and $currAltPos >= $altSeqStart) {
          $priSeqStart = $priSeqPos;
          if ($op eq "M") {
              $priSeqStart -= ($currAltPos-$altSeqStart);
          }
          if ($op eq "I" or $op eq "S") {
              # The following sets $priSeqStart to the primary contig
              # position sync'ed to the initial base of the next 'M'
              # or 'D' op. 
              # This works if the next op is a 'M' or 'D'; otherwise,
              # may cause problems. 

              $priSeqStart += 1;
          }
          if (!defined($altSeqEnd)) {
              $priSeqEnd = 0;
              last;
          }
          #print STDERR "$altContigName:  M op at altSeqStart is $len$op\n";
      }
      if (not defined $priSeqEnd and $currAltPos >= $altSeqEnd) {
          # If an 'I' or 'S' op caused us to be here, then the following will
          # set $priSeqEnd to the primary contig position sync'ed to the
          # last base of the most recent 'M' or 'D' op.

          $priSeqEnd = $priSeqPos;
          if ($op eq "M") {
              $priSeqEnd -= ($currAltPos-$altSeqEnd);
          }
          last;
      }
   }

   if (!defined($priSeqStart) or !defined($priSeqEnd) ) {
      print STDERR "Could not determine priSeqStart or priSeqEnd position.\n";
      exit(1);
   }

   return ($priSeqStart,$priSeqEnd);
}

open my $samFH, "$samFN" or die $!;
while(<$samFH>)
{
   next if /^@/;
   my @t = split;
   next if ($t[0] ne $altContigName or $t[1] > 200);

   my $priSeqName = $t[2];
   my $priSeqPos1 = $t[3];
   my $cigar = $t[5];

   if (!defined($startPos)) {
       print STDOUT "$priSeqName\t0\t0\t$altContigName\t0\t0\n";

   } elsif (!defined($endPos)) {
       my ($priSeqStartPos, $priSeqEndPos) = pos_from_cigar($cigar, $priSeqPos1, $startPos);
       print STDOUT "$priSeqName\t$priSeqStartPos\t0\t$altContigName\t$startPos\t0\n";

   } else {
       my ($priSeqStartPos, $priSeqEndPos) = pos_from_cigar($cigar, $priSeqPos1, $startPos, $endPos);
       print STDOUT "$priSeqName\t$priSeqStartPos\t$priSeqEndPos\t$altContigName\t$startPos\t$endPos\n";
   }
   last;
}
close $samFH;

