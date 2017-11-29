#!/usr/bin/perl

# Usage:
#   xform_read_info.pl pirs_simulation.read.info liftover_basename > read.info.xform

# Given the "truth" alignments in *.read.info from a Pirs simulation,
# based on input references sra[12].fa which were transformed from hg19
# by bcftools using a VCF such as from GiaB, this script transforms the
# truth alignments to hg19 coordinates, and does some reformatting.

# Output format, tab-separated:
#   QNAME  FLAG  SEQ  POS  TLEN  INDEL
# 
#   QNAME = read name, same as in *.read.info
#   FLAG = sum of segment and R-C FLAG bits
#   SEQ = hg19 sequence name
#   POS = hg19 sequence offset, 1-based, for this mate
#   TLEN = Pirs insert length, with sign equal to Pirs strand, as in SAM
#   INDEL = 1 if this mate has a possible indel relative to hg19, else 0
# 
# Note the lists of substitutions and indels in *.read.info do not appear
# in the output.  This is mainly because the "liftover" files are ambiguous
# about exact alignments around block boundaries.  A downstream tool should
# re-align the reads to hg19 at their transformed positions to determine
# detailed alignments if needed.  Simple base comparison or gapless alignment
# may be used when INDEL=0; when INDEL=1, Smith-Waterman should be used, with
# the reference window extended some adequate distance each way, such as 50bp.

# Key notes:
#  1. The read length is hard coded below
#  2. Liftover files with the names below must be present

my ($readinfoFile, $liftoverBasename) = @ARGV;

my $readLen = 100;
my @liftoverNames = (join("_", $liftoverBasename, "1.liftover.txt"), 
                     join("_", $liftoverBasename, "2.liftover.txt") );

# This section parses each "liftover" file, and saves the information in a data structure.
# The data structure is as folows:
#   $xfmBlocks[$s]{$seq}[$n] = ($refPos, $qryPos, $blockLen, $refExtra, $qryExtra)
#     * $s is the strand, 0 or 1, corresponding to the [12] in sra[12].fa and sra[12]_chain_liftover.txt
#     * $seq is the reference sequence name, e.g. "chr1"
#     * $n is an incrementing array index
#     * $refPos is the block start position in hg19
#     * $qryPos is the block start position in sra[12].fa
#     * $blockLen is the length perfectly aligned starting at $refPos and $qryPos
#     * $refExtra is the following length of the unaligned portion in hg19
#     * $qryExtra is the following length of the unaligned portion in sra[12].fa
# Normally, except for final blocks, at least one of $refExtra and $qryExtra is nonzero.
# When both $refExtra and $qryExtra are nonzero, the relative alignment of these two segments is undefined.
for my $s (0,1) {
  open LIFT, "<", $liftoverNames[$s] or die "Can't open '$liftoverNames[$s]'";
  print STDERR "Loading '$liftoverNames[$s]'\n";
  my $seq = undef;
  my $refPos = 0;
  my $qryPos = 0;
  my $ap = undef;
  my $src = undef;
  my $nLines = 0;
  while(my $line = <LIFT>) {
    $nLines++;
    if($line =~ /^chain \d+ (\w+)/) {
      if(defined $ap) {
        my @cols = ($refPos, $qryPos, 0, 0, 0);
        push @{$ap}, \@cols;
      }
      $seq = $1;
      $refPos = 0;
      $qryPos = 0;
      my @a = ();
      $ap = \@a;
      $xfmBlocks[$s]{$seq} = $ap;
      next;
    }
    next unless $line =~ /^\d/;
    my ($blockLen, $refExtra, $qryExtra) = split(/\s/,$line);
    my @cols = ($refPos, $qryPos, $blockLen, $refExtra+0, $qryExtra+0);
    push @{$ap}, \@cols;
    $refPos += $blockLen + $refExtra;
    $qryPos += $blockLen + $qryExtra;
  }
  my @cols = ($refPos, $qryPos, 0, 0, 0);
  push @{$ap}, \@cols;
  close LIFT;
}

# This section parses the Pirs "*.read.info" truth alignments, and transforms to hg19 coordinates.
printf STDERR "Transforming truth alignments\n";
my $numReads = 0;
my $infoFH;
open $infoFH, "$readinfoFile" or die $!;
while(my $line = <$infoFH>) {
  next if $line =~ /^\W/;
  chomp $line;
  my ($qname, $ref, $seq, $pos, $str, $insLen, $junk, $sub, $ins, $del) = split(/\t/,$line);
  # Map "/staging/tmp/theoh/simReads/sra[12].fa" to 0 or 1
  $ref =~ /([12])\D*$/ or die "Read info reference string should contain '1' or '2' indicating source FASTA";
  $ref = $1 - 1;
  # Detect if any indels are listed
  my $indel = ($ins ne '-' || $del ne '-') ? 1 : 0;
  # Adjust the position for this mate.  Pirs records the position of the left ('+') mate
  # for both mates.  For the right ('-') mate, add the insert length minus the read length.
  $pos += $insLen - $readLen if $str eq "-";
  # Put strand flag into TLEN-ish field
  my $tlen = ($str eq "-") ? -$insLen : $insLen;
  # Initialize a 1-based alignment position
  my $alnPos = $pos + 1;
  # Grab the appropriate transform block array
  my $bp = $xfmBlocks[$ref]{$seq};
  if(defined $bp) {
    # Perform a binary search for the last transform block with $qryPos <= $pos,
    # which also implies $pos < $qryPos+$blockLen+$qryExtra
    my $nb = scalar(@$bp);
    my $b = int($nb / 2);
    my $d = int($nb / 4);
    while($d) {
      $b += ($bp->[$b][1] > $pos) ? -$d : $d;
      $d = int($d/2);
    }
    $b++ while $bp->[$b][1] < $pos;
    $b-- while $bp->[$b][1] > $pos;
    my ($refPos, $qryPos, $blockLen, $refExtra, $qryExtra) = @{$bp->[$b]};
    # Measure offset into this block or its following query-extra zone
    my $intoBlock = $pos - $qryPos;
    my $intoIns = ($intoBlock > $blockLen) ? $intoBlock - $blockLen : 0;
    $indel = 1 if $intoBlock + $readLen > $blockLen;
    !$intoIns || $intoIns < $qryExtra;
    # Scale query-extra offset linearly to an approximate reference-extra offset
    my $intoDel = $intoIns ? int($intoIns * $refExtra / $qryExtra) : 0;
    $intoBlock -= $intoIns;
    # Calculate transformed position in hg19 sequence, 1-based
    $alnPos = $refPos + $intoBlock + $intoDel + 1;
  }

  my $flag = 1; 
  if (substr($qname,-2) eq '/1') {
     $flag += 0x40;
  } elsif (substr($qname,-2) eq '/2') {
     $flag += 0x80;
  }
  if ($str eq "-") {
     $flag += 0x10;
  } else {
     $flag += 0x20;
  }

  # Output
  # Strip off trailing /1 or /2 in qname, if present
  $qname =~ s{/[12]$}{};  
  #print "$qname\t$flag\t$seq\t$alnPos\t$tlen\t$indel\n";
  
  my $dLen = $readLen;
  $dLen = -$dLen if ( $tlen > 0 );

  my $pnext = $alnPos + $tlen + $dLen;
  my $indel_tag="yi:i:$indel";
  print "$qname\t$flag\t$seq\t$alnPos\t254\t$readLen","M\t=\t$pnext\t$tlen\t\*\t\*\t$indel_tag\n";

  printf STDERR "\r$numReads" if ++$numReads % 100000 == 0;
}
close $infoFH;
printf STDERR "\r$numReads\n";

