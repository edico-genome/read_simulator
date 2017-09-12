#!/usr/bin/perl

use strict; use warnings;

use File::Basename qw(dirname);
my $exec_dir;
BEGIN { $exec_dir = dirname($0); }
use lib "$exec_dir";
use AltAware::Coords qw(cigar_pos_map);

sub usage {
   print STDERR "Usage:\n";
   print STDERR "\t xform_read_info.pl pirs_simulation.read.info alt.sam [chr6_s::chr6::29775250] > read.info.xform\n";
   print STDERR "\t    (See inline documentation for more details.)\n";
   exit(1);
}

# Given the "truth" alignments in *.read.info from a Pirs simulation
# involving ALT contigs and/or sub-seq's of primary assembly contigs, add
# primary assembly mapping coord as SAM tag to any ALT alignments. Also add
# a user-specified offset value to the truth mapping POS of all reads 
# simulated from shortened primary contigs. Also, the alignment RNAME field
# for reads simulated from a shortened primary contig can be changed to
# the primary contig name. The format of the 3rd (optional) arg on the cmd
# line is "shortened_RNAME:pri_assembly_RNAME:pri_assembly_POS_offset". This
# POS offset value corresponds to the mapping POS on the original contig for
# the first base in the shortened contig.

# Output format, tab-separated:
#   QNAME  FLAG  RNAME  POS  254  *  =  PNEXT  TLEN  *  *  indel.tag  pri.contig.tag
# 
#   QNAME = read name, same as in *.read.info
#   FLAG = sum of segment and R-C FLAG bits
#   RNAME = hg19 sequence name
#   POS = 1-based mapping position for this read
#   PNEXT = mapping position of mate
#   TLEN = insert length, with sign equal to Pirs strand, as in SAM
#   indel.tag (yi) = 1 if this mate has a possible indel relative to ref, else 0
#   pri.contig.tag = mapping pos on primary contig, if read simulated from ALT contig
# 
# Note: CIGAR string not emitted.

# Key notes:
#  The read length is hard coded below
#

my $readLen;

my ($readinfoFile, $altSam, $offsetOpt) = @ARGV;

if (scalar(@ARGV) < 2 or scalar(@ARGV) > 3) {
   print STDERR "Invalid command line.\n";
   usage;
}

# Process the POS offset optional arg if present
my $mod_sqname;
my $pri_sqname;
my $sq_offset;
if (defined $offsetOpt) {
   ($mod_sqname, $pri_sqname, $sq_offset) = split "::", $offsetOpt;
}

# This section parses the header of the read.info file
my $infoFH;
open $infoFH, "$readinfoFile" or die $!;
while(my $line = <$infoFH>) {
  chomp $line;
  if ($line =~ /^#/) {
     last if (substr($line,0,8) eq "# readId");
     $line =~ s/ //g;
     my @t = split ":", $line;
     if (scalar(@t) == 2 and $t[0] eq "#Readlength") {
         $readLen = $t[1];
     }
  } else {
     print "$readLen    $line\n";
     die "Unexpectedly read past the end of the header section in $readinfoFile.";
  }
}
die "Read length setting not found in $readinfoFile." if (not defined $readLen);

# This section parses the Pirs "*.read.info" truth alignments and converts them to SAM
printf "Transforming truth alignments\n";
my $numReads = 0;
my %altToPriXform = ();
my %cigPosMap = (); 
while(my $line = <$infoFH>) {
  chomp $line;

  my ($qname, $ref, $seq, $pos, $str, $insLen, $junk, $sub, $ins, $del) = split(/\t/,$line);
  # $pos is 0-based mapping position

  $seq =~ s/ .*//;  # Discard the first space char in the pirs contig name and all chars after.

  my $pri_contig_tag = "";

  # Detect if any indels are listed
  my $indel = ($ins ne '-' || $del ne '-') ? 1 : 0;

  # Adjust the position if this is the (-) strand mate. Pirs records the position of the left
  # ('+') mate for both mates.  For the ('-') mate, add the insert length minus the read length.
  $pos += $insLen - $readLen if $str eq "-";

  # Put strand flag into TLEN-ish field
  my $tlen = ($str eq "-") ? -$insLen : $insLen;

  # Initialize a 1-based alignment position
  my $alnPos = $pos + 1;

  my $flag = 1;  # paired-end reads
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

  # For reads sim'ed from an ALT contig, add primary assembly mapping position
  # as a SAM tag. If a user-specified position-offset arg was provided on cmd
  # line, check if this $seq matches the seq.name in the optional arg. If it
  # does, then add the offset to pirs' reported mapping position.
  if ($seq =~ /^HLA/ or $seq =~ /_alt$/) {
     # This read was simulated from an ALT contig, so get corresponding
     # mapping POS on primary contig.
     #print STDOUT "Before\t$seq:$pos\n";

     if (not defined $altToPriXform{$seq}) {

        # Only have to do the following _ONCE_ for each ALT contig encountered.

        # We need primary contig name, mapping POS, and CIGAR from ALT SAM file
        # for this ALT contig $seq.
        open my $samFH, "$altSam" or die $!;
        while(<$samFH>)
        {
           next if /^@/;
           my @t = split;
           next if ($t[0] ne $seq or $t[1] > 200);
           my $alt_cigar = $t[5];

           $altToPriXform{$seq}{'FLAG'} = $t[1];
           $altToPriXform{$seq}{'RNAME'} = $t[2];
           $altToPriXform{$seq}{'POS'} = $t[3];
           $altToPriXform{$seq}{'CMAP'} = cigar_pos_map($alt_cigar);  # Returns CMAP array ref.
                                                                      # CMAP maps ALT pos to a
                                                                      # CIGAR-adjusted pos.
           $altToPriXform{$seq}{'LENGTH'} = scalar(@{$altToPriXform{$seq}{'CMAP'}});

           last;
        }
        close $samFH;
     }

     # At this point it's an error if our table is not defined
     if (not defined $altToPriXform{$seq}) {
        die "Unable to create alt_to_pri_map for $seq.";
     }

     my $pri_contig_flag = $flag;
     my $pri_contig_name = $altToPriXform{$seq}{'RNAME'};

     # CMAP->[$pos] maps ALT $pos to a alt->pri CIGAR ajusted offset
     my $delta_p = $altToPriXform{$seq}{'CMAP'}->[$pos];

     if ($altToPriXform{$seq}{'FLAG'} & 0x10) {
         # This ALT contig is RC-oriented w.r.t. its primary assembly contig,
         # so our ALT $pos has an opposite effect on the primary contig mapping pos.
         $delta_p = $altToPriXform{$seq}{'LENGTH'} - $delta_p - $readLen;

         # The map orientation bits in primary contig FLAG must be switched such
         # that 0x10 becomes 0x20 and vice versa.
         $pri_contig_flag += ($flag & 0x10) ? 0x10 : -0x10;
     }
     my $pri_contig_pos = $altToPriXform{$seq}{'POS'} + $delta_p;
     if ($pri_contig_pos <= 0) {
        die "Invalid primary contig mapping tag for $seq:$pos -> $pri_contig_name";
     }

     #print STDOUT "ZM:Z:$pri_contig_name;$pri_contig_pos;$pri_contig_flag\n";
     $pri_contig_tag = "ZM:Z:$pri_contig_name;$pri_contig_pos;$pri_contig_flag";

  } elsif (defined $mod_sqname and ($seq eq $mod_sqname) ) {
     # This read was simulated from a short section of a primary contig, so
     # use the offset from our optional cmd line arg to adjust the mapping POS.
     $alnPos += ($sq_offset - 1);
     $seq = $pri_sqname;
  }

  # Output
  # Strip off trailing /1 or /2 in qname, if present
  $qname =~ s{/[12]$}{};  
  #print "$qname\t$flag\t$seq\t$alnPos\t$tlen\t$indel\n";
  
  my $dLen = $readLen;
  $dLen = -$dLen if ( $tlen > 0 );
     
  my $pnext = $alnPos + $tlen + $dLen;
  my $indel_tag="yi:i:$indel";
  my $tags = $indel_tag;
  $tags .= "\t$pri_contig_tag" if ($pri_contig_tag ne "");

  print "$qname\t$flag\t$seq\t$alnPos\t254\t\*\t=\t$pnext\t$tlen\t\*\t\*\t$tags\n";

  print "\r$numReads" if ++$numReads % 100000 == 0;
}
close $infoFH;
print "\rTotal reads processed = $numReads\n";
exit(0);

