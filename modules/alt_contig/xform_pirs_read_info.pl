#!/usr/bin/perl

use strict; use warnings;

use File::Basename qw(dirname);
my $exec_dir;
BEGIN { $exec_dir = dirname($0); }
use lib "$exec_dir";
use 5.010;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use AltAware::Coords qw(cigar_pos_map);

sub usage {
   print STDERR "Usage:\n";
   print STDERR "\t xform_read_info.pl --read-info pirs_simulation.read.info --alt-sam alt.sam --contig1 chr6:28510120-33383765 --contig1-rename chr6 --contig1-offset 28510120 --contig2 chr6_GL000251v2_alt:1-4795265 --contig2-rename chr6_GL000251v2_alt --contig2-offset 1 > read.info.xform\n";
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

my ($readinfoFile, $altSam);
my ($contig1, $contig1Rename, $contig1Offset);
my ($contig2, $contig2Rename, $contig2Offset);
GetOptions(
    'read-info=s' => \$readinfoFile,
    'alt-sam=s' => \$altSam,
    'contig1=s' => \$contig1,
    'contig1-rename=s' => \$contig1Rename,
    'contig1-offset=s' => \$contig1Offset,
    'contig2=s' => \$contig2,
    'contig2-rename=s' => \$contig2Rename,
    'contig2-offset=s' => \$contig2Offset,
) or usage;

my @offsetOpt = ();
my $readLen;
my %altToPriXform = ();

if (!defined($readinfoFile) or !defined($altSam)) {
   #print STDERR "readinfoFile=$readinfoFile    altSam=$altSam\n";
   print STDERR "Invalid command line.\n";
   usage;
}

push (@offsetOpt, "$contig1,$contig1Rename,$contig1Offset") if (defined $contig1);
push (@offsetOpt, "$contig2,$contig2Rename,$contig2Offset") if (defined $contig2);

# Process the mapping offset optional arg if present
my %map_offset = ();
for my $item (@offsetOpt) {
      my ($seqname, $conv_name, $sq_offset) = split ",", $item;
      $map_offset{$seqname}{NEWNAME} = $conv_name;
      $map_offset{$seqname}{OFFSET} = $sq_offset-1;   

      # Check in altSam if this is an ALT contig
      open my $samFH, "$altSam" or die $!;
      while(<$samFH>)
      {
          next if /^@/;
          my @t = split;
          next if ($t[0] ne $conv_name or $t[1] > 200);
          my $alt_cigar = $t[5];

          # We have an ALT contig, so we need primary contig name, mapping POS, FLAG
          # and CIGAR for it.

          $altToPriXform{$conv_name}{'FLAG'} = $t[1];
          $altToPriXform{$conv_name}{'RNAME'} = $t[2];
          $altToPriXform{$conv_name}{'POS'} = $t[3];
          $altToPriXform{$conv_name}{'CMAP'} = cigar_pos_map($alt_cigar);  # Returns CMAP array ref.
                                                                          # CMAP maps ALT pos to a
                                                                          # CIGAR-adjusted pos.
          $altToPriXform{$conv_name}{'LENGTH'} = scalar(@{$altToPriXform{$conv_name}{'CMAP'}});

          last;
      }
      close $samFH;
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
printf STDERR "Transforming truth alignments\n";
my $numReads = 0;
my %cigPosMap = (); 
while(my $line = <$infoFH>) {
  chomp $line;

  my ($qname, $ref, $seq, $pos, $str, $insLen, $junk, $sub, $ins, $del) = split(/\t/,$line);
  # $pos is 0-based mapping position

  $seq =~ s/ .*//;  # Discard the first space char in the pirs contig name and all chars after.
  $pos += 0;

  if (defined $map_offset{$seq}) {
     $pos += $map_offset{$seq}{OFFSET};  # $pos is still 0-based after this
     $seq = $map_offset{$seq}{NEWNAME};  # can't use $map_offset for current read after this step!
  }

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
  if (defined($altToPriXform{$seq})) {
     # This read was simulated from an ALT contig, so get corresponding
     # mapping POS on primary contig and populate our ZM tag with it.
     #print STDOUT "Before\t$seq:$pos\n";

     my $pri_contig_flag = $flag;
     my $pri_contig_name = $altToPriXform{$seq}{'RNAME'};

     my $delta_p;
     if ($altToPriXform{$seq}{'FLAG'} & 0x10) {
         # This ALT contig is RC-oriented w.r.t. its primary assembly contig,
         # so our ALT $pos has an opposite effect on the primary contig mapping pos.
         my $rc_pos = $altToPriXform{$seq}{'LENGTH'} - $pos - $readLen;
         $delta_p = $altToPriXform{$seq}{'CMAP'}->[$rc_pos];

         # The map orientation bits in primary contig FLAG must be switched such
         # that 0x10 becomes 0x20 and vice versa.
         $pri_contig_flag += ($flag & 0x10) ? 0x10 : -0x10;
     } else {
         # This ALT contig is aligned FWD w.r.t. its primary assembly contig.
         $delta_p = $altToPriXform{$seq}{'CMAP'}->[$pos];
     }

     my $pri_contig_pos = $altToPriXform{$seq}{'POS'} + $delta_p;
     if ($pri_contig_pos <= 0) {
        die "Invalid primary contig mapping tag for $seq:$pos -> $pri_contig_name";
     }

     #print STDOUT "ZM:Z:$pri_contig_name;$pri_contig_pos;$pri_contig_flag\n";
     $pri_contig_tag = "ZM:Z:$pri_contig_name;$pri_contig_pos;$pri_contig_flag";

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

  printf STDERR "\r$numReads" if ++$numReads % 100000 == 0;
}
close $infoFH;
printf STDERR "\rTotal reads processed = $numReads\n";

